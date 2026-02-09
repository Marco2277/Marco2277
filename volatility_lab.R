rm(list = ls())

suppressPackageStartupMessages({
  library(xts); library(qrmdata); library(qrmtools)
  library(rugarch); library(PerformanceAnalytics); library(quantmod); library(zoo)
  library(ggplot2); library(dplyr); library(tidyr); library(MCS)
  library(nortest); library(ADGofTest); library(TTR); library(vars)
})

CONFIG <- list(
  seed = 123,
  end_date = as.Date("2026-01-19"),
  regime = list(start = as.Date("2025-06-01"), end = as.Date("2026-01-19")),
  A = list(tickers = c("GLD", "SPY", "TLT"), oos_target = 30L, oos_min = 0.10, oos_max = 0.30, refit_every = 10),
  B = list(start = as.Date("2006-01-03"), end = as.Date("2025-12-15"), ticker = "GLD", oos_frac = 0.25, refit_every = 50),
  MCS = list(alpha = 0.10, B = 1000),
  loss = list(min_sigma2 = 1e-10),
  output_dir = "output",
  verbose = FALSE
)

set.seed(CONFIG$seed)
options(getSymbols.warning4.0 = FALSE)

output_figures <- file.path(CONFIG$output_dir, "figures")
output_tables <- file.path(CONFIG$output_dir, "tables")
dir.create(output_figures, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tables, showWarnings = FALSE, recursive = TRUE)

# ------------------------------ Funzioni base -----------------------
ewma_sigma2 <- function(r, lambda = 0.94) {
  r2 <- as.numeric(r)^2
  s2 <- numeric(length(r2))
  s2[1] <- r2[1]
  for (i in 2:length(r2)) s2[i] <- lambda * s2[i - 1] + (1 - lambda) * r2[i - 1]
  xts(s2, order.by = index(r))
}

GARCH_nll <- function(params, r) {
  omega <- params[1]; alpha <- params[2]; beta <- params[3]
  r2 <- as.numeric(r)^2
  Tn <- length(r2)
  c0 <- omega / (1 - alpha - beta)
  sigma2 <- matrix(data = c0, nrow = Tn, ncol = 1)
  logL <- matrix(0, nrow = Tn, ncol = 1)
  for (i in 2:Tn) {
    sigma2[i] <- omega + alpha * r2[i - 1] + beta * sigma2[i - 1]
    logL[i] <- -0.5 * log(2 * pi) - 0.5 * log(sigma2[i]) - 0.5 * (r2[i] / sigma2[i])
  }
  -sum(logL)
}

roll_sigma2 <- function(spec, series, oos_frac = 0.25, refit_every = 10, refit_window = "moving") {
  N <- NROW(series)
  oos <- max(1, round(oos_frac * N))
  rr <- ugarchroll(spec = spec, data = series, n.ahead = 1, forecast.length = oos,
                   refit.every = refit_every, refit.window = refit_window,
                   solver = "hybrid", keep.coef = TRUE)
  df <- as.data.frame(rr)
  idx <- if ("Index" %in% names(df)) {
    as.Date(df$Index)
  } else if ("Date" %in% names(df)) {
    as.Date(df$Date)
  } else {
    as.Date(rownames(df))
  }
  xts((df$Sigma)^2, order.by = idx)
}

qlike_from <- function(rv, s2) {
  eps <- .Machine$double.eps
  rv <- pmax(as.numeric(rv), eps); s2 <- pmax(as.numeric(s2), eps)
  rv / s2 - log(rv / s2) - 1
}

enforce_sigma2_floor <- function(sigma2_xts, floor_value, model_name, pipeline_label) {
  sigma2_xts <- xts::as.xts(sigma2_xts)
  sigma2_raw <- as.numeric(sigma2_xts)
  sigma2_clean <- sigma2_raw
  sigma2_clean[!is.finite(sigma2_clean)] <- NA_real_
  sigma2_used <- ifelse(is.na(sigma2_clean), NA_real_, pmax(sigma2_clean, floor_value))
  floored <- is.finite(sigma2_clean) & sigma2_clean < floor_value
  dates <- as.character(index(sigma2_xts))

  summary_row <- data.frame(
    pipeline = pipeline_label,
    model = model_name,
    n_oos = length(sigma2_raw),
    floor_value_used = floor_value,
    n_floored = sum(floored, na.rm = TRUE),
    pct_floored = mean(floored, na.rm = TRUE),
    min_before = min(sigma2_clean, na.rm = TRUE),
    min_after = min(sigma2_used, na.rm = TRUE),
    q01_before = quantile(sigma2_clean, 0.01, na.rm = TRUE, names = FALSE),
    q01_after = quantile(sigma2_used, 0.01, na.rm = TRUE, names = FALSE),
    stringsAsFactors = FALSE
  )

  if (!any(floored, na.rm = TRUE)) {
    floor_dates <- data.frame(
      pipeline = character(0),
      model = character(0),
      date = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    floor_dates <- data.frame(
      pipeline = pipeline_label,
      model = model_name,
      date = dates[floored],
      stringsAsFactors = FALSE
    )
  }

  list(
    sigma2 = xts(sigma2_used, order.by = index(sigma2_xts)),
    summary_row = summary_row,
    floor_dates = floor_dates
  )
}

run_mcs <- function(loss_matrix, alpha, B) {
  mcs_obj <- tryCatch(
    MCS::MCSprocedure(Loss = loss_matrix, alpha = alpha, B = B, statistic = "Tmax"),
    error = function(e) NULL
  )
  if (is.null(mcs_obj)) {
    list(obj = MCS::MCSprocedure(Loss = loss_matrix, alpha = alpha, B = B, statistic = "TR"),
         statistic = "TR")
  } else {
    list(obj = mcs_obj, statistic = "Tmax")
  }
}

get_mcs_models <- function(mcs_obj) {
  out <- tryCatch(MCS::getMCS(mcs_obj), error = function(e) NULL)
  if (!is.null(out) && !is.null(out$models)) return(out$models)
  if (isS4(mcs_obj)) {
    slots <- methods::slotNames(mcs_obj)
    if ("MCS" %in% slots) {
      mcs_tbl <- methods::slot(mcs_obj, "MCS")
      if (!is.null(mcs_tbl)) return(rownames(mcs_tbl))
    }
    if ("SuperiorModel" %in% slots) {
      sup <- methods::slot(mcs_obj, "SuperiorModel")
      if (!is.null(sup)) return(as.character(sup))
    }
    return(character(0))
  }
  if (!is.null(mcs_obj$models)) return(mcs_obj$models)
  if (!is.null(mcs_obj$MCS)) return(rownames(mcs_obj$MCS))
  character(0)
}

make_sigma_panel <- function(s2_list) {
  s2_list2 <- lapply(s2_list, function(x) {
    x <- xts::as.xts(x)
    storage.mode(x) <- "double"
    colnames(x) <- NULL
    x
  })
  v <- Reduce(function(a, b) merge(a, b, all = FALSE), s2_list2)
  colnames(v) <- names(s2_list)
  df <- data.frame(Date = zoo::index(v), zoo::coredata(v), check.names = FALSE)
  num_cols <- setdiff(names(df), "Date")
  for (nm in num_cols) {
    vals <- as.numeric(df[[nm]])
    vals[!is.finite(vals)] <- NA_real_
    vals[vals <= 0] <- NA_real_
    df[[nm]] <- vals
  }
  df_long <- tidyr::pivot_longer(df, -Date,
                                 names_to = "Model",
                                 values_to = "sigma2",
                                 values_transform = list(sigma2 = as.double))
  ggplot(df_long, aes(x = Date, y = sigma2, color = Model)) +
    geom_line(linewidth = 0.5) +
    scale_y_log10() +
    labs(title = "Conditional Variance (daily, log scale)",
         y = expression(hat(sigma)^2), x = NULL, color = NULL) +
    theme_minimal(12)
}

plot_nic <- function(nic_list, title) {
  df <- bind_rows(lapply(names(nic_list), function(nm) {
    data.frame(z = nic_list[[nm]]$zx, nic = nic_list[[nm]]$zy, Model = nm)
  }))
  ggplot(df, aes(x = z, y = nic, color = Model)) +
    geom_line(linewidth = 0.7) +
    labs(title = title, x = "z", y = "News Impact") +
    theme_minimal(12)
}

summarize_mcs <- function(mcs_obj, statistic, alpha, B, path) {
  models <- get_mcs_models(mcs_obj)
  summary_df <- data.frame(
    statistic = statistic,
    alpha = alpha,
    B = B,
    superior_set = paste(models, collapse = ", "),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, path, row.names = FALSE)
  summary_df
}

prepare_loss_matrix <- function(S2, rv, pipeline_label, floor_value) {
  rv_vec <- as.numeric(rv)
  sigma2_floor <- floor_value
  oos_dates <- as.character(index(S2))
  if (length(oos_dates) == 0) {
    oos_dates <- as.character(seq_along(rv_vec))
  }

  Loss_cols <- lapply(seq_len(ncol(S2)), function(j) {
    sigma2_raw <- as.numeric(S2[, j])
    sigma2_used <- pmax(sigma2_raw, sigma2_floor)
    qlike_from(rv_vec, sigma2_used)
  })
  Loss <- do.call(cbind, Loss_cols)
  colnames(Loss) <- colnames(S2)

  floor_summary <- lapply(seq_len(ncol(S2)), function(j) {
    model_name <- if (!is.null(colnames(S2))) colnames(S2)[j] else paste0("Model_", j)
    sigma2_raw <- as.numeric(S2[, j])
    sigma2_used <- pmax(sigma2_raw, sigma2_floor)
    floored <- sigma2_raw < sigma2_floor
    dates_floored <- oos_dates[floored]
    data.frame(
      pipeline = pipeline_label,
      model = model_name,
      n_oos = length(sigma2_raw),
      floor_value_used = sigma2_floor,
      n_floored = sum(floored, na.rm = TRUE),
      pct_floored = mean(floored, na.rm = TRUE),
      first_floored_date = if (length(dates_floored) > 0) dates_floored[1] else NA_character_,
      last_floored_date = if (length(dates_floored) > 0) dates_floored[length(dates_floored)] else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  floor_summary <- dplyr::bind_rows(floor_summary)

  floor_dates <- lapply(seq_len(ncol(S2)), function(j) {
    model_name <- if (!is.null(colnames(S2))) colnames(S2)[j] else paste0("Model_", j)
    sigma2_raw <- as.numeric(S2[, j])
    sigma2_used <- pmax(sigma2_raw, sigma2_floor)
    floored <- sigma2_raw < sigma2_floor
    if (!any(floored, na.rm = TRUE)) {
      return(data.frame(
        pipeline = character(0),
        model = character(0),
        date = character(0),
        sigma2_raw = numeric(0),
        sigma2_used = numeric(0),
        rv = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      pipeline = pipeline_label,
      model = model_name,
      date = oos_dates[floored],
      sigma2_raw = sigma2_raw[floored],
      sigma2_used = sigma2_used[floored],
      rv = rv_vec[floored],
      stringsAsFactors = FALSE
    )
  })
  floor_dates <- dplyr::bind_rows(floor_dates)

  rows_ok <- apply(Loss, 1, function(r) all(is.finite(r)))
  Loss <- Loss[rows_ok, , drop = FALSE]
  cols_ok <- apply(Loss, 2, function(c) all(is.finite(c)) && sd(c) > 0)
  Loss <- as.matrix(Loss[, cols_ok, drop = FALSE])
  if (ncol(Loss) < 2L) stop(sprintf("MCS %s: meno di 2 modelli validi dopo pulizia.", pipeline_label))
  dup_cols <- which(duplicated(as.data.frame(round(Loss, 12))))
  if (length(dup_cols)) Loss <- Loss[, -dup_cols, drop = FALSE]
  cols_ok2 <- apply(Loss, 2, function(c) sd(c) > 1e-10)
  Loss <- as.matrix(Loss[, cols_ok2, drop = FALSE])
  if (ncol(Loss) < 2L) stop(sprintf("MCS %s: meno di 2 modelli validi dopo pulizia.", pipeline_label))
  if (nrow(Loss) < 15) stop(sprintf("MCS %s: OOS troppo corto (<15).", pipeline_label))

  list(loss = Loss, floor_summary = floor_summary, floor_dates = floor_dates, floor_value = sigma2_floor)
}

save_qlike_tables <- function(loss_matrix, pipeline_label, mcs_models, path) {
  means <- apply(loss_matrix, 2, mean, na.rm = TRUE)
  sds <- apply(loss_matrix, 2, sd, na.rm = TRUE)
  n_oos <- nrow(loss_matrix)
  df <- data.frame(
    model = names(means),
    mean_qlike = as.numeric(means),
    sd_qlike = as.numeric(sds),
    n_oos = n_oos,
    included_in_MCS = names(means) %in% mcs_models,
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(path, paste0("qlike_means_", pipeline_label, ".csv")), row.names = FALSE)
  df
}

save_qlike_boxplot <- function(loss_matrix, title, path, use_log_y = FALSE, winsor_p = NULL) {
  df <- as.data.frame(loss_matrix)
  df_long <- tidyr::pivot_longer(df, cols = everything(), names_to = "Model", values_to = "QLIKE")
  if (!is.null(winsor_p)) {
    cap <- quantile(df_long$QLIKE, winsor_p, na.rm = TRUE)
    df_long$QLIKE <- pmin(df_long$QLIKE, cap)
  }
  if (use_log_y) {
    eps <- .Machine$double.eps
    df_long$QLIKE <- ifelse(df_long$QLIKE <= 0, eps, df_long$QLIKE)
  }
  p <- ggplot(df_long, aes(x = Model, y = QLIKE, fill = Model)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal(12) +
    theme(legend.position = "none") +
    labs(title = title, x = NULL, y = "QLIKE")
  if (use_log_y) {
    p <- p + ggplot2::scale_y_log10()
  }
  ggsave(path, p, width = 9, height = 5, dpi = 150)
  p
}

# ------------------------------ Download dati -----------------------
startA <- CONFIG$regime$start
endA <- CONFIG$regime$end
getSymbols(CONFIG$A$tickers, from = startA, to = CONFIG$end_date, auto.assign = TRUE)
getSymbols(CONFIG$B$ticker, from = CONFIG$B$start, to = CONFIG$end_date, auto.assign = TRUE)
getSymbols(c("^GVZ", "^VIX", "^MOVE", "^OVX"),
           src = "yahoo",
           from = "2020-11-01", to = CONFIG$end_date,
           auto.assign = TRUE)
GVZ <- na.locf(GVZ, na.rm = FALSE)
VIX <- na.locf(VIX, na.rm = FALSE)
MOVE <- na.locf(MOVE, na.rm = FALSE)
OVX <- na.locf(OVX, na.rm = FALSE)

# ------------------------------ Pipeline A (regime mid-2025 → 2025-12-15) -----------------------
GLDpx <- Ad(GLD); SPYpx <- Ad(SPY); TLTpx <- Ad(TLT)
GLDpx_A <- GLDpx[paste0(startA, "/", endA)]
ret_A <- diff(log(GLDpx_A))[-1, ]

sigma2_floor_value <- max(CONFIG$loss$min_sigma2, .Machine$double.eps)
sigma2_RM_A_raw <- ewma_sigma2(ret_A, 0.94)
rm_floor_A <- enforce_sigma2_floor(sigma2_RM_A_raw, sigma2_floor_value, "RiskMetrics", "A")
sigma2_RM_A <- rm_floor_A$sigma2

x0 <- c(0.000003, 0.2, 0.6)
A_mat <- t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, -1), nrow = 3, ncol = 4))
b_vec <- c(rep(10^(-16), 3), -1)
resA <- constrOptim(x0, function(p) GARCH_nll(p, ret_A), NULL, ui = A_mat, ci = b_vec)
theta_garch_A <- setNames(resA$par, c("omega", "alpha", "beta"))

specs <- list(
  RiskMetrics = NA,
  sGARCH = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                      distribution.model = "norm"),
  eGARCH = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                      distribution.model = "norm"),
  apARCH = ugarchspec(variance.model = list(model = "apARCH", garchOrder = c(1, 1)),
                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                      distribution.model = "norm"),
  iGARCH = ugarchspec(variance.model = list(model = "iGARCH", garchOrder = c(1, 1)),
                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                      distribution.model = "norm"),
  fiGARCH = ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                       distribution.model = "norm")
)

fitsA <- lapply(specs[-1], function(sp) ugarchfit(data = ret_A, spec = sp, solver = "hybrid"))
sigma2_list_A <- list(
  RiskMetrics = sigma2_RM_A,
  sGARCH = sigma(fitsA$sGARCH)^2,
  eGARCH = sigma(fitsA$eGARCH)^2,
  apARCH = sigma(fitsA$apARCH)^2,
  iGARCH = sigma(fitsA$iGARCH)^2,
  fiGARCH = sigma(fitsA$fiGARCH)^2
)

z_grid <- seq(-5, 5, length.out = 501)
nic_list_A <- lapply(fitsA[c("sGARCH", "eGARCH", "apARCH")], function(f) newsimpact(f, z = z_grid))
roll_models_A <- c("sGARCH", "eGARCH", "apARCH")

OOS_TARGET <- CONFIG$A$oos_target
N_A <- NROW(ret_A)
oos_frac_A <- max(CONFIG$A$oos_min, min(CONFIG$A$oos_max, OOS_TARGET / N_A))

roll_A_raw <- lapply(specs[roll_models_A], function(sp) {
  roll_sigma2(sp, ret_A, oos_frac = oos_frac_A, refit_every = CONFIG$A$refit_every, refit_window = "moving")
})
floor_summaries_A <- list(rm_floor_A$summary_row)
floor_dates_A <- list(rm_floor_A$floor_dates)
roll_A <- lapply(seq_along(roll_A_raw), function(i) {
  model_name <- roll_models_A[i]
  floored <- enforce_sigma2_floor(roll_A_raw[[i]], sigma2_floor_value, model_name, "A")
  floor_summaries_A[[length(floor_summaries_A) + 1]] <<- floored$summary_row
  floor_dates_A[[length(floor_dates_A) + 1]] <<- floored$floor_dates
  floored$sigma2
})

cat("OOS fraction A =", round(oos_frac_A, 3), " (N =", N_A, ")\n")

sizes_A <- sapply(roll_A, NROW)
na_all_A <- sapply(roll_A, function(x) if (NROW(x) == 0) TRUE else all(!is.finite(coredata(x))))
keep_A <- (sizes_A > 0) & (!na_all_A)
if (!any(keep_A)) stop("Pipeline A: tutti i forecast OOS sono vuoti/NA.")
roll_A <- roll_A[keep_A]
models_kept_A <- roll_models_A[keep_A]
cat("Modelli tenuti nel rolling A:", paste(models_kept_A, collapse = ", "), "\n")
cat("NROW per modello A:", paste(paste(models_kept_A, sizes_A[keep_A], sep = ":"), collapse = " | "), "\n")
if (length(roll_A) < 2L) stop("Pipeline A: meno di 2 modelli con forecast OOS validi.")

S2_core_A <- Reduce(function(a, b) merge(a, b, join = "inner"), roll_A)
if (NROW(S2_core_A) == 0) stop("Pipeline A: nessuna data in comune tra forecast OOS.")
colnames(S2_core_A) <- models_kept_A

RM_A_oos <- sigma2_RM_A[index(S2_core_A)]
S2_A <- merge(S2_core_A, RiskMetrics = RM_A_oos, join = "inner")

RV_all_A <- xts(as.numeric(ret_A)^2, order.by = index(ret_A))
RV_shift_A <- xts::lag.xts(RV_all_A, k = -1)
colnames(RV_shift_A) <- "RV"
combo_A <- merge(S2_A, RV_shift_A, join = "inner")
combo_A <- na.omit(combo_A)
if (!"RV" %in% colnames(combo_A) || NROW(combo_A) == 0) {
  stop("Pipeline A: nessuna sovrapposizione tra forecast e RV.")
}

rv_A <- combo_A[, "RV", drop = FALSE]
S2_A <- combo_A[, setdiff(colnames(combo_A), "RV"), drop = FALSE]
cat("OOS righe A dopo allineamento =", NROW(S2_A),
    " | colonne modelli =", NCOL(S2_A), "\n")

loss_A_res <- prepare_loss_matrix(S2_A, rv_A, "A", sigma2_floor_value)
Loss_QLIKE_A <- loss_A_res$loss
mcs_A_res <- run_mcs(Loss_QLIKE_A, CONFIG$MCS$alpha, CONFIG$MCS$B)
mcs_A <- mcs_A_res$obj

oos_dates_A <- rownames(Loss_QLIKE_A)
if (is.null(oos_dates_A)) {
  oos_dates_A <- as.character(index(S2_A))
}
oos_dates_A <- as.character(oos_dates_A)
diag_rows <- lapply(seq_len(ncol(S2_A)), function(j) {
  model_name <- colnames(S2_A)[j]
  sigma2_vals <- as.numeric(S2_A[, j])
  loss_vals <- as.numeric(Loss_QLIKE_A[, model_name])
  ord <- order(loss_vals, decreasing = TRUE, na.last = NA)
  top_idx <- ord[seq_len(min(3, length(ord)))]
  top_dates <- rep(NA_character_, 3)
  top_vals <- rep(NA_real_, 3)
  if (length(top_idx) > 0) {
    top_dates[seq_along(top_idx)] <- oos_dates_A[top_idx]
    top_vals[seq_along(top_idx)] <- loss_vals[top_idx]
  }
  data.frame(
    model = model_name,
    min_sigma2 = min(sigma2_vals, na.rm = TRUE),
    q01_sigma2 = quantile(sigma2_vals, 0.01, na.rm = TRUE, names = FALSE),
    top1_date = top_dates[1],
    top1_qlike = top_vals[1],
    top2_date = top_dates[2],
    top2_qlike = top_vals[2],
    top3_date = top_dates[3],
    top3_qlike = top_vals[3],
    stringsAsFactors = FALSE
  )
})
pipelineA_outliers <- dplyr::bind_rows(diag_rows)
print(pipelineA_outliers)
write.csv(pipelineA_outliers,
          file.path(output_tables, "pipelineA_outlier_diagnostics.csv"),
          row.names = FALSE)
write.csv(loss_A_res$floor_summary,
          file.path(output_tables, "sigma2_floor_diagnostics_A.csv"),
          row.names = FALSE)
write.csv(loss_A_res$floor_dates,
          file.path(output_tables, "sigma2_floor_dates_A.csv"),
          row.names = FALSE)

lb_A <- bind_rows(lapply(names(fitsA), function(m) {
  z <- residuals(fitsA[[m]], standardize = TRUE)
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  z2 <- z^2
  lbz_p <- if (length(z) > 20) Box.test(z, lag = 20, type = "Ljung-Box")$p.value else NA_real_
  lbz2_p <- if (length(z2) > 20) Box.test(z2, lag = 20, type = "Ljung-Box")$p.value else NA_real_
  tibble(Model = m,
         LBz_p = lbz_p,
         LBz2_p = lbz2_p)
}))

# ------------------------------ Pipeline B (GLD long sample) -----------------------
GLDpx_B <- Ad(GLD)
GLDpx_B <- GLDpx_B[paste0(CONFIG$B$start, "/", CONFIG$B$end)]
ret_B <- diff(log(GLDpx_B))[-1]

sigma2_RM_B_raw <- ewma_sigma2(ret_B, 0.94)
rm_floor_B <- enforce_sigma2_floor(sigma2_RM_B_raw, sigma2_floor_value, "RiskMetrics", "B")
sigma2_RM_B <- rm_floor_B$sigma2

resB <- constrOptim(x0, function(p) GARCH_nll(p, ret_B), NULL, ui = A_mat, ci = b_vec)
theta_garch_B <- setNames(resB$par, c("omega", "alpha", "beta"))

fitsB <- lapply(specs[-1], function(sp) ugarchfit(data = ret_B, spec = sp, solver = "hybrid"))
sigma2_list_B <- list(
  RiskMetrics = sigma2_RM_B,
  sGARCH = sigma(fitsB$sGARCH)^2,
  eGARCH = sigma(fitsB$eGARCH)^2,
  apARCH = sigma(fitsB$apARCH)^2,
  iGARCH = sigma(fitsB$iGARCH)^2,
  fiGARCH = sigma(fitsB$fiGARCH)^2
)

nic_list_B <- lapply(fitsB[c("sGARCH", "eGARCH", "apARCH")], function(f) newsimpact(f, z = z_grid))

roll_models_B <- names(specs[-1])
roll_B_raw <- lapply(specs[roll_models_B], function(sp) {
  roll_sigma2(sp, ret_B, oos_frac = CONFIG$B$oos_frac, refit_every = CONFIG$B$refit_every, refit_window = "moving")
})
floor_summaries_B <- list(rm_floor_B$summary_row)
floor_dates_B <- list(rm_floor_B$floor_dates)
roll_B <- lapply(seq_along(roll_B_raw), function(i) {
  model_name <- roll_models_B[i]
  floored <- enforce_sigma2_floor(roll_B_raw[[i]], sigma2_floor_value, model_name, "B")
  floor_summaries_B[[length(floor_summaries_B) + 1]] <<- floored$summary_row
  floor_dates_B[[length(floor_dates_B) + 1]] <<- floored$floor_dates
  floored$sigma2
})

sizes_B <- sapply(roll_B, NROW)
na_all_B <- sapply(roll_B, function(x) if (NROW(x) == 0) TRUE else all(!is.finite(coredata(x))))
keep_B <- (sizes_B > 0) & (!na_all_B)
if (!any(keep_B)) stop("Pipeline B: tutti i forecast OOS sono vuoti/NA.")
roll_B <- roll_B[keep_B]
models_kept_B <- roll_models_B[keep_B]
cat("Modelli tenuti nel rolling B:", paste(models_kept_B, collapse = ", "), "\n")
cat("NROW per modello B:", paste(paste(models_kept_B, sizes_B[keep_B], sep = ":"), collapse = " | "), "\n")
if (length(roll_B) < 2L) stop("Pipeline B: meno di 2 modelli con forecast OOS validi.")

S2_core_B <- Reduce(function(a, b) merge(a, b, join = "inner"), roll_B)
if (NROW(S2_core_B) == 0) stop("Pipeline B: nessuna data in comune tra forecast OOS.")
colnames(S2_core_B) <- models_kept_B

RM_B_oos <- sigma2_RM_B[index(S2_core_B)]
S2_B <- merge(S2_core_B, RiskMetrics = RM_B_oos, join = "inner")

RV_all_B <- xts(as.numeric(ret_B)^2, order.by = index(ret_B))
RV_shift_B <- xts::lag.xts(RV_all_B, k = -1)
colnames(RV_shift_B) <- "RV"
combo_B <- merge(S2_B, RV_shift_B, join = "inner")
combo_B <- na.omit(combo_B)
if (!"RV" %in% colnames(combo_B) || NROW(combo_B) == 0) {
  stop("Pipeline B: nessuna sovrapposizione tra forecast e RV.")
}

rv_B <- combo_B[, "RV", drop = FALSE]
S2_B <- combo_B[, setdiff(colnames(combo_B), "RV"), drop = FALSE]
cat("OOS righe B dopo allineamento =", NROW(S2_B),
    " | colonne modelli =", NCOL(S2_B), "\n")

loss_B_res <- prepare_loss_matrix(S2_B, rv_B, "B", sigma2_floor_value)
Loss_QLIKE_B <- loss_B_res$loss
mcs_B_res <- run_mcs(Loss_QLIKE_B, CONFIG$MCS$alpha, CONFIG$MCS$B)
mcs_B <- mcs_B_res$obj

lb_B <- bind_rows(lapply(names(fitsB), function(m) {
  z <- residuals(fitsB[[m]], standardize = TRUE)
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  z2 <- z^2
  lbz_p <- if (length(z) > 20) Box.test(z, lag = 20, type = "Ljung-Box")$p.value else NA_real_
  lbz2_p <- if (length(z2) > 20) Box.test(z2, lag = 20, type = "Ljung-Box")$p.value else NA_real_
  tibble(Model = m,
         LBz_p = lbz_p,
         LBz2_p = lbz2_p)
}))
write.csv(loss_B_res$floor_summary,
          file.path(output_tables, "sigma2_floor_diagnostics_B.csv"),
          row.names = FALSE)
write.csv(loss_B_res$floor_dates,
          file.path(output_tables, "sigma2_floor_dates_B.csv"),
          row.names = FALSE)

floor_summary_source_A <- dplyr::bind_rows(floor_summaries_A)
floor_dates_source_A <- dplyr::bind_rows(floor_dates_A)
floor_summary_source_B <- dplyr::bind_rows(floor_summaries_B)
floor_dates_source_B <- dplyr::bind_rows(floor_dates_B)

# ------------------------------ Output principali (plots & tabelle) -----------------------
sigma2_list_A_plot <- sigma2_list_A
aparch_vals <- sigma2_list_A_plot$apARCH
exclude_aparch_A <- FALSE
if (!is.null(aparch_vals)) {
  aparch_vals_num <- as.numeric(aparch_vals)
  p_bad <- mean(!is.finite(aparch_vals_num) | aparch_vals_num <= 0, na.rm = TRUE)
  p_tiny <- mean(aparch_vals_num < 1e-12, na.rm = TRUE)
  if (is.finite(p_bad) && is.finite(p_tiny) && (p_bad > 0.05 || p_tiny > 0.50)) {
    sigma2_list_A_plot$apARCH <- NULL
    exclude_aparch_A <- TRUE
    message("Pipeline A: apARCH excluded from in-sample variance plot due to numerical degeneracy.")
  }
}
pA <- make_sigma_panel(sigma2_list_A_plot)
pB <- make_sigma_panel(sigma2_list_B)

ggsave(file.path(output_figures, "conditional_variance_A.png"), pA, width = 9, height = 5, dpi = 150)
ggsave(file.path(output_figures, "conditional_variance_B.png"), pB, width = 9, height = 5, dpi = 150)

if (exclude_aparch_A && "apARCH" %in% names(nic_list_A)) {
  nic_list_A <- nic_list_A[names(nic_list_A) != "apARCH"]
}
nic_plot_A <- plot_nic(nic_list_A, "News Impact Curve (A)")
nic_plot_B <- plot_nic(nic_list_B, "News Impact Curve (B)")
ggsave(file.path(output_figures, "nic_A.png"), nic_plot_A, width = 8, height = 5, dpi = 150)
ggsave(file.path(output_figures, "nic_B.png"), nic_plot_B, width = 8, height = 5, dpi = 150)

if (exclude_aparch_A && "apARCH" %in% lb_A$Model) {
  lb_A <- lb_A[lb_A$Model != "apARCH", , drop = FALSE]
  message("Pipeline A: apARCH excluded from Ljung-Box diagnostics due to numerical degeneracy.")
}
write.csv(lb_A, file.path(output_tables, "lb_A.csv"), row.names = FALSE)
write.csv(lb_B, file.path(output_tables, "lb_B.csv"), row.names = FALSE)

mcs_models_A <- get_mcs_models(mcs_A)
mcs_models_B <- get_mcs_models(mcs_B)

mcs_summary_A <- summarize_mcs(mcs_A, mcs_A_res$statistic, CONFIG$MCS$alpha, CONFIG$MCS$B,
                               file.path(output_tables, "mcs_summary_A.csv"))
mcs_summary_B <- summarize_mcs(mcs_B, mcs_B_res$statistic, CONFIG$MCS$alpha, CONFIG$MCS$B,
                               file.path(output_tables, "mcs_summary_B.csv"))

qlike_A_tbl <- save_qlike_tables(Loss_QLIKE_A, "A", mcs_models_A, output_tables)
qlike_B_tbl <- save_qlike_tables(Loss_QLIKE_B, "B", mcs_models_B, output_tables)

qlike_plot_A_log <- save_qlike_boxplot(
  Loss_QLIKE_A,
  "QLIKE Distribution - Pipeline A (log scale)",
  file.path(output_figures, "qlike_boxplot_A_log.png"),
  use_log_y = TRUE,
  winsor_p = NULL
)
qlike_plot_A_winsor <- save_qlike_boxplot(
  Loss_QLIKE_A,
  "QLIKE Distribution - Pipeline A (winsorized 97.5%)",
  file.path(output_figures, "qlike_boxplot_A_winsor.png"),
  use_log_y = FALSE,
  winsor_p = 0.975
)
qlike_plot_B <- save_qlike_boxplot(
  Loss_QLIKE_B,
  "QLIKE Distribution - Pipeline B",
  file.path(output_figures, "qlike_boxplot_B.png"),
  use_log_y = FALSE,
  winsor_p = NULL
)
print(qlike_plot_A_log)
print(qlike_plot_A_winsor)
print(qlike_plot_B)

summary_all <- bind_rows(
  mutate(qlike_B_tbl, pipeline = "B"),
  mutate(qlike_A_tbl, pipeline = "A")
) %>%
  dplyr::select(pipeline, model, mean_qlike, sd_qlike, n_oos, included_in_MCS)
write.csv(summary_all, file.path(output_tables, "qlike_summary_all.csv"), row.names = FALSE)

print(pA)
print(pB)
print(lb_A); print(lb_B)
print(mcs_A); print(mcs_B)
if (isTRUE(CONFIG$verbose)) {
  invisible(lapply(fitsA, function(f) { show(f) }))
  invisible(lapply(fitsB, function(f) { show(f) }))
} else {
  capture.output(lapply(names(fitsA), function(m) {
    tryCatch(show(fitsA[[m]]),
             error = function(e) cat(sprintf("ERROR for %s: %s\n", m, e$message)))
  }), file = file.path(output_tables, "fitsA_summary.txt"))
  capture.output(lapply(names(fitsB), function(m) {
    tryCatch(show(fitsB[[m]]),
             error = function(e) cat(sprintf("ERROR for %s: %s\n", m, e$message)))
  }), file = file.path(output_tables, "fitsB_summary.txt"))
}

# ------------------------------ Regime evidence (levels) -----------------------
regime_range <- paste0(CONFIG$regime$start, "/", CONFIG$regime$end)

GLD_regime <- GLDpx[regime_range]
gld_df <- data.frame(Date = index(GLD_regime), Price = as.numeric(GLD_regime))
gld_plot <- ggplot(gld_df, aes(x = Date, y = Price)) +
  geom_line(color = "steelblue", linewidth = 0.7) +
  labs(title = "GLD Adjusted Price (Regime Window)", x = NULL, y = "Adjusted Price") +
  theme_minimal(12)
ggsave(file.path(output_figures, "gld_regime_price.png"), gld_plot, width = 9, height = 4.5, dpi = 150)

GVZ <- na.locf(Cl(GVZ))
VIX <- na.locf(Cl(VIX))
MOVE <- na.locf(Cl(MOVE))
OVX <- na.locf(Cl(OVX))
ds <- na.omit(merge(GVZ, VIX, MOVE, OVX))
colnames(ds) <- c("GVZ", "VIX", "MOVE", "OVX")
ds_regime <- ds[regime_range]

levels_df <- data.frame(Date = index(ds_regime), coredata(ds_regime))
levels_long <- pivot_longer(levels_df, -Date, names_to = "Index", values_to = "Level")
levels_plot <- ggplot(levels_long, aes(x = Date, y = Level, color = Index)) +
  geom_line(linewidth = 0.7) +
  labs(title = "GVZ, VIX, MOVE Levels (Regime Window)", x = NULL, y = "Level", color = NULL) +
  theme_minimal(12)
ggsave(file.path(output_figures, "gvz_vix_move_regime.png"), levels_plot, width = 9, height = 4.5, dpi = 150)

zscore_from <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) return(rep(NA_real_, length(x)))
  (x - mu) / sdv
}

z_df <- levels_df
z_df$GVZ <- zscore_from(z_df$GVZ)
z_df$VIX <- zscore_from(z_df$VIX)
z_df$MOVE <- zscore_from(z_df$MOVE)
z_df$OVX <- zscore_from(z_df$OVX)
z_long <- pivot_longer(z_df, -Date, names_to = "Index", values_to = "Z")
z_plot <- ggplot(z_long, aes(x = Date, y = Z, color = Index)) +
  geom_line(linewidth = 0.7) +
  labs(title = "Z-Score Levels: GVZ, VIX, MOVE, OVX", x = NULL, y = "Z-score", color = NULL) +
  theme_minimal(12)
ggsave(file.path(output_figures, "zscore_gvz_vix_move_ovx.png"), z_plot, width = 9, height = 4.5, dpi = 150)
print(z_plot)

rc30_GVZ_MOVE <- xts(runCor(ds$GVZ, ds$MOVE, n = 30), order.by = index(ds))
rc30_GVZ_VIX <- xts(runCor(ds$GVZ, ds$VIX, n = 30), order.by = index(ds))
rc60_GVZ_MOVE <- xts(runCor(ds$GVZ, ds$MOVE, n = 60), order.by = index(ds))
rc60_GVZ_VIX <- xts(runCor(ds$GVZ, ds$VIX, n = 60), order.by = index(ds))

rc30_df <- data.frame(Date = index(rc30_GVZ_MOVE),
                      GVZ_MOVE = as.numeric(rc30_GVZ_MOVE),
                      GVZ_VIX = as.numeric(rc30_GVZ_VIX)) %>%
  na.omit()
rc30_df <- rc30_df[rc30_df$Date >= CONFIG$regime$start & rc30_df$Date <= CONFIG$regime$end, ]
rc30_long <- pivot_longer(rc30_df, -Date, names_to = "Pair", values_to = "Corr")
rc30_plot <- ggplot(rc30_long, aes(x = Date, y = Corr, color = Pair)) +
  geom_line(linewidth = 0.7) +
  labs(title = "Rolling Corr 30d: GVZ vs MOVE/VIX", x = NULL, y = "Correlation", color = NULL) +
  theme_minimal(12)
ggsave(file.path(output_figures, "rolling_corr_30d.png"), rc30_plot, width = 9, height = 4.5, dpi = 150)

rc60_df <- data.frame(Date = index(rc60_GVZ_MOVE),
                      GVZ_MOVE = as.numeric(rc60_GVZ_MOVE),
                      GVZ_VIX = as.numeric(rc60_GVZ_VIX)) %>%
  na.omit()
rc60_df <- rc60_df[rc60_df$Date >= CONFIG$regime$start & rc60_df$Date <= CONFIG$regime$end, ]
rc60_long <- pivot_longer(rc60_df, -Date, names_to = "Pair", values_to = "Corr")
rc60_plot <- ggplot(rc60_long, aes(x = Date, y = Corr, color = Pair)) +
  geom_line(linewidth = 0.7) +
  labs(title = "Rolling Corr 60d: GVZ vs MOVE/VIX", x = NULL, y = "Correlation", color = NULL) +
  theme_minimal(12)
ggsave(file.path(output_figures, "rolling_corr_60d.png"), rc60_plot, width = 9, height = 4.5, dpi = 150)

pct_30 <- mean(rc30_df$GVZ_MOVE > rc30_df$GVZ_VIX, na.rm = TRUE)
pct_60 <- mean(rc60_df$GVZ_MOVE > rc60_df$GVZ_VIX, na.rm = TRUE)
metric_tbl <- data.frame(
  window = c("30d", "60d"),
  pct_gvz_move_gt_gvz_vix = c(pct_30, pct_60)
)
write.csv(metric_tbl, file.path(output_tables, "rolling_corr_metric.csv"), row.names = FALSE)
cat(sprintf("Percentuale (30d) GVZ-MOVE > GVZ-VIX: %.2f%%\n", 100 * pct_30))
cat(sprintf("Percentuale (60d) GVZ-MOVE > GVZ-VIX: %.2f%%\n", 100 * pct_60))

gld_overlay_df <- data.frame(
  Date = index(GLD_regime),
  GLD = as.numeric(GLD_regime),
  GVZ = as.numeric(ds_regime$GVZ)
)
gld_overlay_df <- na.omit(gld_overlay_df)
gld_min <- min(gld_overlay_df$GLD, na.rm = TRUE)
gld_max <- max(gld_overlay_df$GLD, na.rm = TRUE)
gvz_min <- min(gld_overlay_df$GVZ, na.rm = TRUE)
gvz_max <- max(gld_overlay_df$GVZ, na.rm = TRUE)
scale_factor <- (gld_max - gld_min) / (gvz_max - gvz_min)
gld_overlay_df$GVZ_scaled <- gld_min + (gld_overlay_df$GVZ - gvz_min) * scale_factor
gld_overlay_plot <- ggplot(gld_overlay_df, aes(x = Date)) +
  geom_line(aes(y = GLD, color = "GLD")) +
  geom_line(aes(y = GVZ_scaled, color = "GVZ (scaled)")) +
  scale_y_continuous(
    name = "GLD Adjusted Price",
    sec.axis = sec_axis(~ (. - gld_min) / scale_factor + gvz_min, name = "GVZ")
  ) +
  labs(title = "GLD with GVZ Overlay (Regime Window)", x = NULL, color = NULL) +
  theme_minimal(12)
ggsave(file.path(output_figures, "gld_overlay_gvz.png"), gld_overlay_plot, width = 9, height = 4.5, dpi = 150)
summary_stats <- levels_df %>%
  dplyr::select(-Date) %>%
  dplyr::summarise(dplyr::across(everything(),
                   list(mean = ~ mean(.x, na.rm = TRUE),
                        median = ~ median(.x, na.rm = TRUE),
                        p25 = ~ quantile(.x, 0.25, na.rm = TRUE),
                        p75 = ~ quantile(.x, 0.75, na.rm = TRUE))))
summary_stats <- as.data.frame(t(summary_stats))
summary_stats$Metric <- rownames(summary_stats)
summary_stats <- summary_stats %>%
  tidyr::separate(Metric, into = c("Series", "Statistic"), sep = "_") %>%
  tidyr::pivot_wider(names_from = Statistic, values_from = V1)

gld_start <- as.numeric(first(GLD_regime))
gld_end <- as.numeric(last(GLD_regime))
gld_return <- (gld_end / gld_start) - 1
gld_summary <- data.frame(
  Series = "GLD",
  mean = NA_real_,
  median = NA_real_,
  p25 = NA_real_,
  p75 = NA_real_,
  start_price = gld_start,
  end_price = gld_end,
  pct_change = gld_return,
  stringsAsFactors = FALSE
)

summary_table <- summary_stats %>%
  mutate(start_price = NA_real_, end_price = NA_real_, pct_change = NA_real_) %>%
  bind_rows(gld_summary)

write.csv(summary_table, file.path(output_tables, "regime_summary.csv"), row.names = FALSE)

# ------------------------------ Appendice: Analisi GVZ–MOVE–VIX & EGARCH-X -----------------------

reg_df <- na.omit(merge(ds$GVZ, lag(ds$MOVE, 1), lag(ds$VIX, 1)))
reg_df <- reg_df[regime_range]
colnames(reg_df) <- c("GVZ", "MOVE_l1", "VIX_l1")
reg1 <- lm(GVZ ~ MOVE_l1 + VIX_l1, data = as.data.frame(reg_df))
print(summary(reg1))

ret_GLD_APP <- diff(log(GLDpx))[-1]
ret_GLD_APP <- ret_GLD_APP[regime_range]

X_APP <- merge(lag(ds$MOVE, 1), lag(ds$VIX, 1))
colnames(X_APP) <- c("MOVE_l1", "VIX_l1")
tmp_APP <- na.omit(merge(ret_GLD_APP, X_APP, join = "inner"))
ret_GLD_APP <- tmp_APP[, 1]
X_APP <- as.matrix(tmp_APP[, -1])

spec_EGARCHX_APP <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1),
                        external.regressors = X_APP),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
fit_EGARCHX_APP <- ugarchfit(spec_EGARCHX_APP, data = ret_GLD_APP)
show(fit_EGARCHX_APP)

forecast_len_APP <- round(0.25 * NROW(ret_GLD_APP))
roll_EGARCHX_APP <- ugarchroll(
  spec_EGARCHX_APP, data = ret_GLD_APP,
  n.ahead = 1,
  forecast.length = forecast_len_APP,
  refit.every = 10, refit.window = "moving"
)

s2_EGARCHX_APP <- (function(rr) {
  df <- as.data.frame(rr)
  idx <- if ("Index" %in% names(df)) as.Date(df$Index) else
    if ("Date" %in% names(df)) as.Date(df$Date) else
      as.Date(rownames(df))
  xts((df$Sigma)^2, order.by = idx)
})(roll_EGARCHX_APP)
s2_EGARCHX_APP <- enforce_sigma2_floor(s2_EGARCHX_APP, sigma2_floor_value, "EGARCH_X", "APP")$sigma2

RV_APP <- (ret_GLD_APP)^2
colnames(RV_APP) <- "RV"
combo_APP <- na.omit(merge(s2_EGARCHX_APP, RV_APP, join = "inner"))
s2_vec <- as.numeric(combo_APP[, 1])
rv_vec <- as.numeric(combo_APP[, "RV"])
QLIKE_EGARCHX_APP <- qlike_from(rv_vec, s2_vec)

spec_EGARCH_base <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
fit_EGARCH_base <- ugarchfit(spec_EGARCH_base, data = ret_GLD_APP)
show(fit_EGARCH_base)

roll_EGARCH_base <- ugarchroll(
  spec_EGARCH_base, data = ret_GLD_APP,
  n.ahead = 1,
  forecast.length = forecast_len_APP,
  refit.every = 10, refit.window = "moving"
)

s2_EGARCH_base <- (function(rr) {
  df <- as.data.frame(rr)
  idx <- if ("Index" %in% names(df)) as.Date(df$Index) else
    if ("Date" %in% names(df)) as.Date(df$Date) else
      as.Date(rownames(df))
  xts((df$Sigma)^2, order.by = idx)
})(roll_EGARCH_base)
s2_EGARCH_base <- enforce_sigma2_floor(s2_EGARCH_base, sigma2_floor_value, "EGARCH_base", "APP")$sigma2

combo_base <- na.omit(merge(s2_EGARCH_base, RV_APP, join = "inner"))
s2_base_vec <- as.numeric(combo_base[, 1])
rv_base_vec <- as.numeric(combo_base[, "RV"])
QLIKE_EGARCH_base <- qlike_from(rv_base_vec, s2_base_vec)

egarch_summary <- data.frame(
  model = c("EGARCH_X", "EGARCH_base"),
  mean_qlike = c(mean(QLIKE_EGARCHX_APP, na.rm = TRUE), mean(QLIKE_EGARCH_base, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
egarch_summary$diff_base_minus_x <- egarch_summary$mean_qlike[2] - egarch_summary$mean_qlike[1]
write.csv(egarch_summary, file.path(output_tables, "egarchx_vs_base.csv"), row.names = FALSE)

cat("QLIKE medio EGARCH-X =", mean(QLIKE_EGARCHX_APP, na.rm = TRUE), "\n")
cat("QLIKE medio EGARCH base =", mean(QLIKE_EGARCH_base, na.rm = TRUE), "\n")

reg_coef <- summary(reg1)$coefficients
reg_table <- data.frame(
  term = rownames(reg_coef),
  estimate = reg_coef[, 1],
  p_value = reg_coef[, 4],
  stringsAsFactors = FALSE
)
reg_table$AIC <- AIC(reg1)
reg_table$BIC <- BIC(reg1)

egarch_ic_x <- infocriteria(fit_EGARCHX_APP)
egarch_ic_base <- infocriteria(fit_EGARCH_base)
egarch_table <- data.frame(
  model = c("EGARCH_X", "EGARCH_base"),
  AIC = c(egarch_ic_x[1], egarch_ic_base[1]),
  BIC = c(egarch_ic_x[2], egarch_ic_base[2]),
  mean_qlike = c(mean(QLIKE_EGARCHX_APP, na.rm = TRUE), mean(QLIKE_EGARCH_base, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

conclusions <- c(
  sprintf("30d corr(GVZ,MOVE) > corr(GVZ,VIX) in %.1f%% of days.", 100 * pct_30),
  sprintf("60d corr(GVZ,MOVE) > corr(GVZ,VIX) in %.1f%% of days.", 100 * pct_60),
  sprintf("EGARCH-X mean QLIKE: %.4f.", mean(QLIKE_EGARCHX_APP, na.rm = TRUE)),
  sprintf("EGARCH base mean QLIKE: %.4f.", mean(QLIKE_EGARCH_base, na.rm = TRUE)),
  "Pipeline B (long GLD sample) is the main benchmark; Pipeline A is a short-window stress test.",
  "Regime window shows elevated GVZ/MOVE versus equity volatility.",
  "Gold price dynamics align with macro-risk regime narrative."
)

pdf(file.path(CONFIG$output_dir, "one_pager.pdf"), width = 11, height = 8.5)
print(z_plot)
print(rc30_plot)
print(gld_overlay_plot)
print(qlike_plot_A_log)
grid::grid.newpage()
grid::grid.text(sprintf("OOS size A: %d | OOS size B: %d", NROW(S2_A), NROW(S2_B)),
                x = 0.02, y = 0.98, just = c("left", "top"))
grid::grid.text("Sigma^2 floored percentage by model (A):",
                x = 0.02, y = 0.90, just = c("left", "top"))
grid::grid.text(capture.output(print(floor_summary_source_A[, c("model", "pct_floored")])),
                x = 0.02, y = 0.86, just = c("left", "top"))
grid::grid.text("Sigma^2 floored percentage by model (B):",
                x = 0.02, y = 0.66, just = c("left", "top"))
grid::grid.text(capture.output(print(floor_summary_source_B[, c("model", "pct_floored")])),
                x = 0.02, y = 0.62, just = c("left", "top"))
grid::grid.newpage()
grid::grid.text("Regression summary (GVZ ~ MOVE_l1 + VIX_l1):",
                x = 0.02, y = 0.95, just = c("left", "top"))
grid::grid.text(capture.output(print(reg_table)),
                x = 0.02, y = 0.90, just = c("left", "top"))
grid::grid.text("EGARCH summary (AIC/BIC + mean QLIKE):",
                x = 0.02, y = 0.55, just = c("left", "top"))
grid::grid.text(capture.output(print(egarch_table)),
                x = 0.02, y = 0.50, just = c("left", "top"))
grid::grid.text("Pipeline A sigma^2 floor diagnostics:",
                x = 0.02, y = 0.38, just = c("left", "top"))
grid::grid.text(capture.output(print(loss_A_res$floor_summary)),
                x = 0.02, y = 0.34, just = c("left", "top"))
grid::grid.newpage()
grid::grid.text(paste0(
  "Conclusions:\n- ",
  paste(conclusions, collapse = "\n- "),
  "\n- Pipeline A is a short-window stress test; eGARCH/apARCH can produce near-zero variance forecasts, which inflates QLIKE. Main robustness comes from Pipeline B (long sample)."
), x = 0.02, y = 0.98, just = c("left", "top"))
dev.off()

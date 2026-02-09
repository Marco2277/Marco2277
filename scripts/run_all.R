# Usage: Rscript scripts/run_all.R

source("R/utils_vol.R")
source("R/plots.R")
source("R/pipelines.R")
source("reports/report.R")

config <- list(
  ticker_A = "GLD",
  ticker_B = "GLD",
  start_A = "2025-06-01",
  end_A = "2026-01-19",
  start_B = "2006-01-03",
  end_B = "2025-12-15",
  oos_frac_A = 0.30,
  oos_frac_B = 0.25,
  refit_every_A = 10,
  refit_every_B = 50,
  refit_window_A = "moving",
  refit_window_B = "moving",
  loss = list(min_sigma2 = 1e-10),
  MCS = list(alpha = 0.10, B = 1000)
)

require_pkgs(c("quantmod", "xts", "zoo"))

output_dir <- "output"
fig_dir <- file.path(output_dir, "figures")
tab_dir <- file.path(output_dir, "tables")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}
if (!dir.exists(tab_dir)) {
  dir.create(tab_dir, recursive = TRUE)
}
log_path <- file.path(tab_dir, "log.txt")
if (file.exists(log_path)) {
  file.remove(log_path)
}
show_plots <- interactive() &&
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable()

ensure_rstudio_plots <- function() {
  if (!show_plots) return(FALSE)

  dl <- grDevices::dev.list()
  if (is.null(dl)) {
    grDevices::dev.new()  # in RStudio apre RStudioGD nel Plots
    return(TRUE)
  }

  nms <- names(dl)
  ids <- as.integer(dl)

  if ("RStudioGD" %in% nms) {
    grDevices::dev.set(which = ids[match("RStudioGD", nms)])
    return(TRUE)
  }

  # id does not exist RStudioGD, do ot open other device (windows/quartz/X11)
  FALSE
}

run_plot_safely <- function(plot_fun) {
  if (!show_plots) return(invisible(NULL))

  # Force plots into RStudio Plots (RStudioGD)
  if (!ensure_rstudio_plots()) return(invisible(NULL))

  old_par <- par(no.readonly = TRUE)
  old_par$new <- FALSE
  on.exit(try(par(old_par), silent = TRUE), add = TRUE)

  tryCatch(
    {
      par(
        mfrow = c(1, 1),
        mar = c(4, 4, 2, 1) + 0.1,
        oma = c(0, 0, 0, 0),
        mgp = c(2.2, 0.7, 0),
        las = 1,
        xpd = FALSE,
        new = FALSE
      )

      plot_fun()

      if (exists("dev.flush", envir = asNamespace("grDevices"), inherits = FALSE)) {
        try(grDevices::dev.flush(), silent = TRUE)
      }
    },
    error = function(e) message("Plot failed: ", conditionMessage(e))
  )

  invisible(NULL)
}

save_png <- function(file_path, width, height, plot_fun, res = 150) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename = file_path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_fun()
  invisible(NULL)
}

show_saved_figures <- function(fig_dir = "output/figures") {
  if (!dir.exists(fig_dir)) {
    return(invisible(NULL))
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    warning("Package 'png' is required to display saved figures.")
    return(invisible(NULL))
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    warning("Package 'grid' is required to display saved figures.")
    return(invisible(NULL))
  }
  figs <- list.files(fig_dir, pattern = "\\.png$", full.names = TRUE)
  if (length(figs) == 0) {
    return(invisible(NULL))
  }
  for (path in figs) {
    run_plot_safely(function() {
      img <- png::readPNG(path)
      grid::grid.raster(img)
      grid::grid.text(basename(path), x = 0.01, y = 0.99, just = c("left", "top"))
    })
    readline(prompt = "Press [Enter] to continue...")
  }
}

add_legend <- function(where, labels, cols) {
  legend(
    where,
    legend = labels,
    col = cols,
    lty = 1,
    bty = "n",
    cex = 0.9,
    inset = 0.015
  )
}

winsorize <- function(x, p = 0.99) {
  x <- as.numeric(x)
  qs <- stats::quantile(x, probs = c(1 - p, p), na.rm = TRUE, names = FALSE)
  pmax(pmin(x, qs[2]), qs[1])
}

# External regressors for Pipeline B
reg_env <- new.env()
suppressWarnings(quantmod::getSymbols(c("^VIX", "^MOVE"), src = "yahoo", from = config$start_B, to = config$end_B, env = reg_env))
get_symbol <- function(env, primary, fallback) {
  if (exists(primary, envir = env)) {
    return(env[[primary]])
  }
  if (exists(fallback, envir = env)) {
    return(env[[fallback]])
  }
  stop(sprintf("Symbol not found in env: %s or %s", primary, fallback))
}
VIX <- quantmod::Cl(get_symbol(reg_env, "VIX", "^VIX"))
MOVE <- quantmod::Cl(get_symbol(reg_env, "MOVE", "^MOVE"))
VIX <- zoo::na.locf(VIX, na.rm = FALSE)
MOVE <- zoo::na.locf(MOVE, na.rm = FALSE)
VIX <- stats::na.omit(VIX)
MOVE <- stats::na.omit(MOVE)
regressors <- build_regressors(VIX, MOVE)

# Run pipelines
tmpA <- stats::na.omit(fetch_returns(config$ticker_A, config$start_A, config$end_A))
nA <- NROW(tmpA)
config$oos_target_A <- 30
config$oos_frac_A <- min(max(config$oos_target_A / nA, 0.10), 0.30)
pipeline_A <- run_pipeline_A(config)
pipeline_B <- run_pipeline_B(config, regressors = regressors, log_path = log_path)

# Save sigma2 panels
floor_value <- max(config$loss$min_sigma2, .Machine$double.eps)
save_sigma_panel(pipeline_A$S2, file.path(fig_dir, "conditional_variance_A.png"),
                 "Conditional Variance (Pipeline A)", floor_value)
save_sigma_panel(pipeline_B$S2, file.path(fig_dir, "conditional_variance_B.png"),
                 "Conditional Variance (Pipeline B)", floor_value)
# cat("dev.cur =", dev.cur(), "name =", device_name_for_id(dev.cur()), "\n")
# print(grDevices::dev.list())
if (show_plots) {
  run_plot_safely(function() {
    make_sigma_panel(pipeline_A$S2, "Conditional Variance (Pipeline A)", floor_value, compact = TRUE)
  })
  run_plot_safely(function() {
    make_sigma_panel(pipeline_B$S2, "Conditional Variance (Pipeline B)", floor_value, compact = TRUE)
  })
}

# NIC plots
save_nic_panel(pipeline_A$nic_list, file.path(fig_dir, "nic_A.png"), "NIC (Pipeline A)")
save_nic_panel(pipeline_B$nic_list, file.path(fig_dir, "nic_B.png"), "NIC (Pipeline B)")
if (show_plots) {
  run_plot_safely(function() {
    plot_nic(pipeline_A$nic_list, "NIC (Pipeline A)")
  })
  run_plot_safely(function() {
    plot_nic(pipeline_B$nic_list, "NIC (Pipeline B)")
  })
}

# Ljung-Box tables
write.csv(pipeline_A$lb_table, file.path(tab_dir, "lb_A.csv"), row.names = FALSE)
write.csv(pipeline_B$lb_table, file.path(tab_dir, "lb_B.csv"), row.names = FALSE)

# QLIKE tables
qlike_means_A <- data.frame(model = names(pipeline_A$Loss),
                            mean_qlike = vapply(pipeline_A$Loss, function(x) mean(x, na.rm = TRUE), numeric(1)))
qlike_means_B <- data.frame(model = names(pipeline_B$Loss),
                            mean_qlike = vapply(pipeline_B$Loss, function(x) mean(x, na.rm = TRUE), numeric(1)))
write.csv(qlike_means_A, file.path(tab_dir, "qlike_means_A.csv"), row.names = FALSE)
write.csv(qlike_means_B, file.path(tab_dir, "qlike_means_B.csv"), row.names = FALSE)

# Sigma2 floor diagnostics (summary)
fdA <- if (length(pipeline_A$floor_diag) > 0) do.call(rbind, pipeline_A$floor_diag) else data.frame()
fdB <- if (length(pipeline_B$floor_diag) > 0) do.call(rbind, pipeline_B$floor_diag) else data.frame()
write.csv(fdA, file.path(tab_dir, "sigma2_floor_diagnostics_A.csv"), row.names = FALSE)
write.csv(fdB, file.path(tab_dir, "sigma2_floor_diagnostics_B.csv"), row.names = FALSE)

# Sigma2 floor dates (exact days)
fdtA <- if (length(pipeline_A$floor_dates) > 0) do.call(rbind, pipeline_A$floor_dates) else data.frame()
fdtB <- if (length(pipeline_B$floor_dates) > 0) do.call(rbind, pipeline_B$floor_dates) else data.frame()
write.csv(fdtA, file.path(tab_dir, "sigma2_floor_dates_A.csv"), row.names = FALSE)
write.csv(fdtB, file.path(tab_dir, "sigma2_floor_dates_B.csv"), row.names = FALSE)

# QLIKE outlier diagnostics (top 3 worst days per model)
write.csv(pipeline_A$outlier_diag, file.path(tab_dir, "pipelineA_outlier_diagnostics.csv"), row.names = FALSE)
write.csv(pipeline_B$outlier_diag, file.path(tab_dir, "pipelineB_outlier_diagnostics.csv"), row.names = FALSE)

# QLIKE summaries
qlike_summary_all <- rbind(
  cbind(pipeline = "A", qlike_means_A),
  cbind(pipeline = "B", qlike_means_B)
)
write.csv(qlike_summary_all, file.path(tab_dir, "qlike_summary_all.csv"), row.names = FALSE)

# QLIKE boxplots
lossA <- pipeline_A$Loss_matrix
lossA_log <- log(pmax(lossA, 1e-12))
lossA_winsor <- apply(lossA, 2, winsorize, p = 0.99)
save_qlike_boxplot(
  lossA_log,
  file.path(fig_dir, "qlike_boxplot_A_log.png"),
  "QLIKE (A, log)",
  mar = c(14, 4.4, 4.0, 2.0) + 0.1,
  cex_axis = 0.60
)

save_qlike_boxplot(
  lossA_winsor,
  file.path(fig_dir, "qlike_boxplot_A_winsor.png"),
  "QLIKE (A, winsor)",
  mar = c(14, 4.4, 4.0, 2.0) + 0.1,
  cex_axis = 0.60
)

save_qlike_boxplot(
  pipeline_B$Loss_matrix,
  file.path(fig_dir, "qlike_boxplot_B.png"),
  "QLIKE (B)",
  mar = c(16, 4.4, 4.0, 2.0) + 0.1,
  cex_axis = 0.55
)

plot_qlike_boxplot_rstudio <- function(x, main_title) {
  par(
    mar = c(6, 4.2, 2.8, 1.2) + 0.1,  
    las = 2,                         
    tcl = -0.2,
    xpd = NA,
    cex.axis = 0.50,                 
    cex.main = 0.95
  )
  boxplot(x, main = main_title, las = 2, col = "steelblue")
}

if (show_plots) {
  run_plot_safely(function() {
    plot_qlike_boxplot_rstudio(lossA_log, "QLIKE (A, log)")
  })
  run_plot_safely(function() {
    plot_qlike_boxplot_rstudio(lossA_winsor, "QLIKE (A, winsor)")
  })
  run_plot_safely(function() {
    plot_qlike_boxplot_rstudio(pipeline_B$Loss_matrix, "QLIKE (B)")
  })
}

# MCS summaries
mcs_summary_A <- data.frame(model = pipeline_A$MCS_models)
mcs_summary_B <- data.frame(model = pipeline_B$MCS_models)
write.csv(mcs_summary_A, file.path(tab_dir, "mcs_summary_A.csv"), row.names = FALSE)
write.csv(mcs_summary_B, file.path(tab_dir, "mcs_summary_B.csv"), row.names = FALSE)
if (!is.null(pipeline_A$MCS)) {
  mcs_table_A <- try(methods::slot(pipeline_A$MCS, "show"), silent = TRUE)
  if (inherits(mcs_table_A, "try-error") || is.null(mcs_table_A)) {
    mcs_table_A <- as.data.frame(pipeline_A$MCS)
  }
  write.csv(mcs_table_A, file.path(tab_dir, "mcs_table_A.csv"), row.names = TRUE)
}
if (!is.null(pipeline_B$MCS)) {
  mcs_table_B <- try(methods::slot(pipeline_B$MCS, "show"), silent = TRUE)
  if (inherits(mcs_table_B, "try-error") || is.null(mcs_table_B)) {
    mcs_table_B <- as.data.frame(pipeline_B$MCS)
  }
  write.csv(mcs_table_B, file.path(tab_dir, "mcs_table_B.csv"), row.names = TRUE)
}

# Regime plots and regime evidence
suppressWarnings(quantmod::getSymbols(c("GLD", "^GVZ", "^OVX"), src = "yahoo", from = config$start_B, to = config$end_B, env = reg_env))
GLD <- quantmod::Cl(get_symbol(reg_env, "GLD", "GLD"))
GVZ <- quantmod::Cl(get_symbol(reg_env, "GVZ", "^GVZ"))
OVX <- quantmod::Cl(get_symbol(reg_env, "OVX", "^OVX"))
GLD <- zoo::na.locf(GLD, na.rm = FALSE)
GVZ <- zoo::na.locf(GVZ, na.rm = FALSE)
VIX <- zoo::na.locf(VIX, na.rm = FALSE)
MOVE <- zoo::na.locf(MOVE, na.rm = FALSE)
OVX <- zoo::na.locf(OVX, na.rm = FALSE)
GLD <- stats::na.omit(GLD)
GVZ <- stats::na.omit(GVZ)
VIX <- stats::na.omit(VIX)
MOVE <- stats::na.omit(MOVE)
OVX <- stats::na.omit(OVX)

z_gvz <- (GVZ - mean(GVZ, na.rm = TRUE)) / stats::sd(GVZ, na.rm = TRUE)
z_vix <- (VIX - mean(VIX, na.rm = TRUE)) / stats::sd(VIX, na.rm = TRUE)
z_move <- (MOVE - mean(MOVE, na.rm = TRUE)) / stats::sd(MOVE, na.rm = TRUE)
z_ovx <- (OVX - mean(OVX, na.rm = TRUE)) / stats::sd(OVX, na.rm = TRUE)

plot_gld_regime_price <- function() {
  x <- as.Date(zoo::index(GLD))
  y <- as.numeric(GLD)
  ok <- is.finite(y)
  if (!any(ok)) {
    plot.new()
    title(main = "GLD Price")
    return()
  }
  plot(x[ok], y[ok], main = "GLD Price", col = "steelblue", type = "l", xlab = "Date", ylab = "Price")
}
merge_inner <- function(...) {
  Reduce(function(a, b) merge(a, b, join = "inner"), list(...))
}

add_legend_top <- function(labels, cols, ncol = length(labels), cex = 0.9) {
  usr <- par("usr")
  x <- mean(usr[1:2])
  y <- usr[4] + 0.08 * diff(usr[3:4])     # 8% sopra il bordo alto
  legend(x, y, legend = labels, col = cols, lty = 1, bty = "n",
         horiz = TRUE, xjust = 0.5, yjust = 0, xpd = NA,
         ncol = ncol, cex = cex)
}

plot_gvz_vix_move_regime <- function() {
  merged <- merge_inner(GVZ, VIX, MOVE)
  x <- as.Date(zoo::index(merged))
  gvz <- as.numeric(merged[, 1])
  vix <- as.numeric(merged[, 2])
  move <- as.numeric(merged[, 3])

  ok <- is.finite(gvz) & is.finite(vix) & is.finite(move)
  if (!any(ok)) { plot.new(); title(main = "GVZ/VIX/MOVE"); return() }

  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mar = c(5, 4, 5, 5) + 0.1, xpd = NA)

  ylim_left <- range(c(gvz[ok], vix[ok]), na.rm = TRUE)
  plot(x[ok], gvz[ok], type = "l", col = "steelblue",
       main = "GVZ/VIX/MOVE", xlab = "Date", ylab = "GVZ / VIX",
       ylim = ylim_left)
  lines(x[ok], vix[ok], col = "firebrick")

  par(new = TRUE)
  ylim_right <- range(move[ok], na.rm = TRUE)
  plot(x[ok], move[ok], type = "l", col = "darkgreen",
       axes = FALSE, xlab = "", ylab = "", ylim = ylim_right)
  axis(4); mtext("MOVE", side = 4, line = 3)

  add_legend_top(c("GVZ", "VIX", "MOVE"),
                 c("steelblue", "firebrick", "darkgreen"),
                 ncol = 3)
}

plot_zscore_gvz_vix_move_ovx <- function() {
  merged <- merge_inner(z_gvz, z_vix, z_move, z_ovx)
  x <- as.Date(zoo::index(merged))
  gvz <- as.numeric(merged[, 1])
  vix <- as.numeric(merged[, 2])
  move <- as.numeric(merged[, 3])
  ovx <- as.numeric(merged[, 4])

  ok <- is.finite(gvz) & is.finite(vix) & is.finite(move) & is.finite(ovx)
  if (!any(ok)) { plot.new(); title(main = "Z-scores"); return() }

  ylim <- range(c(gvz[ok], vix[ok], move[ok], ovx[ok]), na.rm = TRUE)
  pad <- 0.05 * diff(ylim)
  if (is.finite(pad) && pad > 0) ylim <- c(ylim[1] - pad, ylim[2] + pad)

  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mar = c(5, 4, 5, 2) + 0.1, xpd = FALSE)

  plot(x[ok], gvz[ok], col = "steelblue", main = "Z-scores", ylab = "z", xlab = "Date",
       type = "l", ylim = ylim)
  lines(x[ok], vix[ok], col = "firebrick")
  lines(x[ok], move[ok], col = "darkgreen")
  lines(x[ok], ovx[ok], col = "purple")

  add_legend_top(c("GVZ", "VIX", "MOVE", "OVX"),
                 c("steelblue", "firebrick", "darkgreen", "purple"),
                 ncol = 4)
}

save_and_show_png(file.path(fig_dir, "gld_regime_price.png"), 1200, 700, plot_gld_regime_price)
save_and_show_png(file.path(fig_dir, "gvz_vix_move_regime.png"), 1200, 700, plot_gvz_vix_move_regime)
save_and_show_png(file.path(fig_dir, "zscore_gvz_vix_move_ovx.png"), 1200, 700, plot_zscore_gvz_vix_move_ovx)

rolling_corr <- function(x, y, window) {
  merged <- merge(x, y, join = "inner")
  if (NROW(merged) < 1) {
    return(xts::xts(numeric(0), order.by = zoo::index(merged)))
  }
  vals <- zoo::rollapplyr(
    merged,
    width = window,
    FUN = function(mat) stats::cor(mat[, 1], mat[, 2], use = "complete.obs"),
    by.column = FALSE,
    fill = NA_real_
  )
  xts::xts(as.numeric(vals), order.by = zoo::index(merged))
}

corr_30_move <- rolling_corr(GVZ, MOVE, 30)
corr_30_vix <- rolling_corr(GVZ, VIX, 30)
corr_60_move <- rolling_corr(GVZ, MOVE, 60)
corr_60_vix <- rolling_corr(GVZ, VIX, 60)

plot_rolling_corr_30d <- function() {
  merged <- merge_inner(corr_30_move, corr_30_vix)
  x <- as.Date(zoo::index(merged))
  move <- as.numeric(merged[, 1])
  vix  <- as.numeric(merged[, 2])

  ok <- is.finite(move) & is.finite(vix)
  if (!any(ok)) { plot.new(); title(main = "Rolling Corr (30d)"); return() }

  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mar = c(5, 4, 5, 2) + 0.1, xpd = NA)

  plot(x[ok], move[ok], main = "Rolling Corr (30d)", col = "steelblue", type = "l",
       xlab = "Date", ylab = "Corr", ylim = c(-1, 1))
  lines(x[ok], vix[ok], col = "firebrick")
  abline(h = 0, col = "gray70")

  add_legend_top(c("GVZ-MOVE", "GVZ-VIX"),
                 c("steelblue", "firebrick"),
                 ncol = 2)
}

plot_rolling_corr_60d <- function() {
  merged <- merge_inner(corr_60_move, corr_60_vix)
  x <- as.Date(zoo::index(merged))
  move <- as.numeric(merged[, 1])
  vix  <- as.numeric(merged[, 2])

  ok <- is.finite(move) & is.finite(vix)
  if (!any(ok)) { plot.new(); title(main = "Rolling Corr (60d)"); return() }

  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mar = c(5, 4, 5, 2) + 0.1, xpd = NA)

  plot(x[ok], move[ok], main = "Rolling Corr (60d)", col = "steelblue", type = "l",
       xlab = "Date", ylab = "Corr", ylim = c(-1, 1))
  lines(x[ok], vix[ok], col = "firebrick")
  abline(h = 0, col = "gray70")

  add_legend_top(c("GVZ-MOVE", "GVZ-VIX"),
                 c("steelblue", "firebrick"),
                 ncol = 2)
}

save_and_show_png(file.path(fig_dir, "rolling_corr_30d.png"), 1200, 700, plot_rolling_corr_30d)
save_and_show_png(file.path(fig_dir, "rolling_corr_60d.png"), 1200, 700, plot_rolling_corr_60d)

calc_pct_move_gt_vix <- function(corr_move, corr_vix) {
  merged <- merge(corr_move, corr_vix, join = "inner")
  mean(merged[, 1] > merged[, 2], na.rm = TRUE)
}

rolling_corr_metric <- data.frame(
  window = c(30, 60),
  pct_move_gt_vix = c(calc_pct_move_gt_vix(corr_30_move, corr_30_vix),
                      calc_pct_move_gt_vix(corr_60_move, corr_60_vix)),
  mean_corr_gvz_move = c(mean(corr_30_move, na.rm = TRUE), mean(corr_60_move, na.rm = TRUE)),
  mean_corr_gvz_vix = c(mean(corr_30_vix, na.rm = TRUE), mean(corr_60_vix, na.rm = TRUE))
)
write.csv(rolling_corr_metric, file.path(tab_dir, "rolling_corr_metric.csv"), row.names = FALSE)

plot_gld_overlay_gvz <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4, 5, 5) + 0.1, xpd = NA)

  merged <- merge_inner(GLD, GVZ)
  x <- as.Date(zoo::index(merged))
  gld <- as.numeric(merged[, 1])
  gvz <- as.numeric(merged[, 2])

  ok <- is.finite(gld) & is.finite(gvz)
  if (!any(ok)) { plot.new(); title(main = "GLD vs GVZ"); return() }

  plot(x[ok], gld[ok], col = "steelblue", main = "GLD vs GVZ", type = "l", xlab = "Date", ylab = "GLD")

  par(new = TRUE)
  plot(x[ok], gvz[ok], col = "firebrick", axes = FALSE, xlab = "", ylab = "", type = "l")
  axis(side = 4)
  mtext("GVZ", side = 4, line = 3)

  add_legend_top(c("GLD", "GVZ"), c("steelblue", "firebrick"), ncol = 2)
}
save_and_show_png(file.path(fig_dir, "gld_overlay_gvz.png"), 1200, 700, plot_gld_overlay_gvz)

regression_data <- merge(GVZ, build_regressors(VIX, MOVE), join = "inner")
colnames(regression_data)[1] <- "GVZ"
regression_df <- data.frame(
  GVZ = as.numeric(regression_data$GVZ),
  MOVE_l1 = as.numeric(regression_data$MOVE_l1),
  VIX_l1 = as.numeric(regression_data$VIX_l1)
)
reg_fit <- stats::lm(GVZ ~ MOVE_l1 + VIX_l1, data = regression_df)
reg_summary <- summary(reg_fit)
regime_summary <- data.frame(
  term = rownames(reg_summary$coefficients),
  estimate = reg_summary$coefficients[, "Estimate"],
  std_error = reg_summary$coefficients[, "Std. Error"],
  t_value = reg_summary$coefficients[, "t value"],
  p_value = reg_summary$coefficients[, "Pr(>|t|)"],
  r_squared = reg_summary$r.squared
)
write.csv(regime_summary, file.path(tab_dir, "regime_summary.csv"), row.names = FALSE)

# one-pager
write_one_pager(output_dir, list(pipeline_A = pipeline_A, pipeline_B = pipeline_B))

check_outputs <- function() {
  required_files <- c(
  "output/tables/lb_A.csv",
  "output/tables/lb_B.csv",
  "output/tables/mcs_summary_A.csv",
  "output/tables/mcs_summary_B.csv",
  "output/tables/qlike_means_A.csv",
  "output/tables/qlike_means_B.csv",
  "output/tables/qlike_summary_all.csv",
  "output/tables/sigma2_floor_diagnostics_A.csv",
  "output/tables/sigma2_floor_diagnostics_B.csv",
  "output/tables/sigma2_floor_dates_A.csv",
  "output/tables/sigma2_floor_dates_B.csv",
  "output/tables/pipelineA_outlier_diagnostics.csv",
  "output/tables/pipelineB_outlier_diagnostics.csv",
  "output/tables/rolling_corr_metric.csv",
  "output/tables/regime_summary.csv"
)

  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    warning(sprintf("Missing output files: %s", paste(missing_files, collapse = ", ")))
  } else {
    message("All required output files exist.")
  }

  if (is.null(pipeline_B$MCS)) {
    warning("MCS object is NULL")
  } else if (length(pipeline_B$MCS_models) < 2) {
    warning("MCS model extraction returned <2 models; check get_mcs_models")
  } else {
    message("MCS has at least 2 valid models.")
  }

  if (!"EGARCH_X" %in% names(pipeline_B$S2)) {
    message("EGARCH_X excluded from Pipeline B (see log for reason).")
  } else {
    message("EGARCH_X included in Pipeline B.")
  }
}

check_outputs()

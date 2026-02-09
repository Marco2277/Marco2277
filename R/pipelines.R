# Volatility contest pipelines.

align_series <- function(x, y) {
  merge(x, y, join = "inner")
}

build_regressors <- function(vix, move) {
  vix_l1 <- xts::lag.xts(vix, k = 1)
  move_l1 <- xts::lag.xts(move, k = 1)
  merged <- merge(vix_l1, move_l1, join = "inner")
  colnames(merged) <- c("VIX_l1", "MOVE_l1")
  merged
}

align_returns_regressors <- function(ret_xts, regressors_xts) {
  merged <- merge(ret_xts, regressors_xts, join = "inner")
  ret_aligned <- merged[, "ret", drop = FALSE]
  xreg <- merged[, c("VIX_l1", "MOVE_l1"), drop = FALSE]
  list(returns = ret_aligned, regressors = xreg)
}

build_model_specs <- function(regressors = NULL) {
  require_pkgs(c("rugarch"))
  specs <- list(
    sGARCH = rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                                distribution.model = "norm"),
    eGARCH = rugarch::ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                                distribution.model = "norm"),
    apARCH = rugarch::ugarchspec(variance.model = list(model = "apARCH", garchOrder = c(1, 1)),
                                 mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                                 distribution.model = "norm"),
    iGARCH = rugarch::ugarchspec(variance.model = list(model = "iGARCH", garchOrder = c(1, 1)),
                                 mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                                 distribution.model = "norm"),
    fiGARCH = rugarch::ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1)),
                                  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                                  distribution.model = "norm")
  )
  if (!is.null(regressors)) {
    specs$EGARCH_X <- rugarch::ugarchspec(
      variance.model = list(model = "eGARCH", garchOrder = c(1, 1), external.regressors = as.matrix(regressors)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      distribution.model = "norm"
    )
  }
  specs
}

fit_roll_model <- function(spec, returns, oos, refit_every, refit_window) {
  require_pkgs(c("rugarch"))
  rugarch::ugarchroll(
    spec = spec,
    data = returns,
    n.ahead = 1,
    forecast.length = oos,
    refit.every = refit_every,
    refit.window = refit_window,
    solver = "hybrid",
    calculate.VaR = FALSE,
    keep.coef = TRUE
  )
}

extract_sigma2_roll <- function(roll) {
  if (is.null(roll)) return(NULL)

  sigma <- NULL
  idx <- NULL

  if (methods::is(roll, "uGARCHroll")) {
    fc <- try(roll@forecast, silent = TRUE)
    if (!inherits(fc, "try-error") && is.list(fc) && !is.null(fc$density)) {
      dens <- fc$density
      if (is.matrix(dens) || is.data.frame(dens)) {
        cn <- colnames(dens)
        j <- grep("^Sigma$", cn, ignore.case = TRUE)
        if (length(j) == 0) j <- grep("Sigma", cn, ignore.case = TRUE)
        if (length(j) > 0) sigma <- suppressWarnings(as.numeric(dens[, j[1]]))

        rn <- rownames(dens)
        if (!is.null(rn)) {
          idx <- suppressWarnings(as.Date(rn))
          if (all(is.na(idx))) idx <- suppressWarnings(as.Date(substr(rn, 1, 10)))
        }
      }
    }
  }

  if (is.null(sigma) || all(!is.finite(sigma)) || is.null(idx) || all(is.na(idx))) {
    df <- try(as.data.frame(roll), silent = TRUE)
    if (inherits(df, "try-error") || !is.data.frame(df) || nrow(df) == 0) df <- NULL

    if (is.null(df)) {
      rep <- try(rugarch::report(roll), silent = TRUE)
      if (!inherits(rep, "try-error")) {
        df2 <- try(as.data.frame(rep), silent = TRUE)
        if (!inherits(df2, "try-error") && is.data.frame(df2) && nrow(df2) > 0) df <- df2
      }
    }

    if (!is.null(df) && nrow(df) > 0) {
      cn <- names(df)
      j <- grep("^Sigma$", cn, ignore.case = TRUE)
      if (length(j) == 0) j <- grep("Sigma", cn, ignore.case = TRUE)
      if (length(j) > 0) sigma <- suppressWarnings(as.numeric(df[[j[1]]]))

      if ("Index" %in% cn) {
        idx <- suppressWarnings(as.Date(df$Index))
      } else if ("Date" %in% cn) {
        idx <- suppressWarnings(as.Date(df$Date))
      } else {
        rn <- rownames(df)
        idx <- suppressWarnings(as.Date(rn))
        if (all(is.na(idx))) idx <- suppressWarnings(as.Date(substr(rn, 1, 10)))
      }
    }
  }

  ok <- is.finite(sigma) & !is.na(idx)
  if (!any(ok)) return(NULL)

  xts::xts((sigma[ok])^2, order.by = idx[ok])
}

run_vol_contest_pipeline <- function(config, ticker, from, to, oos_frac, refit_every, refit_window,
                                     include_regressors = FALSE, regressors = NULL,
                                     include_riskmetrics = FALSE, pipeline_label = "",
                                     log_path = NULL) {
  require_pkgs(c("xts", "rugarch", "quantmod", "zoo"))
  ret_xts <- fetch_returns(ticker, from, to)
  ret_xts <- stats::na.omit(ret_xts)
  if (include_regressors) {
    aligned <- align_returns_regressors(ret_xts, regressors)
    ret_xts <- aligned$returns
    regressors <- aligned$regressors
    ret_xts <- stats::na.omit(ret_xts)
    regressors <- stats::na.omit(regressors)
  }
  ret_series <- ret_xts[, "ret", drop = FALSE]
  ret_vec <- as.numeric(ret_series$ret)
  n <- length(ret_vec)
  if (n < 2) {
    msg <- sprintf("%s: not enough returns to run pipeline (n=%d)", pipeline_label, n)
    message(msg)
    if (!is.null(log_path)) {
      cat(msg, file = log_path, sep = "\n", append = TRUE)
    }
    return(list(
      returns = ret_xts,
      rv = ret_xts$ret^2,
      S2 = list(),
      Loss = list(),
      Loss_matrix = matrix(numeric(0), nrow = 0, ncol = 0),
      MCS = NULL,
      MCS_models = character(0),
      fits = list(),
      nic_list = list(),
      lb_table = data.frame(),
      floor_diag = list(),
      floor_dates = list(),
      outlier_diag = data.frame(
        model = character(0),
        top1_date = character(0), top1_qlike = numeric(0),
        top2_date = character(0), top2_qlike = numeric(0),
        top3_date = character(0), top3_qlike = numeric(0),
        stringsAsFactors = FALSE
      ),
      n_start = 0
    ))
  }
  oos <- max(1, round(oos_frac * n))
  oos <- min(oos, n - 1)
  n_start <- n - oos
  refit_every <- min(refit_every, max(1, oos))

  rv_xts <- ret_xts$ret^2
  rv_shift <- xts::lag.xts(rv_xts, k = -1)

  specs <- build_model_specs(regressors = if (include_regressors) regressors else NULL)
  floor_value <- max(config$loss$min_sigma2, .Machine$double.eps)

  sigma2_list <- list()
  fits <- list()
  lb_table <- list()
  nic_list <- list()
  floor_diag <- list()
  floor_dates <- list()
  loss_list <- list()
  status <- list()

  model_names <- names(specs)
  if (pipeline_label == "Pipeline A") {
    model_names <- intersect(model_names, c("sGARCH", "eGARCH", "apARCH"))
  }
  for (nm in model_names) {
    status[[nm]] <- list(
      model = nm,
      fitted = FALSE,
      rolled = FALSE,
      sigma2_nonempty = FALSE,
      loss_nonempty = FALSE,
      included_in_loss_matrix = FALSE,
      included_in_mcs_input = FALSE,
      reason = "other"
    )
    spec <- specs[[nm]]
    fit <- NULL
    roll <- NULL
    if (include_regressors && nm == "EGARCH_X") {
      fit <- try(rugarch::ugarchfit(spec = spec, data = ret_series, solver = "hybrid"), silent = TRUE)
      if (inherits(fit, "try-error")) {
        msg <- sprintf("%s: EGARCH_X failed to fit; excluding", pipeline_label)
        message(msg)
        if (!is.null(log_path)) {
          cat(msg, file = log_path, sep = "\n", append = TRUE)
        }
        status[[nm]]$reason <- "fit_failed"
        next
      }
      status[[nm]]$fitted <- TRUE
      roll <- try(fit_roll_model(spec, ret_series, oos = oos, refit_every = refit_every,
                                 refit_window = refit_window), silent = TRUE)
      if (inherits(roll, "try-error")) {
        msg <- sprintf("%s: EGARCH_X failed to roll; excluding", pipeline_label)
        message(msg)
        if (!is.null(log_path)) {
          cat(msg, file = log_path, sep = "\n", append = TRUE)
        }
        status[[nm]]$reason <- "roll_failed"
        next
      }
      status[[nm]]$rolled <- TRUE
    } else {
      fit <- try(rugarch::ugarchfit(spec = spec, data = ret_series, solver = "hybrid"), silent = TRUE)
      if (inherits(fit, "try-error")) {
        msg <- sprintf("%s: %s failed to fit; excluding", pipeline_label, nm)
        message(msg)
        if (!is.null(log_path)) {
          cat(msg, file = log_path, sep = "\n", append = TRUE)
        }
        status[[nm]]$reason <- "fit_failed"
        next
      }
      status[[nm]]$fitted <- TRUE
      roll <- try(fit_roll_model(spec, ret_series, oos = oos, refit_every = refit_every,
                                 refit_window = refit_window), silent = TRUE)
      if (inherits(roll, "try-error")) {
        msg <- sprintf("%s: %s failed to roll; excluding", pipeline_label, nm)
        message(msg)
        if (!is.null(log_path)) {
          cat(msg, file = log_path, sep = "\n", append = TRUE)
        }
        status[[nm]]$reason <- "roll_failed"
        next
      }
      status[[nm]]$rolled <- TRUE
    }
    s2_xts <- extract_sigma2_roll(roll)
    if (is.null(s2_xts) || NROW(s2_xts) == 0 || all(is.na(s2_xts))) {
      msg <- sprintf("%s: %s produced empty sigma2 series; excluding", pipeline_label, nm)
      message(msg)
      if (!is.null(log_path)) {
        cat(msg, file = log_path, sep = "\n", append = TRUE)
      }
      status[[nm]]$reason <- "empty_sigma2"
      next
    }
    status[[nm]]$sigma2_nonempty <- TRUE
    s2_raw <- as.numeric(s2_xts)
    s2_used <- enforce_sigma2_floor(s2_raw, floor_value = floor_value)
    s2_xts <- xts::xts(s2_used, order.by = zoo::index(s2_xts))
    sigma2_list[[nm]] <- s2_xts
    fits[[nm]] <- fit
    floor_diag[[nm]] <- sigma2_floor_diagnostics(s2_raw, s2_used, floor_value)

    rv_aligned <- merge(rv_shift, s2_xts, join = "right")[, 1]

    floor_dates[[nm]] <- sigma2_floor_dates(
      idx = zoo::index(s2_xts),
      sigma2_raw = s2_raw,
      sigma2_used = s2_used,
      floor_value = floor_value,
      pipeline_label = pipeline_label,
      model_name = nm,
      rv_aligned = rv_aligned
    )

    if (!is.null(fit)) {
      z <- get_std_residuals(fit)

      if (!is.null(z) && length(z) > 0) {
        lb_z  <- try(stats::Box.test(z,   lag = 20, type = "Ljung-Box"), silent = TRUE)
        lb_z2 <- try(stats::Box.test(z^2, lag = 20, type = "Ljung-Box"), silent = TRUE)
  
        if (!inherits(lb_z, "try-error") && !inherits(lb_z2, "try-error")) {
          lb_table[[nm]] <- data.frame(
            model = nm,
            lb_stat_z  = lb_z$statistic,
            lb_p_z     = lb_z$p.value,
            lb_stat_z2 = lb_z2$statistic,
            lb_p_z2    = lb_z2$p.value
          )
      }
    } else {
      msg <- sprintf("%s: %s standardized residuals not available; skipping Ljung-Box (NIC may also fail)", pipeline_label, nm)
      message(msg)
      if (!is.null(log_path)) cat(msg, file = log_path, sep = "\n", append = TRUE)
    }

      nic_tmp <- try(compute_nic(fit, log_path = log_path, model_name = nm), silent = TRUE)
      if (!inherits(nic_tmp, "try-error")) nic_list[[nm]] <- nic_tmp
    }

    merged_loss <- merge(rv_shift, s2_xts, join = "inner")
    merged_loss <- stats::na.omit(merged_loss)
    if (NROW(merged_loss) < 1) {
      msg <- sprintf("%s: %s has no aligned RV/sigma2 rows; excluding", pipeline_label, nm)
      message(msg)
      if (!is.null(log_path)) {
        cat(msg, file = log_path, sep = "\n", append = TRUE)
      }
      status[[nm]]$reason <- "no_aligned_loss"
      next
    }
    loss_series <- qlike_from(merged_loss[, 1], merged_loss[, 2])
    loss_xts <- xts::xts(matrix(loss_series, ncol = 1), order.by = zoo::index(merged_loss))
    colnames(loss_xts) <- nm
    loss_list[[nm]] <- loss_xts
    status[[nm]]$loss_nonempty <- NROW(loss_xts) > 0
    status[[nm]]$reason <- "ok"
  }

  if (include_riskmetrics) {
    status[["RiskMetrics"]] <- list(
      model = "RiskMetrics",
      fitted = TRUE,
      rolled = TRUE,
      sigma2_nonempty = FALSE,
      loss_nonempty = FALSE,
      included_in_loss_matrix = FALSE,
      included_in_mcs_input = FALSE,
      reason = "other"
    )
    s2_rm <- ewma_sigma2(ret_vec, lambda = 0.94)
    s2_rm <- s2_rm[(n_start + 1):n]
    idx_rm <- zoo::index(ret_xts)[(n_start + 1):n]
    s2_raw <- as.numeric(s2_rm)
    s2_used <- enforce_sigma2_floor(s2_raw, floor_value = floor_value)

    s2_rm_xts <- xts::xts(s2_used, order.by = idx_rm)
    sigma2_list[["RiskMetrics"]] <- s2_rm_xts

    floor_diag[["RiskMetrics"]] <- sigma2_floor_diagnostics(s2_raw, s2_used, floor_value)

    rv_aligned <- merge(rv_shift, s2_rm_xts, join = "right")[, 1]

    floor_dates[["RiskMetrics"]] <- sigma2_floor_dates(
      idx = zoo::index(s2_rm_xts),
      sigma2_raw = s2_raw,
      sigma2_used = s2_used,
      floor_value = floor_value,
      pipeline_label = pipeline_label,
      model_name = "RiskMetrics",
      rv_aligned = rv_aligned
    )
    status[["RiskMetrics"]]$sigma2_nonempty <- NROW(s2_rm_xts) > 0
    merged_loss <- merge(rv_shift, s2_rm_xts, join = "inner")
    merged_loss <- stats::na.omit(merged_loss)
    if (NROW(merged_loss) < 1) {
      status[["RiskMetrics"]]$reason <- "no_aligned_loss"
    } else {
      loss_series <- qlike_from(merged_loss[, 1], merged_loss[, 2])
      loss_xts <- xts::xts(matrix(loss_series, ncol = 1), order.by = zoo::index(merged_loss))
      colnames(loss_xts) <- "RiskMetrics"
      loss_list[["RiskMetrics"]] <- loss_xts
      status[["RiskMetrics"]]$loss_nonempty <- NROW(loss_xts) > 0
      status[["RiskMetrics"]]$reason <- "ok"
    }
  }

  loss_matrix <- prepare_loss_matrix(loss_list)
included_before <- colnames(loss_matrix)
if (is.null(included_before)) included_before <- character(0)

loss_matrix_clean <- try(clean_loss_matrix(loss_matrix), silent = TRUE)
if (inherits(loss_matrix_clean, "try-error")) {
  loss_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  included_after <- character(0)
  mcs <- NULL
  mcs_models <- character(0)
} else {
  loss_matrix <- loss_matrix_clean
  included_after <- colnames(loss_matrix)
  mcs <- run_mcs(loss_matrix, alpha = config$MCS$alpha, B = config$MCS$B)
  if (is.null(mcs)) {
    mcs_models <- character(0)
  } else {
    mcs_models <- get_mcs_models(mcs)
  }
}

  if (length(status) > 0) {
    for (nm in names(status)) {
      status[[nm]]$included_in_loss_matrix <- nm %in% included_before
      status[[nm]]$included_in_mcs_input <- nm %in% included_after
      if (status[[nm]]$included_in_loss_matrix && !status[[nm]]$included_in_mcs_input) {
        status[[nm]]$reason <- "dropped_in_cleaning"
      }
      if (is.na(status[[nm]]$reason) || status[[nm]]$reason == "") {
        status[[nm]]$reason <- "other"
      }
    }
  }
  if (pipeline_label == "Pipeline B") {
    tab_dir <- if (!is.null(log_path)) dirname(log_path) else "output/tables"
    status_df <- do.call(rbind, lapply(status, function(x) {
      data.frame(
        model = x$model,
        fitted = isTRUE(x$fitted),
        rolled = isTRUE(x$rolled),
        sigma2_nonempty = isTRUE(x$sigma2_nonempty),
        loss_nonempty = isTRUE(x$loss_nonempty),
        included_in_loss_matrix = isTRUE(x$included_in_loss_matrix),
        included_in_mcs_input = isTRUE(x$included_in_mcs_input),
        reason = ifelse(is.null(x$reason), "other", x$reason),
        stringsAsFactors = FALSE
      )
    }))
    write.csv(status_df, file.path(tab_dir, "model_elimination_B.csv"), row.names = FALSE)
  }

  outlier_diag <- loss_outlier_diagnostics(loss_list, top_k = 3)

  list(
    returns = ret_xts,
    rv = rv_xts,
    S2 = sigma2_list,
    Loss = loss_list,
    Loss_matrix = loss_matrix,
    MCS = mcs,
    MCS_models = mcs_models,
    fits = fits,
    nic_list = nic_list,
    lb_table = do.call(rbind, lb_table),
    floor_diag = floor_diag,
    floor_dates = floor_dates,
    outlier_diag = outlier_diag,
    n_start = n_start
  )
}

run_pipeline_A <- function(config) {
  run_vol_contest_pipeline(
    config = config,
    ticker = config$ticker_A,
    from = config$start_A,
    to = config$end_A,
    oos_frac = config$oos_frac_A,
    refit_every = config$refit_every_A,
    refit_window = config$refit_window_A,
    include_regressors = FALSE,
    include_riskmetrics = TRUE,
    pipeline_label = "Pipeline A"
  )
}

run_pipeline_B <- function(config, regressors, log_path = NULL) {
  run_vol_contest_pipeline(
    config = config,
    ticker = config$ticker_B,
    from = config$start_B,
    to = config$end_B,
    oos_frac = config$oos_frac_B,
    refit_every = config$refit_every_B,
    refit_window = config$refit_window_B,
    include_regressors = TRUE,
    regressors = regressors,
    include_riskmetrics = TRUE,
    pipeline_label = "Pipeline B",
    log_path = log_path
  )
}

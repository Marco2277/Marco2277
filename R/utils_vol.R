# Utility functions for volatility modeling, diagnostics, and evaluation.

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required packages: %s", paste(missing, collapse = ", ")))
  }
}

fetch_returns <- function(ticker, from, to) {
  require_pkgs(c("quantmod", "xts"))
  env <- new.env()
  suppressWarnings(quantmod::getSymbols(ticker, src = "yahoo", from = from, to = to, env = env, auto.assign = TRUE))
  px <- env[[ticker]]
  adj <- try(quantmod::Ad(px), silent = TRUE)
  if (inherits(adj, "try-error") || is.null(adj)) {
    adj <- quantmod::Cl(px)
  }
  ret <- diff(log(adj))
  ret <- ret[!is.na(ret)]
  colnames(ret) <- "ret"
  ret
}

ewma_sigma2 <- function(returns, lambda = 0.94) {
  returns <- as.numeric(returns)
  n <- length(returns)
  sigma2 <- rep(NA_real_, n)
  if (n == 0) {
    return(sigma2)
  }
  sigma2[1] <- returns[1]^2
  if (n > 1) {
    for (i in 2:n) {
      sigma2[i] <- lambda * sigma2[i - 1] + (1 - lambda) * returns[i - 1]^2
    }
  }
  sigma2
}

GARCH_nll <- function(params, returns) {
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  n <- length(returns)
  sigma2 <- rep(NA_real_, n)
  sigma2[1] <- var(returns, na.rm = TRUE)
  for (i in 2:n) {
    sigma2[i] <- omega + alpha * returns[i - 1]^2 + beta * sigma2[i - 1]
  }
  -sum(dnorm(returns, mean = 0, sd = sqrt(sigma2), log = TRUE), na.rm = TRUE)
}

roll_sigma2 <- function(returns, window = 22) {
  returns <- as.numeric(returns)
  n <- length(returns)
  sigma2 <- rep(NA_real_, n)
  if (n < window) {
    return(sigma2)
  }
  for (i in seq_len(n)) {
    if (i >= window) {
      sigma2[i] <- mean(returns[(i - window + 1):i]^2, na.rm = TRUE)
    }
  }
  sigma2
}

get_std_residuals <- function(fit) {
  if (is.null(fit)) return(NULL)
  if (inherits(fit, "try-error")) return(NULL)

  z <- try(rugarch::residuals(fit, standardize = TRUE), silent = TRUE)
  if (!inherits(z, "try-error")) {
    z <- as.numeric(z)
    z <- z[is.finite(z)]
    if (length(z) > 0) return(z)
  }

  r <- try(rugarch::residuals(fit, standardize = FALSE), silent = TRUE)
  s <- try(rugarch::sigma(fit), silent = TRUE)
  if (!inherits(r, "try-error") && !inherits(s, "try-error")) {
    r <- as.numeric(r); s <- as.numeric(s)
    z2 <- r / s
    z2 <- z2[is.finite(z2)]
    if (length(z2) > 0) return(z2)
  }

  NULL
}

compute_nic <- function(fit,
                        grid = NULL,
                        n_grid = 501,
                        use_empirical_grid = TRUE,
                        p = c(0.01, 0.99),
                        fallback_grid = seq(-3, 3, length.out = 501),
                        log_path = NULL,
                        model_name = NULL) {
  if (is.null(fit)) return(NULL)

  # grid "realistico"
  if (is.null(grid)) {
    if (use_empirical_grid) {
      z <- try(get_std_residuals(fit), silent = TRUE)
      if (!inherits(z, "try-error") && !is.null(z)) {
        z <- z[is.finite(z)]
      } else {
        z <- numeric(0)
      }

      if (length(z) >= 50) {
        qs <- stats::quantile(z, probs = p, na.rm = TRUE, names = FALSE)
        if (all(is.finite(qs)) && qs[1] < qs[2]) {
          grid <- seq(qs[1], qs[2], length.out = n_grid)
        } else {
          grid <- fallback_grid
        }
      } else {
        grid <- fallback_grid
      }
    } else {
      grid <- fallback_grid
    }
  }

  nic_obj <- try(rugarch::newsimpact(fit, z = grid), silent = TRUE)
  if (inherits(nic_obj, "try-error") || is.null(nic_obj)) {
    if (!is.null(log_path)) {
      msg <- sprintf("NIC failed for %s: %s",
                     ifelse(is.null(model_name), "model", model_name),
                     ifelse(inherits(nic_obj, "try-error"), as.character(nic_obj), "NULL"))
      cat(msg, file = log_path, sep = "\n", append = TRUE)
    }
    return(NULL)
  }

  # estrazione robusta x/y
  x <- y <- NULL

  if (is.list(nic_obj) && !is.null(nic_obj$zx) && !is.null(nic_obj$zy)) {
    x <- as.numeric(nic_obj$zx); y <- as.numeric(nic_obj$zy)
  } else if (is.matrix(nic_obj) && ncol(nic_obj) >= 2) {
    x <- as.numeric(nic_obj[, 1]); y <- as.numeric(nic_obj[, 2])
  } else if (is.data.frame(nic_obj) && ncol(nic_obj) >= 2) {
    x <- as.numeric(nic_obj[[1]]); y <- as.numeric(nic_obj[[2]])
  } else if (is.numeric(nic_obj) && length(nic_obj) == length(grid)) {
    x <- as.numeric(grid); y <- as.numeric(nic_obj)
  }

  if (is.null(x) || is.null(y)) return(NULL)

  data.frame(x = x, y = y)
}

qlike_from <- function(rv, sigma2, eps = 1e-12) {
  if (length(rv) != length(sigma2)) {
    stop("rv and sigma2 must be the same length")
  }
  rv <- as.numeric(rv)
  sigma2 <- as.numeric(sigma2)
  out <- rep(NA_real_, length(rv))
  ok <- is.finite(rv) & is.finite(sigma2)
  if (!any(ok)) {
    return(out)
  }
  rv_adj <- pmax(rv[ok], eps)
  s2_adj <- pmax(sigma2[ok], eps)
  ratio <- rv_adj / s2_adj
  out[ok] <- ratio - log(ratio) - 1
  out
}

enforce_sigma2_floor <- function(sigma2, floor_value) {
  sigma2 <- as.numeric(sigma2)
  sigma2_clean <- sigma2
  sigma2_clean[!is.finite(sigma2_clean)] <- NA_real_
  sigma2_used <- ifelse(is.na(sigma2_clean), floor_value, pmax(sigma2_clean, floor_value))
  sigma2_used
}

sigma2_floor_diagnostics <- function(sigma2_raw, sigma2_used, floor_value) {
  sigma2_raw <- as.numeric(sigma2_raw)
  sigma2_used <- as.numeric(sigma2_used)

  raw_clean <- sigma2_raw
  raw_clean[!is.finite(raw_clean)] <- NA_real_

  bad <- is.na(raw_clean)
  floored <- !bad & raw_clean < floor_value

  data.frame(
    total = length(sigma2_raw),
    floor_value_used = floor_value,
    n_bad = sum(bad),
    n_floored = sum(floored),
    pct_bad = mean(bad),
    pct_floored = mean(floored),
    min_before = suppressWarnings(min(raw_clean, na.rm = TRUE)),
    min_after = suppressWarnings(min(sigma2_used, na.rm = TRUE)),
    q01_before = suppressWarnings(quantile(raw_clean, 0.01, na.rm = TRUE, names = FALSE)),
    q01_after = suppressWarnings(quantile(sigma2_used, 0.01, na.rm = TRUE, names = FALSE)),
    stringsAsFactors = FALSE
  )
}

sigma2_floor_dates <- function(idx, sigma2_raw, sigma2_used, floor_value, pipeline_label, model_name, rv_aligned = NULL) {
  idx <- as.Date(idx)
  sigma2_raw <- as.numeric(sigma2_raw)
  sigma2_used <- as.numeric(sigma2_used)

  raw_clean <- sigma2_raw
  raw_clean[!is.finite(raw_clean)] <- NA_real_

  bad <- is.na(raw_clean)
  floored <- !bad & raw_clean < floor_value
  keep <- bad | floored

  if (!any(keep)) {
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

  rv_vec <- if (!is.null(rv_aligned)) as.numeric(rv_aligned) else rep(NA_real_, length(idx))

  data.frame(
    pipeline = pipeline_label,
    model = model_name,
    date = as.character(idx[keep]),
    sigma2_raw = sigma2_raw[keep],
    sigma2_used = sigma2_used[keep],
    rv = rv_vec[keep],
    stringsAsFactors = FALSE
  )
}

prepare_loss_matrix <- function(loss_list) {
  if (length(loss_list) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  merged <- Reduce(function(x, y) merge(x, y, join = "inner"), loss_list)
  merged <- stats::na.omit(merged)
  as.matrix(merged)
}

clean_loss_matrix <- function(loss_matrix, min_rows = 15, min_cols = 2, sd_eps = 1e-10) {
  if (!is.matrix(loss_matrix)) {
    loss_matrix <- as.matrix(loss_matrix)
  }
  loss_matrix <- loss_matrix[apply(loss_matrix, 1, function(x) all(is.finite(x))), , drop = FALSE]
  if (nrow(loss_matrix) == 0) {
    stop("MCS: no finite loss rows after cleaning.")
  }
  sds <- apply(loss_matrix, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds >= sd_eps
  loss_matrix <- loss_matrix[, keep, drop = FALSE]
  if (ncol(loss_matrix) == 0) {
    stop("MCS: no models left after removing near-zero variance columns.")
  }
  rounded <- apply(loss_matrix, 2, function(x) round(x, 12))
  rounded_df <- if (is.vector(rounded)) {
    data.frame(V1 = rounded)
  } else {
    as.data.frame(rounded)
  }
  col_keys <- vapply(rounded_df, function(col) paste(col, collapse = "|"), character(1))
  dup <- duplicated(col_keys)
  loss_matrix <- loss_matrix[, !dup, drop = FALSE]
  if (ncol(loss_matrix) < min_cols || nrow(loss_matrix) < min_rows) {
    stop(sprintf("MCS: insufficient data after cleaning (rows=%d, cols=%d).", nrow(loss_matrix), ncol(loss_matrix)))
  }
  loss_matrix
}

loss_outlier_diagnostics <- function(loss_list, top_k = 3) {
  out <- list()

  for (nm in names(loss_list)) {
    x <- loss_list[[nm]]
    if (is.null(x) || NROW(x) == 0) next
    vals <- as.numeric(x[, 1])
    idx <- as.Date(zoo::index(x))
    ok <- is.finite(vals)
    if (!any(ok)) next

    ord <- order(vals[ok], decreasing = TRUE)
    k <- min(top_k, length(ord))
    top_i <- which(ok)[ord[1:k]]

    row <- data.frame(
      model = nm,
      stringsAsFactors = FALSE
    )
    for (j in 1:top_k) {
      row[[paste0("top", j, "_date")]] <- if (j <= k) as.character(idx[top_i[j]]) else NA_character_
      row[[paste0("top", j, "_qlike")]] <- if (j <= k) vals[top_i[j]] else NA_real_
    }
    out[[nm]] <- row
  }

  if (length(out) == 0) {
    return(data.frame(
      model = character(0),
      top1_date = character(0), top1_qlike = numeric(0),
      top2_date = character(0), top2_qlike = numeric(0),
      top3_date = character(0), top3_qlike = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, out)
}

run_mcs <- function(loss_matrix, alpha = 0.10, B = 1000) {
  require_pkgs(c("MCS"))
  if (!is.matrix(loss_matrix)) {
    loss_matrix <- as.matrix(loss_matrix)
  }
  loss_matrix <- loss_matrix[complete.cases(loss_matrix), , drop = FALSE]
  if (nrow(loss_matrix) == 0) {
    return(NULL)
  }
  mcs_result <- NULL
  mcs_error <- NULL
  try({
    mcs_result <- MCS::MCSprocedure(Loss = loss_matrix, alpha = alpha, B = B, statistic = "Tmax")
  }, silent = TRUE)
  if (is.null(mcs_result)) {
    try({
      mcs_result <- MCS::MCSprocedure(Loss = loss_matrix, alpha = alpha, B = B, statistic = "TR")
    }, silent = TRUE)
  }
  if (is.null(mcs_result)) {
    return(NULL)
  }
  mcs_result
}

get_mcs_models <- function(mcs_result) {
  if (is.null(mcs_result)) {
    return(character(0))
  }
  table <- NULL
  if (isS4(mcs_result)) {
    slots <- methods::slotNames(mcs_result)
    if ("show" %in% slots) {
      table <- methods::slot(mcs_result, "show")
    }
  }
  if (is.null(table)) {
    table <- try(as.data.frame(mcs_result), silent = TRUE)
    if (inherits(table, "try-error")) {
      table <- NULL
    }
  }
  if (is.null(table)) {
    return(character(0))
  }
  table <- as.data.frame(table)
  model_names <- rownames(table)
  if ("MCS_M" %in% colnames(table)) {
    return(model_names[table$MCS_M >= 1 - 1e-12])
  }
  if ("MCS_R" %in% colnames(table)) {
    return(model_names[table$MCS_R >= 1 - 1e-12])
  }
  model_names
}

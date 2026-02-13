library(quantmod)
library(xts)
library(zoo)
library(dplyr)
library(tidyr)
library(purrr)

#' Backtest a long/short strategy using forecasted mean and volatility
#' @param prices xts object of adjusted close prices (dates x tickers)
#' @param mu optional xts of expected cumulative returns over horizon h
#' @param sigma optional xts of expected cumulative standard deviations over horizon h
#' @param h holding horizon in trading days
#' @param theta_S score threshold
#' @param kappa transaction cost per side as fraction of notional
#' @param C0 initial capital
#' @return list containing equity_curve, cycle_returns, trades, stats
backtest_longshort <- function(prices, mu = NULL, sigma = NULL,
                               h = 5, theta_S = 0.5,
                               kappa = 0.0005, C0 = 100000) {
  stopifnot(inherits(prices, "xts"))

  # Compute returns if forecasts not provided (simple returns)
  returns <- diff(prices) / lag(prices)

  if (is.null(mu) || is.null(sigma)) {
    # Rendimento cumulato a h giorni (somma dei rendimenti semplici)
    mu <- rollapply(returns, width = h, FUN = function(x) apply(x, 2, sum, na.rm = TRUE),
                    by.column = FALSE, align = "right", fill = NA)
    # Volatilità del rendimento cumulato: deviazione standard giornaliera * sqrt(h)
    sigma <- rollapply(returns, width = h, FUN = function(x) {
      sd_daily <- apply(x, 2, sd, na.rm = TRUE)
      sd_daily * sqrt(h)
    }, by.column = FALSE, align = "right", fill = NA)
  }

  # Reindex mu/sigma to the price dates to prevent dimname/length mismatches
  ref_index <- index(prices)
  tickers <- colnames(prices)

  template <- xts(matrix(NA_real_, nrow = length(ref_index), ncol = length(tickers)), order.by = ref_index)
  colnames(template) <- tickers

  fill_with_template <- function(x) {
    # Keep only overlapping columns, then place them into the template aligned to price dates
    if (!is.null(colnames(x))) {
      x <- x[, intersect(colnames(x), tickers), drop = FALSE]
    } else {
      x <- x[, 0]
    }

    out <- template
    if (ncol(x) > 0) {
      out[, colnames(x)] <- x[ref_index, colnames(x)]
    }
    out
  }

  mu <- fill_with_template(mu)
  sigma <- fill_with_template(sigma)

  n_dates <- nrow(prices)
  tickers <- colnames(prices)

  # Identify valid decision rows where signals and prices exist for entry/exit
  candidate_idx <- if (n_dates > h) seq_len(n_dates - h) else integer(0)
  has_signal <- rowSums(!is.na(mu)) > 0 & rowSums(!is.na(sigma)) > 0
  entry_ok <- rowSums(!is.na(prices[candidate_idx + 1, , drop = FALSE])) > 0
  exit_ok <- rowSums(!is.na(prices[candidate_idx + h, , drop = FALSE])) > 0
  has_prices <- entry_ok & exit_ok
  valid_rows <- candidate_idx[has_signal[candidate_idx] & has_prices]

  if (length(valid_rows) == 0) stop("Not enough data to run backtest. Ensure mu/sigma windows and price data overlap.")

  start_idx <- min(valid_rows)
  decision_idx <- seq(from = start_idx, to = n_dates - h, by = h)

  capital <- C0
  equity_dates <- index(prices)[decision_idx + h]
  equity_values <- numeric(length(decision_idx))
  cycle_returns <- numeric(length(decision_idx))
  trades_list <- list()
  trade_counter <- 1

  for (step_idx in seq_along(decision_idx)) {
    t_idx <- decision_idx[step_idx]
    entry_idx <- t_idx + 1
    exit_idx <- t_idx + h

    mu_t <- as.numeric(mu[t_idx, ])
    sigma_t <- as.numeric(sigma[t_idx, ])
    names(mu_t) <- names(sigma_t) <- tickers

    score <- mu_t / sigma_t
    valid <- !is.na(score) & !is.na(mu_t) & !is.na(sigma_t) & sigma_t != 0

    # Ensure prices exist for entry and exit
    tradable <- valid & !is.na(as.numeric(prices[entry_idx, ])) & !is.na(as.numeric(prices[exit_idx, ]))
    tradable <- tradable & abs(score) >= theta_S

    I_t <- tickers[tradable]

    if (length(I_t) == 0) {
      # No trades
      cycle_returns[step_idx] <- 0
      equity_values[step_idx] <- capital
      next
    }

    mu_sel <- mu_t[I_t]
    sigma_sel <- sigma_t[I_t]
    score_sel <- score[I_t]

    d_i <- sign(mu_sel)
    w_tilde <- abs(mu_sel)
    w_i <- w_tilde / sum(w_tilde)

    pnl_vec <- numeric(length(I_t))

    for (j in seq_along(I_t)) {
      ticker <- I_t[j]
      direction <- ifelse(d_i[j] > 0, "long", "short")
      P_in <- as.numeric(prices[entry_idx, ticker])
      P_out <- as.numeric(prices[exit_idx, ticker])

      C_abs <- w_i[j] * capital
      n_shares <- floor(C_abs / P_in)

      if (n_shares <= 0) next

      pnl_gross <- d_i[j] * n_shares * (P_out - P_in)
      V_in <- n_shares * P_in
      V_out <- n_shares * P_out
      cost <- kappa * (V_in + V_out)
      pnl_net <- pnl_gross - cost
      pnl_vec[j] <- pnl_net

      trades_list[[trade_counter]] <- tibble(
        ticker = ticker,
        t_entry = index(prices)[entry_idx],
        t_exit = index(prices)[exit_idx],
        direction = direction,
        n_shares = n_shares,
        P_in = P_in,
        P_out = P_out,
        PnL_net = pnl_net
      )
      trade_counter <- trade_counter + 1
    }

    pnl_port <- sum(pnl_vec, na.rm = TRUE)
    R_port <- pnl_port / capital
    capital <- capital * (1 + R_port)

    cycle_returns[step_idx] <- R_port
    equity_values[step_idx] <- capital
  }

  equity_curve <- xts(equity_values, order.by = equity_dates)
  names(equity_curve) <- "equity"

  # Performance stats
  mean_ret <- mean(cycle_returns, na.rm = TRUE)
  sd_ret <- sd(cycle_returns, na.rm = TRUE)
  sharpe <- ifelse(sd_ret == 0, NA, mean_ret / sd_ret * sqrt(252 / h))

  cummax_equity <- cummax(as.numeric(equity_curve))
  drawdown <- as.numeric(equity_curve) / cummax_equity - 1
  max_dd <- min(drawdown, na.rm = TRUE)

  stats <- list(
    annualized_sharpe = sharpe,
    max_drawdown = max_dd,
    total_return = as.numeric(last(equity_curve)) / C0 - 1
  )

  trades_df <- if (length(trades_list) > 0) bind_rows(trades_list) else tibble()

  list(
    equity_curve = equity_curve,
    cycle_returns = xts(cycle_returns, order.by = equity_dates),
    trades = trades_df,
    stats = stats
  )
}

# ----------------------------
# Example usage
if (sys.nframe() == 0) {
  tickers <- c(
    "AAPL", "MSFT", "AMZN", "GOOGL", "META", "NVDA", "TSLA",
    "JPM", "BAC", "C", "WFC", "GS", "MS", "V", "MA",
    "JNJ", "PFE", "MRK", "ABBV", "UNH", "LLY", "TMO",
    "XOM", "CVX", "COP", "SLB",
    "KO", "PEP", "MCD", "WMT", "HD", "COST", "NKE", "DIS", "NFLX",
    "T", "VZ", "IBM", "CSCO", "ORCL", "INTC", "TXN", "QCOM",
    "HON", "UPS", "CAT", "BA", "UNP", "LMT", "ADP"
  )
  start_date <- "2014-01-01"
  end_date <- "2024-12-31"

  getSymbols(tickers, src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE, warnings = FALSE)
  price_list <- lapply(tickers, function(tk) Ad(get(tk)))
  prices <- do.call(merge, price_list)
  colnames(prices) <- tickers

  # Align on common dates (drop leading NAs)
  prices <- na.omit(prices)

  res <- backtest_longshort(prices, h = 5, theta_S = 0.5, kappa = 0.0005, C0 = 100000)

  cat("\nPrime date dell'equity curve (uscita dai cicli):\n")
  print(head(res$equity_curve))

  cat("\nUltime date dell'equity curve (uscita dai cicli):\n")
  print(tail(res$equity_curve))

  cat("\nStatistiche riassuntive:\n")
  print(res$stats)

  cat("\nCapitale finale: ", formatC(as.numeric(last(res$equity_curve)), format = "f", digits = 2), "\n", sep = "")

  cat("\nOperazioni eseguite (tutte le righe del backtest):\n")
  if (nrow(res$trades) > 0) {
    print(res$trades, n = nrow(res$trades))
  } else {
    cat("Nessuna operazione eseguita.\n")
  }

  if (NROW(res$equity_curve) > 1) {
    plot(res$equity_curve, main = "Equity curve della strategia", ylab = "Capitale", major.ticks = "years", grid.ticks.on = "years")
  }
}

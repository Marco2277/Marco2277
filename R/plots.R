# Plot builders and helpers.

make_sigma_panel <- function(sigma2_list, title, floor_value, compact = FALSE) {
  op <- par(no.readonly = TRUE)
  op$new <- FALSE
  on.exit(try(par(op), silent = TRUE), add = TRUE)

  n <- length(sigma2_list)
  if (n == 0) {
    plot.new()
    title(main = title)
    return(invisible(NULL))
  }

  nrow <- ceiling(sqrt(n))
  ncol <- ceiling(n / nrow)

  # Default = come prima (per PNG). Compact = per RStudio Plots (evita "figure margins too large").
  if (isTRUE(compact)) {
    mar <- c(3.2, 3.4, 1.6, 0.6)
    oma <- c(1.4, 1.4, 2.0, 0.2)
    mgp <- c(2.0, 0.55, 0)
    cex_axis <- 0.65
    cex_lab  <- 0.75
    cex_main <- 0.80
    title_cex <- 1.0
    title_line <- 0.5
    xlab_line <- 0.6
    ylab_line <- 0.6
  } else {
    mar <- c(5.2, 4.6, 2.0, 0.8)
    oma <- c(2.2, 2.2, 3.4, 0.2)
    mgp <- c(2.3, 0.75, 0)
    cex_axis <- 0.75
    cex_lab  <- 0.85
    cex_main <- 0.85
    title_cex <- 1.0
    title_line <- 0.8
    xlab_line <- 0.8
    ylab_line <- 0.9
  }

  par(
    mfrow = c(nrow, ncol),
    mar = mar,
    oma = oma,
    mgp = mgp,
    las = 1,
    tcl = -0.2,
    xpd = FALSE,
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main
  )

  for (nm in names(sigma2_list)) {
    s2 <- sigma2_list[[nm]]
    if (is.null(s2) || NROW(s2) < 1) {
      plot.new(); title(main = nm); next
    }

    x <- zoo::index(s2)
    y <- as.numeric(s2)
    ok <- !is.na(x) & is.finite(y) & y > 0
    if (!any(ok)) {
      plot.new(); title(main = nm); next
    }

    x_ok <- as.Date(x[ok])
    y_ok <- y[ok]

    log_y <- log10(y_ok)
    log_y <- log_y[is.finite(log_y)]
    if (length(log_y) < 2) {
      plot.new(); title(main = nm); next
    }

    qs <- stats::quantile(log_y, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    if (!all(is.finite(qs)) || qs[1] == qs[2]) qs <- range(log_y, na.rm = TRUE)

    y_clip <- 10^pmax(pmin(log10(y_ok), qs[2]), qs[1])

    plot(x_ok, y_clip, type = "l", main = nm, xlab = "", ylab = "", log = "y", xaxt = "n")
    at <- pretty(x_ok, n = 4)
    axis.Date(1, at = at, format = "%b %y", las = 2)
    abline(h = floor_value, lty = 2)
  }

  mtext("Date",    side = 1, outer = TRUE, line = xlab_line, cex = cex_lab)
  mtext(title,     side = 3, outer = TRUE, line = title_line, cex = title_cex)

  invisible(NULL)
}

plot_nic <- function(nic_list, title,
                     normalize_at_zero = TRUE,
                     log10_scale = TRUE,
                     clip_probs = c(0.01, 0.99),
                     eps = 1e-12,
                     show_global_ylab = FALSE,
                     min_ylim_span = 1e-3) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  n <- length(nic_list)
  if (n == 0) {
    plot.new(); title(main = title)
    return(invisible(NULL))
  }

  nrow <- ceiling(sqrt(n))
  ncol <- ceiling(n / nrow)

  global_ylab <- if (normalize_at_zero && log10_scale) "log10(rel. response)" else "Response"

  par(
    mfrow = c(nrow, ncol),
    mar = c(3.4, 4.4, 1.8, 0.8),
    oma = c(0, if (isTRUE(show_global_ylab)) 2.8 else 0, 2.2, 0),
    mgp = c(2.3, 0.75, 0),
    las = 1,
    tcl = -0.2,
    cex.axis = 0.85,
    cex.lab  = 0.90,
    cex.main = 0.95,
    xpd = FALSE
  )

  for (nm in names(nic_list)) {
    nic <- nic_list[[nm]]

    if (is.null(nic) || !all(c("x", "y") %in% names(nic))) {
      plot.new(); title(main = nm)
      text(0.5, 0.5, "NIC not available")
      next
    }

    x <- as.numeric(nic$x)
    y <- as.numeric(nic$y)
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) {
      plot.new(); title(main = nm)
      text(0.5, 0.5, "NIC not available")
      next
    }

    x_ok <- x[ok]
    y_ok <- y[ok]

    # Heuristica originale: se ci sono valori <=0, assumo che y sia su scala log e torno su scala positiva.
    if (any(y_ok <= 0, na.rm = TRUE)) {
      y_ok <- exp(y_ok)
    }

    if (normalize_at_zero) {
      i0 <- which.min(abs(x_ok))
      base <- max(y_ok[i0], eps)
      y_ok <- y_ok / base
    }

    if (log10_scale) {
      y_ok <- log10(pmax(y_ok, eps))
    }

    ord <- order(x_ok)
    x_ok <- x_ok[ord]
    y_ok <- y_ok[ord]

    qs <- stats::quantile(y_ok, probs = clip_probs, na.rm = TRUE, names = FALSE)
    if (!all(is.finite(qs)) || qs[1] == qs[2]) qs <- range(y_ok, na.rm = TRUE)
    y_clip <- pmax(pmin(y_ok, qs[2]), qs[1])

    # Stabilizza l’asse Y quando la curva è quasi piatta (evita scale ridicole tipo e-07 e “spike” visivi).
    ylim <- range(y_clip, na.rm = TRUE)
    span <- diff(ylim)
    if (!is.finite(span) || span < min_ylim_span) {
      mid <- mean(ylim)
      pad <- max(min_ylim_span, 0.05 * (abs(mid) + 1e-6))
      ylim <- c(mid - pad, mid + pad)
    } else {
      pad <- 0.04 * span
      ylim <- c(ylim[1] - pad, ylim[2] + pad)
    }

    plot(x_ok, y_clip, type = "n", main = nm, xlab = "Shock", ylab = "", ylim = ylim)

    if (normalize_at_zero) {
      abline(h = if (log10_scale) 0 else 1, col = "gray70")
    } else {
      abline(h = 0, col = "gray70")
    }
    lines(x_ok, y_clip, lwd = 1)
  }

  if (isTRUE(show_global_ylab)) {
    mtext(global_ylab, side = 2, outer = TRUE, line = 1.0, cex = 0.95)
  }
  mtext(title, side = 3, outer = TRUE, line = 0.5, cex = 1.1)
}

save_qlike_boxplot <- function(loss_matrix, file_path, title,
                               mar = c(8, 4.4, 4.0, 2.0) + 0.1,
                               cex_axis = 1.0) {
  save_png(file_path, 1800, 1000, function() {
    op <- par(no.readonly = TRUE)
    op$new <- FALSE
    on.exit(try(par(op), silent = TRUE), add = TRUE)

    if (is.null(loss_matrix) || length(loss_matrix) == 0 ||
        !is.matrix(loss_matrix) || ncol(loss_matrix) < 1 || nrow(loss_matrix) < 1) {
      plot.new()
      title(main = title)
      text(0.5, 0.5, "No data to plot")
      return(invisible(NULL))
    }

    par(mar = mar, las = 1, tcl = -0.2, xpd = NA, cex.axis = cex_axis)
    boxplot(loss_matrix, las = 2, main = title)
  })
  invisible(NULL)
}

save_sigma_panel <- function(sigma2_list, file_path, title, floor_value) {
  save_png(file_path, 1800, 1000, function() {
    make_sigma_panel(sigma2_list, title, floor_value)
  })
  invisible(NULL)
}

save_nic_panel <- function(nic_list, file_path, title) {
  save_png(file_path, 1800, 1000, function() {
    plot_nic(nic_list, title)
  })
  invisible(NULL)
}

save_and_show_png <- function(file_path, width, height, plot_fun, open_new_device = FALSE) {
  save_png(file_path, width, height, plot_fun)
  run_plot_safely(plot_fun)
  invisible(NULL)
}

#' Create Cross-Correlogram Plot Matrix
#'
#' @param data A data frame or matrix where each column is a time series
#' @param max_lag Maximum lag to display (if NULL, uses length/4)
#' @param ci Confidence interval level (0.95 for 95% CI)
#' @param main_title Overall title for the plot
#' @param cex.main Size of the main title
#' @param mar Margins for individual plots
#' @param oma Outer margins for the entire figure
#' @return NULL (creates plot as side effect)
#' @export
cross_correlogram <- function(data, max_lag = NULL, ci = 0.95,
                              main_title = NULL, cex.main = 1.2,
                              mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 2)) {
  # Convert to matrix if data frame
  if (is.data.frame(data)) {
    data_matrix <- as.matrix(data)
  } else {
    data_matrix <- data
  }

  # Get dimensions and names
  n_series <- ncol(data_matrix)
  n_obs <- nrow(data_matrix)

  # Get column names
  col_names <- colnames(data_matrix)
  if (is.null(col_names)) {
    col_names <- paste0("Y", 1:n_series)
  }

  # Set default max_lag
  if (is.null(max_lag)) {
    max_lag <- floor(n_obs / 4)
  }

  # Calculate confidence interval bounds
  ci_bound <- qnorm((1 + ci) / 2) / sqrt(n_obs)

  # Set up the plotting layout
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_series, n_series),
      mar = mar,
      oma = oma)

  # Create each subplot
  for (i in 1:n_series) {
    for (j in 1:n_series) {

      if (i == j) {
        # Diagonal: Auto-correlation
        series <- data_matrix[, i]
        valid_idx <- !is.na(series)

        if (sum(valid_idx) > 3) {
          series_clean <- series[valid_idx]
          acf_result <- acf(series_clean, lag.max = max_lag,
                            plot = FALSE, na.action = na.pass)

          # Plot ACF
          plot(acf_result$lag, acf_result$acf, type = "h",
               xlim = c(0, max_lag), ylim = c(-1, 1),
               xlab = "", ylab = "", main = "",
               col = "blue", lwd = 1.5)

          # Add confidence intervals
          abline(h = c(ci_bound, -ci_bound), col = "red", lty = 2, lwd = 1)
          abline(h = 0, col = "black", lwd = 1)

          # Add labels for diagonal
          if (i == 1) {
            title(main = col_names[i], cex.main = 0.8, line = 0.5)
          }
          if (i == n_series) {
            mtext(col_names[i], side = 1, line = 2, cex = 0.7)
          }
          if (j == 1) {
            mtext(col_names[i], side = 2, line = 2, cex = 0.7)
          }

        } else {
          # Empty plot if insufficient data
          plot.new()
          text(0.5, 0.5, "Insufficient\nData", cex = 0.8)
        }

      } else {
        if (i>j){
          d_sign = -1
          u_sign = 0
        } else{
          d_sign = 0
          u_sign = 1

        }
        # Off-diagonal: Cross-correlation
        series1 <- data_matrix[, i]
        series2 <- data_matrix[, j]

        # Remove missing values
        valid_idx <- complete.cases(series1, series2)

        if (sum(valid_idx) > 3) {
          series1_clean <- series1[valid_idx]
          series2_clean <- series2[valid_idx]

          ccf_result <- stats::ccf(series1_clean, series2_clean,
                            lag.max = max_lag, plot = FALSE,
                            na.action = na.pass)

          # Plot CCF
          plot(ccf_result$lag, ccf_result$acf, type = "h",
               xlim = c(d_sign * max_lag, u_sign* max_lag), ylim = c(-1, 1),
               xlab = "", ylab = "", main = "",
               col = "blue", lwd = 1.5)

          # Add confidence intervals
          abline(h = c(ci_bound, -ci_bound), col = "red", lty = 2, lwd = 1)
          abline(h = 0, col = "black", lwd = 1)

          # Add column labels on top row
          if (i == 1) {
            title(main = col_names[j], cex.main = 0.8, line = 0.5)
          }

          # Add row labels on left column
          if (j == 1) {
            mtext(col_names[i], side = 2, line = 2, cex = 0.7)
          }

          # Add x-axis labels on bottom row
          if (i == n_series) {
            mtext("Lag", side = 1, line = 2, cex = 0.7)
          }

        } else {
          # Empty plot if insufficient data
          plot.new()
          text(0.5, 0.5, "Insufficient\nData", cex = 0.8)
        }
      }
    }
  }
  # Add overall title
  if (is.null(main_title)) {
    main_title <- paste("Cross correlogram of the", n_series, "time series")
  }

  mtext(main_title, outer = TRUE, cex = cex.main, line = 2)

  # Add axis labels
  mtext("Lag", side = 1, outer = TRUE, cex = 1, line = 2)
  mtext("Cross-Correlation", side = 2, outer = TRUE, cex = 1, line = 2)
}


#' Plot Maximum Cross-Correlation Across All Series Pairs
#'
#' For each lag, computes cross-correlation between all pairs of series
#' and plots the maximum absolute correlation value
#'
#' @param data Data frame or matrix where each column is a time series
#' @param max_lag Maximum lag to compute and plot
#' @param main Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param col Color for points/lines
#' @param pch Point character
#' @param lwd Line width
#' @param cex Point size
#' @param add_line Whether to connect points with lines
#' @param grid Whether to add grid
#' @param cex.lab Label text size
#' @param cex.axis Axis text size
#' @param cex.main Title text size
#' @return List with lag values and maximum correlations
plot_max_ccf_all_pairs <- function(data, max_lag = 20,
                                   main = "Maximum Cross-Correlation Across All Pairs",
                                   xlab = "Lag", ylab = "Maximum Cross-Correlation",
                                   col = "blue", pch = 19, lwd = 2, cex = 1,
                                   add_line = TRUE, grid = TRUE,
                                   cex.lab = 1.2, cex.axis = 1, cex.main = 1.2) {

  # Convert to matrix
  if (is.data.frame(data)) {
    data_matrix <- as.matrix(data)
  } else {
    data_matrix <- data
  }

  n_series <- ncol(data_matrix)
  n_obs <- nrow(data_matrix)

  if (n_series < 2) {
    stop("Need at least 2 time series")
  }

  # Ensure max_lag doesn't exceed series length
  max_lag <- min(max_lag, floor(n_obs / 2))

  # Initialize arrays to store maximum correlations at each lag
  lags <- 0:max_lag
  max_correlations <- numeric(length(lags))

  cat("Computing cross-correlations for", choose(n_series, 2), "pairs...\n")

  # For each lag, find maximum correlation across all pairs
  for (lag_idx in seq_along(lags)) {
    current_lag <- lags[lag_idx]
    max_corr_at_lag <- 0

    # Compare all pairs of series
    for (i in 1:(n_series-1)) {
      for (j in (i+1):n_series) {

        series1 <- data_matrix[, i]
        series2 <- data_matrix[, j]

        # Remove missing values
        valid_idx <- complete.cases(series1, series2)
        if (sum(valid_idx) < 10) next

        series1_clean <- series1[valid_idx]
        series2_clean <- series2[valid_idx]

        # Calculate cross-correlation for this pair
        ccf_result <- stats::ccf(series1_clean, series2_clean,
                          lag.max = max_lag, plot = FALSE)

        # Find correlation at current lag
        lag_position <- which(ccf_result$lag == current_lag)
        if (length(lag_position) > 0) {
          corr_value <- abs(ccf_result$acf[lag_position])
          max_corr_at_lag <- max(max_corr_at_lag, corr_value, na.rm = TRUE)
        }
      }
    }

    max_correlations[lag_idx] <- max_corr_at_lag
  }

  # Create the plot
  ylim <- c(0, min(1, max(max_correlations, na.rm = TRUE) * 1.1))

  if (add_line) {
    plot(lags, max_correlations,
         type = "b", pch = pch, col = col, cex = cex, lwd = lwd,
         main = main, xlab = xlab, ylab = ylab,
         xlim = c(0, max_lag), ylim = ylim,
         cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
  } else {
    plot(lags, max_correlations,
         type = "p", pch = pch, col = col, cex = cex,
         main = main, xlab = xlab, ylab = ylab,
         xlim = c(0, max_lag), ylim = ylim,
         cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
  }

  # Add grid
  if (grid) {
    grid(col = "lightgray", lty = 3)
  }

  # Add some summary statistics
  cat("Maximum correlation at lag 0:", round(max_correlations[1], 3), "\n")
  cat("Overall maximum correlation:", round(max(max_correlations, na.rm = TRUE), 3),
      "at lag", lags[which.max(max_correlations)], "\n")

  # Return results
  result <- list(
    lag = lags,
    max_correlation = max_correlations,
    summary = data.frame(
      lag = lags,
      max_corr = max_correlations
    )
  )

  invisible(result)
}


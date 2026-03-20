# ============================================================
# File: 40_Sup_smoothing.R
# Purpose:
#   LOESS smoothing of bin-wise odds ratios for the
#   superiority comparison (Treatment I vs Control C).
#
# Contents:
#   1. Smooth_OR_Sup_tune()
#      Computes raw odds ratios within each Decipher bin,
#      applies continuity correction when needed, and fits
#      a LOESS smoother on the log-odds-ratio scale.
#
# Notes:
#   The smoothed odds ratios are used in Stage 2 as a
#   gatekeeping step before likelihood ratio testing of
#   expanded thresholds.
# ============================================================



Smooth_OR_Sup_tune <- function(data, 
                               span = 0.85, 
                               combine_bins = c(8, 9, 10),
                               combined_x = 8,
                               plot = TRUE) {
  
  # Create analysis bin: combine upper bins into one pooled subgroup
  data$analysis_bin <- ifelse(data$decipher_bin %in% combine_bins,
                              paste0(min(combine_bins), "+"),
                              as.character(data$decipher_bin))
  
  original_bins <- sort(unique(data$decipher_bin))
  kept_bins <- original_bins[!original_bins %in% combine_bins]
  analysis_bins <- c(as.character(kept_bins), paste0(min(combine_bins), "+"))
  bin_x <- c(kept_bins, combined_x)
  
  raw_ORs <- numeric(length(analysis_bins))
  eC_vec <- neC_vec <- nC_vec <- numeric(length(analysis_bins))
  eI_vec <- neI_vec <- nI_vec <- numeric(length(analysis_bins))
  event_rate_C <- event_rate_I <- rep(NA_real_, length(analysis_bins))
  
  calc_or <- function(eC, neC, eI, neI) {
    cells <- c(eC, neC, eI, neI)
    if (any(cells == 0)) {
      eC  <- eC  + 0.5
      neC <- neC + 0.5
      eI  <- eI  + 0.5
      neI <- neI + 0.5
    }
    (eI / neI) / (eC / neC)
  }
  
  for (i in seq_along(analysis_bins)) {
    this_bin <- analysis_bins[i]
    
    if (this_bin == paste0(min(combine_bins), "+")) {
      bin_data <- subset(data, decipher_bin %in% combine_bins)
    } else {
      bin_data <- subset(data, decipher_bin == as.numeric(this_bin))
    }
    
    eC <- sum(bin_data$treatment == "C" & bin_data$status_5yr == 1L, na.rm = TRUE)
    nC <- sum(bin_data$treatment == "C", na.rm = TRUE)
    eI <- sum(bin_data$treatment == "I" & bin_data$status_5yr == 1L, na.rm = TRUE)
    nI <- sum(bin_data$treatment == "I", na.rm = TRUE)
    
    neC <- nC - eC
    neI <- nI - eI
    
    eC_vec[i] <- eC
    neC_vec[i] <- neC
    nC_vec[i] <- nC
    eI_vec[i] <- eI
    neI_vec[i] <- neI
    nI_vec[i] <- nI
    
    event_rate_C[i] <- ifelse(nC > 0, eC / nC, NA_real_)
    event_rate_I[i] <- ifelse(nI > 0, eI / nI, NA_real_)
    
    if (nC > 0 && nI > 0) {
      raw_ORs[i] <- calc_or(eC, neC, eI, neI)
    } else {
      raw_ORs[i] <- NA_real_
    }
  }
  
  # Fit LOESS to all available analysis-level ORs
  valid_idx <- !is.na(raw_ORs)
  smoothed_ORs <- rep(NA_real_, length(analysis_bins))
  
  if (sum(valid_idx) >= 2) {
    suppressWarnings({
      loess_fit <- loess(log(raw_ORs[valid_idx]) ~ bin_x[valid_idx],
                         span = span,
                         degree = 1)
      smoothed_pred <- predict(loess_fit, newdata = bin_x)
      smoothed_ORs <- exp(smoothed_pred)
      
      if (anyNA(smoothed_ORs)) {
        fitted_valid <- exp(predict(loess_fit, newdata = bin_x[valid_idx]))
        smoothed_ORs <- approx(x = bin_x[valid_idx],
                               y = fitted_valid,
                               xout = bin_x,
                               rule = 2)$y
      }
    })
  } else {
    smoothed_ORs <- raw_ORs
  }
  
  if (plot) {
    ymin <- min(c(raw_ORs, smoothed_ORs), na.rm = TRUE)
    ymax <- max(c(raw_ORs, smoothed_ORs), na.rm = TRUE)
    
    plot(bin_x, raw_ORs, type = "b", pch = 16, col = "blue", lwd = 2,
         ylim = c(ymin - 0.05, ymax + 0.05),
         xaxt = "n",
         xlab = "Decipher Bin",
         ylab = "Odds Ratio (I vs C)",
         main = "Raw vs LOESS-Smoothed Odds Ratios")
    
    axis(1, at = bin_x, labels = analysis_bins)
    
    lines(bin_x, smoothed_ORs, type = "b", pch = 17, col = "darkorange",
          lty = 2, lwd = 2)
    abline(h = 1, col = "gray", lty = 3)
    
    legend("topright",
           legend = c("Raw ORs", "LOESS Smoothed ORs", "OR = 1"),
           col = c("blue", "darkorange", "gray"),
           pch = c(16, 17, NA),
           lty = c(1, 2, 3),
           lwd = 2)
  }
  
  data.frame(
    bin = analysis_bins,
    x_for_smoothing = bin_x,
    eC = eC_vec,
    neC = neC_vec,
    nC = nC_vec,
    event_rate_C = event_rate_C,
    eI = eI_vec,
    neI = neI_vec,
    nI = nI_vec,
    event_rate_I = event_rate_I,
    raw_OR = raw_ORs,
    smoothed_OR = round(smoothed_ORs, 2)
  )
}
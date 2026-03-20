# ============================================================
# File: 30_NI_smoothing.R
# Purpose:
#   LOESS smoothing of bin-wise odds ratios for the
#   non-inferiority comparison (Treatment D vs Control C).
#
# Contents:
#   1. Smooth_OR_NI_Sep()
#      Computes raw odds ratios within each Decipher bin,
#      applies continuity correction when needed, and fits
#      a LOESS smoother on the log-odds-ratio scale.
#
# Notes:
#   The smoothed odds ratios are used in Stage 2 as a
#   gatekeeping step before likelihood ratio testing of
#   expanded thresholds.
# ============================================================


Smooth_OR_NI_tune <- function(data, 
                              span = 0.85, 
                              combine_bins = c(1, 2, 3),
                              combined_x = 3,
                              plot = TRUE) {
  
  # Create analysis bin: combine lower bins into one pooled subgroup
  data$analysis_bin <- ifelse(data$decipher_bin %in% combine_bins,
                              paste0(min(combine_bins), "-", max(combine_bins)),
                              as.character(data$decipher_bin))
  
  original_bins <- sort(unique(data$decipher_bin))
  kept_bins <- original_bins[!original_bins %in% combine_bins]
  analysis_bins <- c(paste0(min(combine_bins), "-", max(combine_bins)),
                     as.character(kept_bins))
  bin_x <- c(combined_x, kept_bins)
  
  raw_ORs <- numeric(length(analysis_bins))
  eC_vec <- neC_vec <- nC_vec <- numeric(length(analysis_bins))
  eD_vec <- neD_vec <- nD_vec <- numeric(length(analysis_bins))
  event_rate_C <- event_rate_D <- rep(NA_real_, length(analysis_bins))
  
  calc_or <- function(eC, neC, eD, neD) {
    cells <- c(eC, neC, eD, neD)
    if (any(cells == 0)) {
      eC  <- eC  + 0.5
      neC <- neC + 0.5
      eD  <- eD  + 0.5
      neD <- neD + 0.5
    }
    (eD / neD) / (eC / neC)
  }
  
  for (i in seq_along(analysis_bins)) {
    this_bin <- analysis_bins[i]
    
    if (this_bin == paste0(min(combine_bins), "-", max(combine_bins))) {
      bin_data <- subset(data, decipher_bin %in% combine_bins)
    } else {
      bin_data <- subset(data, decipher_bin == as.numeric(this_bin))
    }
    
    eC <- sum(bin_data$treatment == "C" & bin_data$status_5yr == 1L, na.rm = TRUE)
    nC <- sum(bin_data$treatment == "C", na.rm = TRUE)
    eD <- sum(bin_data$treatment == "D" & bin_data$status_5yr == 1L, na.rm = TRUE)
    nD <- sum(bin_data$treatment == "D", na.rm = TRUE)
    
    neC <- nC - eC
    neD <- nD - eD
    
    eC_vec[i] <- eC
    neC_vec[i] <- neC
    nC_vec[i] <- nC
    eD_vec[i] <- eD
    neD_vec[i] <- neD
    nD_vec[i] <- nD
    
    event_rate_C[i] <- ifelse(nC > 0, eC / nC, NA_real_)
    event_rate_D[i] <- ifelse(nD > 0, eD / nD, NA_real_)
    
    if (nC > 0 && nD > 0) {
      raw_ORs[i] <- calc_or(eC, neC, eD, neD)
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
  
  # Expand back to original bins so bins 1,2,3 share the same pooled value
  smoothed_OR_full <- rep(NA_real_, length(original_bins))
  raw_OR_full      <- rep(NA_real_, length(original_bins))
  
  pooled_label <- paste0(min(combine_bins), "-", max(combine_bins))
  pooled_idx   <- which(analysis_bins == pooled_label)
  
  smoothed_OR_full[original_bins %in% combine_bins] <- smoothed_ORs[pooled_idx]
  raw_OR_full[original_bins %in% combine_bins]      <- raw_ORs[pooled_idx]
  
  for (b in kept_bins) {
    a_idx <- which(analysis_bins == as.character(b))
    o_idx <- which(original_bins == b)
    smoothed_OR_full[o_idx] <- smoothed_ORs[a_idx]
    raw_OR_full[o_idx]      <- raw_ORs[a_idx]
  }
  
  if (plot) {
    ymin <- min(c(raw_OR_full, smoothed_OR_full), na.rm = TRUE)
    ymax <- max(c(raw_OR_full, smoothed_OR_full), na.rm = TRUE)
    
    plot(original_bins, raw_OR_full, type = "b", pch = 16, col = "blue", lwd = 2,
         ylim = c(ymin - 0.05, ymax + 0.05),
         xlab = "Decipher Bin",
         ylab = "Odds Ratio (D vs C)",
         main = "Raw vs LOESS-Smoothed Odds Ratios")
    
    lines(original_bins, smoothed_OR_full, type = "b", pch = 17, col = "darkorange",
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
    bin = original_bins,
    x_for_smoothing = original_bins,
    raw_OR = raw_OR_full,
    smoothed_OR = round(smoothed_OR_full, 2)
  )
}

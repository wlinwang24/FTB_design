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
                               min_events = 1, 
                               span = 0.85, 
                               max_rawOR = 3,
                               plot = TRUE) {
  bin_ids <- sort(unique(data$decipher_bin))
  bin_mids <- bin_ids
  raw_ORs <- numeric(length(bin_ids))
  used_in_smoothing <- logical(length(bin_ids))
  
  # Initialize vectors
  eC_vec <- neC_vec <- nC_vec <- numeric(length(bin_ids))
  eI_vec <- neI_vec <- nI_vec <- numeric(length(bin_ids))
  event_rate_C <- event_rate_I <- rep(NA, length(bin_ids))
  
  for (i in seq_along(bin_ids)) {
    bin_data <- subset(data, decipher_bin == bin_ids[i])
    
    # Counts
    eC <- sum(bin_data$treatment == "C" & bin_data$status_5yr == 1L)
    nC <- sum(bin_data$treatment == "C")
    eI <- sum(bin_data$treatment == "I" & bin_data$status_5yr == 1L)
    nI <- sum(bin_data$treatment == "I")
    
    neC <- nC - eC
    neI <- nI - eI
    
    # Store counts
    eC_vec[i] <- eC; neC_vec[i] <- neC; nC_vec[i] <- nC
    eI_vec[i] <- eI; neI_vec[i] <- neI; nI_vec[i] <- nI
    event_rate_C[i] <- ifelse(nC > 0, eC / nC, NA)
    event_rate_I[i] <- ifelse(nI > 0, eI / nI, NA)
    
    # --- Continuity correction: If any cell is zero, add 0.5 to ALL ---
    cells <- c(eC, neC, eI, neI)
    if (any(cells == 0)) {
      adj_eC  <- eC  + 0.5
      adj_neC <- neC + 0.5
      adj_eI  <- eI  + 0.5
      adj_neI <- neI + 0.5
    } else {
      adj_eC  <- eC
      adj_neC <- neC
      adj_eI  <- eI
      adj_neI <- neI
    }
    raw_ORs[i] <- (adj_eI / adj_neI) / (adj_eC / adj_neC)
    
    # --- Decide if this bin is used in smoothing ---
    if (eC >= min_events && eI >= min_events && raw_ORs[i] < max_rawOR) {
      used_in_smoothing[i] <- TRUE
    } else {
      used_in_smoothing[i] <- FALSE
    }
  }
  
  # --- Fit LOESS on eligible bins ---
  valid_idx <- !is.na(raw_ORs) & used_in_smoothing
  smoothed_ORs <- rep(NA, length(bin_ids))
  
  if (sum(valid_idx) >= 3) {
    suppressWarnings({
      loess_fit <- loess(log(raw_ORs[valid_idx]) ~ bin_mids[valid_idx], span = span)
      smoothed_pred <- predict(loess_fit, newdata = bin_mids)
      
      # Back-transform
      smoothed_ORs <- exp(smoothed_pred)
      
      # If prediction has NAs (e.g., edge bins), linearly extrapolate
      if (anyNA(smoothed_ORs)) {
        smoothed_ORs <- approx(
          x = bin_mids[valid_idx],
          y = exp(predict(loess_fit, newdata = bin_mids[valid_idx])),
          xout = bin_mids,
          rule = 2  # linear extrapolation at edges
        )$y
      }
    })
  } else {
    smoothed_ORs <- raw_ORs  # fallback
  }
  
  # --- Plot ---
  if (plot) {
    plot(bin_mids, raw_ORs, type = "b", pch = 16, col = "blue", lwd = 2,
         ylim = c(min(raw_ORs, smoothed_ORs, na.rm = TRUE) - 0.05,
                  max(raw_ORs, smoothed_ORs, na.rm = TRUE) + 0.05),
         xlab = "Decipher Bin", ylab = "Odds Ratio (I vs C)",
         main = "Raw vs LOESS-Smoothed Odds Ratios")
    lines(bin_mids, smoothed_ORs, type = "b", pch = 17, col = "darkorange", lty = 2, lwd = 2)
    abline(h = 1, col = "gray", lty = 3)
    legend("topright", legend = c("Raw ORs", "LOESS Smoothed ORs", "OR = 1"),
           col = c("blue", "darkorange", "gray"), pch = c(16, 17, NA),
           lty = c(1, 2, 3), lwd = 2)
  }
  
  # --- Return full dataset ---
  return(data.frame(
    bin = bin_ids,
    eC = eC_vec, neC = neC_vec, nC = nC_vec, event_rate_C = event_rate_C,
    eI = eI_vec, neI = neI_vec, nI = nI_vec, event_rate_I = event_rate_I,
    raw_OR = raw_ORs,
    smoothed_OR = round(smoothed_ORs, 2),
    used_in_smoothing = used_in_smoothing
  ))
}











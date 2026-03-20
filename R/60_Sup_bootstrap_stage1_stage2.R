# ============================================================
# File: 60_Sup_bootstrap_stage1_stage2.R
# Purpose:
#   Two-stage superiority threshold selection procedure
#   for the FTB design.
#
# Contents:
#   1. c2_bootstrap_LRT_Smooth()
#      Implements:
#      - Stage 1: bootstrap-based threshold selection using
#        the minimum re-centered Wald-type Z-statistic across
#        candidate c2 values
#      - Stage 2: fixed-order sequential expansion using
#        smoothed odds ratios for gatekeeping and one-sided
#        likelihood ratio tests at the nominal level
#
# Notes:
#   The procedure first identifies the most favorable c2 and
#   performs a global test to control family-wise error across
#   the initial threshold search. If significant, the selected
#   threshold is expanded sequentially to smaller c2 values.
# ============================================================


c2_bootstrap_LRT_Smooth <- function(data,
                                    c2_values,
                                    alpha = 0.025,
                                    n_bootstrap = 1000) {
  # 1) Get observed Z-stats for each c2
  observed_results <- lapply(c2_values, function(c2) {
    res <- analyze_c2_zstat_ordinal(data, c2)
    if (is.null(res) || is.na(res$z_stat)) {
      return(list(c2     = c2,
                  z_stat = Inf,  # store Inf so it won't pass
                  p_value = NA,
                  or = NA,
                  log_or = NA,
                  se_log_or = NA))
    } else {
      return(res)
    }
  })
  
  # If everything is NA or null, bail out
  if (all(sapply(observed_results, function(x) is.na(x$z_stat)))) {
    stop("No valid results for any c2 value.")
  }
  
  # Extract the observed log(OR) and SE
  log_ors    <- sapply(observed_results, function(x) x$log_or)
  se_log_ors <- sapply(observed_results, function(x) x$se_log_or)
  
  # Observed Z-stats (superiority):
  # H0: log(OR) >= 0 => Z = log(OR)/SE
  observed_z_stats <- log_ors / se_log_ors
  
  
  # Order them from smallest to largest
  original_order <- order(observed_z_stats, decreasing = FALSE)
  sorted_observed_z_stats <- observed_z_stats[original_order]
  
  n_cutoffs <- length(c2_values)
  smoothed_OR <- Smooth_OR_Sup_tune(data = data, span = 1, plot = FALSE)$smoothed_OR
  
  
  bootstrap_studentized_z <- matrix(NA, nrow = n_bootstrap, ncol = n_cutoffs)
  
  # 2) For the bootstrap, we typically sample from patients who meet "decipher_bin >= min(c2_values)"
  #    so that all subsets are valid. 
  subset_bootstrap <- data[data$decipher >= 0.3 & data$treatment %in% c("I","C"), ]
  
  # Split by treatment
  group_list <- split(subset_bootstrap, subset_bootstrap$treatment)
  
  # 3) Bootstrap loop
  for (m in seq_len(n_bootstrap)) {
    # Resample within each group
    resampled_data <- do.call(rbind, lapply(group_list, function(g) {
      g[sample(seq_len(nrow(g)), nrow(g), replace = TRUE), ]
    }))
    
    # For each c2, run analyze_c2_zstat_ordinal
    resampled_results <- lapply(c2_values, function(c2) {
      tmp <- tryCatch({
        analyze_c2_zstat_ordinal(resampled_data, c2)
      }, error = function(e) {
        return(list(log_or = NA, se_log_or = NA))
      })
      
      # Return minimal list structure
      if (is.null(tmp) || is.na(tmp$log_or) || is.na(tmp$se_log_or)) {
        return(list(log_or = NA, se_log_or = NA))
      } else {
        return(tmp)
      }
    })
    
    # Gather the log_or and se for each c2
    resampled_log_ors    <- sapply(resampled_results, `[[`, "log_or")
    resampled_se_log_ors <- sapply(resampled_results, `[[`, "se_log_or")
    

    studentized_z_stats <- (resampled_log_ors - log_ors) / resampled_se_log_ors
    
    # Replace NA/Inf with Inf so they won't pass min
    studentized_z_stats[is.na(studentized_z_stats)] <- Inf
    studentized_z_stats[is.infinite(studentized_z_stats)] <- Inf
    
    # Reorder to match sorted observed stats
    bootstrap_studentized_z[m, ] <- studentized_z_stats[original_order]
  }
  
  
  #### ---- 4)  GLOBAL & PER-THRESHOLD CRITICAL VALUES  -------------------- ####
  
  ## -- 4A.  Global (“min across thresholds”) distribution – used ONCE ---- ##
  global_min_vec <- apply(bootstrap_studentized_z, 1L, min, na.rm = TRUE)
  crit_global    <- stats::quantile(global_min_vec, probs = alpha, na.rm = TRUE)
  
  
  #### ---- 5)  STEP-DOWN RULE ------------------------------------------- ####
  
  # Observed Z-stats (already computed earlier)
  obs_z <- observed_z_stats        # numeric vector length = n_cutoffs
  
  # Identify most “significant” c1  (smallest Z => strongest NI evidence)
  idx_star <- which.min(obs_z)     # idx_star is position in c2_values
  c2_star  <- c2_values[idx_star]
  
  rejected <- accepted <- NULL      # to collect outcomes
  
  
  ## ----- 5A.  First test: compare to global critical value --------------
  if (obs_z[idx_star] < crit_global) {
    rejected <- c(rejected, c2_star)          # reject Sup null at c2*
    
    ## ---- 5B.  Test *smaller* thresholds sequentially at nominal α  ---- ##
    # thresholds ordered by their numeric value
    smaller_idx <- which(c2_values < c2_star) 
    smaller_idx <- smaller_idx[order(c2_values[smaller_idx], decreasing = TRUE)]   # descending
    # smaller_idx: Positions in c2_values
    
    for (j in seq_along(smaller_idx)) {
      ## ----  Gatekeeper: OR of the next bin(s) must be < 1 ----------- ##
      #  “smaller” = index of the first biomarker level *not yet included*
      
      next_idx   <- smaller_idx[j] # position in c2_values
      or_next    <- smoothed_OR[next_idx]
      this_c2    <- c2_values[next_idx] # actual biomarker bin, like 6 or 5.
      
      
      
      if (is.na(or_next) || or_next >= 0.8) {
        accepted <- this_c2   # stop: gatekeeper failed
        break
      }
      
      
      p_val <- lrt_Sup_pval(data, this_c2)
      if (!is.na(p_val) && p_val < alpha) {
        rejected <- c(rejected, this_c2)   # keep stepping down
      } else {
        accepted <- this_c2               # first acceptance stops the chain
        break
      }
    }
    
  } else {
    accepted <- c2_star         # failed at first hurdle → stop
  }
  
  
  
  return(list(
    accepted_c2_value  = accepted,
    rejected_c2_values = rejected,
    observed_z_stats   = obs_z,
    critical_global    = crit_global,
    OR                 = smoothed_OR
  ))
}




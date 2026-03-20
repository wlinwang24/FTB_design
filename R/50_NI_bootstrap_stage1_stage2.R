# ============================================================
# File: 50_NI_bootstrap_stage1_stage2.R
# Purpose:
#   Two-stage non-inferiority threshold selection procedure
#   for the FTB design.
#
# Contents:
#   1. c1_bootstrap_LRT_Smooth()
#      Implements:
#      - Stage 1: bootstrap-based threshold selection using
#        the minimum re-centered Wald-type Z-statistic across
#        candidate c1 values
#      - Stage 2: fixed-order sequential expansion using
#        smoothed odds ratios for gatekeeping and one-sided
#        likelihood ratio tests at the nominal level
#
# Notes:
#   The procedure first identifies the most favorable c1 and
#   performs a global test to control family-wise error across
#   the initial threshold search. If significant, the selected
#   threshold is expanded sequentially to larger c1 values.
# ============================================================



c1_bootstrap_LRT_Smooth <- function(data,
                                    c1_values,
                                    non_inferiority_margin,
                                    alpha = 0.05,
                                    n_bootstrap = 1000) {
  # 1) Get observed Z-stats for each c1
  observed_results <- lapply(c1_values, function(c1) {
    res <- analyze_c1_zstat_ordinal(
      data                = data,
      c1                  = c1,
      non_inferiority_margin = non_inferiority_margin
    )
    if (is.null(res) || is.na(res$z_stat)) {
      # If the result is missing or z_stat is NA, store Inf for the z_stat
      list(c1     = c1,
           z_stat = Inf,
           p_value= NA,
           or     = NA,
           log_or = NA,
           se_log_or = NA)
    } else {
      res
    }
  })
  
  # If everything is NA, bail out
  if (all(sapply(observed_results, function(x) is.na(x$z_stat)))) {
    stop("No valid results for any c1 value.")
  }
  
  # Extract the observed log(OR) and standard errors
  log_ors    <- sapply(observed_results, function(x) x$log_or)
  se_log_ors <- sapply(observed_results, function(x) x$se_log_or)
  
  # Compute the observed Z-stats for each c1
  observed_z_stats <- (log_ors - log(non_inferiority_margin)) / se_log_ors
  
  
  # Order them from smallest to largest
  original_order <- order(observed_z_stats, decreasing = FALSE)
  sorted_observed_z_stats <- observed_z_stats[original_order]
  
  # 2) Prepare for bootstrap
  n_cutoffs <- length(c1_values)
  smoothed_OR <- Smooth_OR_NI_tune(data = data, span = 1, plot = FALSE)$smoothed_OR
  
  
  # Studentized Z array: one row per bootstrap, one col per c1
  bootstrap_studentized_z <- matrix(NA, nrow = n_bootstrap, ncol = n_cutoffs)
  
  # Subset the data for D or C up to the max(c1_values) if relevant.
  # (Up to you whether to restrict to the largest bin.)
  # Or you may bootstrap from the entire data set (just for demonstration):
  subset_bootstrap <- data[data$decipher <= 0.7 & data$treatment %in% c("D","C"), ]
  
  # Split by treatment for separate sampling (if you want to preserve # of events in each arm)
  group_list <- split(subset_bootstrap, subset_bootstrap$treatment)
  
  # 3) Bootstrap loop
  for (m in seq_len(n_bootstrap)) {
    # Resample within each treatment group
    resampled_data <- do.call(rbind, lapply(group_list, function(g) {
      g[sample(seq_len(nrow(g)), nrow(g), replace = TRUE), ]
    }))
    
    # For each c1, run analyze_c1_zstat_ordinal()
    resampled_results <- lapply(c1_values, function(c1) {
      val <- tryCatch({
        analyze_c1_zstat_ordinal(
          data    = resampled_data,
          c1      = c1,
          non_inferiority_margin = non_inferiority_margin
        )
      }, error = function(e) {
        list(log_or = NA, se_log_or = NA)
      })
      if (is.null(val) || is.na(val$log_or) || is.na(val$se_log_or)) {
        list(log_or = NA, se_log_or = NA)
      } else {
        val
      }
    })
    
    # Gather the log_or and se for each c1
    resampled_log_ors    <- sapply(resampled_results, `[[`, "log_or")
    resampled_se_log_ors <- sapply(resampled_results, `[[`, "se_log_or")
    
    # Studentized difference:
    #    (resampled_log_or - observed_log_or) / resampled_se_log_or
    studentized_z_stats <- (resampled_log_ors - log_ors) / resampled_se_log_ors
    
    # Replace any NA/Inf with Inf so they don't break the min
    studentized_z_stats[is.na(studentized_z_stats)] <- Inf
    
    # Reorder to match the sorted (small->large) order
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
  idx_star <- which.min(obs_z)     
  c1_star  <- c1_values[idx_star]
  
  rejected <- accepted <- NULL      # to collect outcomes
  
  ## ----- 5A.  First test: compare to global critical value --------------
  if (obs_z[idx_star] < crit_global) {
    rejected <- c(rejected, c1_star)          # reject NI null at c1*
    
    ## ---- 5B.  Test *larger* thresholds sequentially at nominal α  ---- ##
    # thresholds ordered by their numeric value
    bigger_idx <- which(c1_values > c1_star)
    bigger_idx <- bigger_idx[order(c1_values[bigger_idx])]   # ascending
    
    for (j in seq_along(bigger_idx)) {
      
      ## ----  Gatekeeper: OR of the next bin(s) must be < 1.5 ----------- ##
      #  “higher” = index of the first Decipher bin *not yet included*
      next_idx <- bigger_idx[j]
      
      this_c1    <- c1_values[next_idx]
      or_next  <- smoothed_OR[this_c1]
      
      if (is.na(or_next) || or_next >= 1.3) {
        accepted <- this_c1   # stop: gatekeeper failed
        break
      }
      
      
      
      p_val <- lrt_NI_pval(data, this_c1, non_inferiority_margin)
      if (!is.na(p_val) && p_val < alpha) {
        rejected <- c(rejected, this_c1)   # keep stepping up
      } else {
        accepted <- this_c1                # first acceptance stops the chain
        break
      }
    }
    
  } else {
    accepted <- c1_star         # failed at first hurdle → stop
  }
  
  return(list(
    accepted_c1_value  = accepted,
    rejected_c1_values = rejected,
    observed_z_stats   = obs_z,
    critical_global    = crit_global,
    OR                 = smoothed_OR
  ))
}

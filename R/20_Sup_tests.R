# ============================================================
# File: 20_Sup_tests.R
# Purpose:
#   Superiority testing functions for the FTB design.
#
# Contents:
#   1. analyze_c2_zstat_ordinal()
#      Computes the one-sided Wald-type Z-statistic for a
#      candidate superiority threshold c2 using logistic
#      regression adjusted for ordinal Decipher bin.
#
#   2. lrt_Sup_pval()
#      Computes the one-degree-of-freedom likelihood ratio test
#      p-value for superiority at a fixed threshold c2.
#
# Notes:
#   These functions are used in:
#   - Stage 1 bootstrap-based threshold selection
#   - Stage 2 fixed-order sequential expansion
# ============================================================


analyze_c2_zstat_ordinal <- function(data, c2) {
  # Subset for decipher_bin >= c2 and treatments I or C
  subset_data <- data[data$decipher_bin >= c2 & data$treatment %in% c("I", "C"), ]
  
  # Make decipher_bin an ordered factor.  
  subset_data$decipher_bin <- factor(
    subset_data$decipher_bin,
    levels = sort(unique(subset_data$decipher_bin)),
    ordered = TRUE
  )
  
  # Force treatment factor with levels "C" (control) then "I" (investigational/test)
  subset_data$treatment <- factor(subset_data$treatment, levels = c("C", "I"))
  
  # Check we have at least 2 events in each group
  events_in_I <- sum(subset_data$treatment == "I" & subset_data$status_5yr == 1, na.rm = TRUE)
  events_in_C <- sum(subset_data$treatment == "C" & subset_data$status_5yr == 1, na.rm = TRUE)
  if (events_in_I < 2 || events_in_C < 2) {
    return(list(c2     = c2,
                z_stat = NA,
                p_value = NA,
                or = NA,
                log_or = NA,
                se_log_or = NA))
  }
  
  if (length(unique(subset_data$decipher_bin)) < 2) {
    ## Only one bin level present – drop the factor
    glm_formula <- status_5yr ~ treatment
  } else {
    glm_formula <- status_5yr ~ decipher_bin + treatment
  }
  
  
  glm_model <- tryCatch(
    glm(glm_formula, data = subset_data, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  
  # If the model fails to converge or has other issues, return NA
  if (is.null(glm_model)) {
    return(list(c2     = c2,
                z_stat = NA,
                p_value = NA,
                or = NA,
                log_or = NA,
                se_log_or = NA))
  }
  
  coefs <- summary(glm_model)$coefficients
  # We need the row for 'treatmentI'
  if (!"treatmentI" %in% rownames(coefs)) {
    return(list(c2     = c2,
                z_stat = NA,
                p_value = NA,
                or = NA,
                log_or = NA,
                se_log_or = NA))
  }
  
  log_or    <- coefs["treatmentI", "Estimate"]
  se_log_or <- coefs["treatmentI", "Std. Error"]
  if (is.na(log_or) || is.na(se_log_or) ||
      is.infinite(log_or) || is.infinite(se_log_or)) {
    return(list(c2     = c2,
                z_stat = NA,
                p_value = NA,
                or = NA,
                log_or = NA,
                se_log_or = NA))
  }
  
  # OR is exp(log_or)
  or <- exp(log_or)
  
  # *** Superiority Z-statistic: 
  #    H0: OR >= 1  <=> log(OR) >= 0
  #    => z = log_or / se_log_or
  z_stat <- log_or / se_log_or
  
  # One-sided p-value (small p => strong evidence OR < 1)
  one_sided_p_value <- pnorm(z_stat)
  
  return(list(c2       = c2,
              z_stat   = z_stat,
              p_value  = one_sided_p_value,
              or       = or,
              log_or   = log_or,
              se_log_or = se_log_or))
}


## ------------------------------------------------------------------ ##
##  one-df likelihood-ratio Sup test at a single cutoff   ##
## ------------------------------------------------------------------ ##
lrt_Sup_pval <- function(data, c2) {
  
  ## 1)  keep Decipher ≤ c1   and   treatments C / D
  sub <- data[data$decipher_bin >= c2 & data$treatment %in% c("C","I"), ]
  
  sub$decipher_bin <- factor(sub$decipher_bin,
                             levels = sort(unique(sub$decipher_bin)),
                             ordered = TRUE)
  sub$treat_num  <- as.numeric(sub$treatment == "I")     # 1 = I, 0 = C
  
  
  
  ## full vs. null (offset = log(margin) * treat_num)
  mod_full <- tryCatch(
    glm(status_5yr ~ decipher_bin + treat_num,
        data = sub, family = binomial()),
    error = function(e) NULL)
  
  mod_null <- tryCatch(
    glm(status_5yr ~ decipher_bin,
        data = sub, family = binomial()),
    error = function(e) NULL)
  
  if (is.null(mod_full) || is.null(mod_null))
    return(NA)
  
  lrt_stat <- 2 * (logLik(mod_full) - logLik(mod_null))
  pchisq(lrt_stat, df = 1, lower.tail = FALSE)   # one-df χ² p-value
  
  
}











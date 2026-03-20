# ============================================================
# File: 10_NI_tests.R
# Purpose:
#   Non-inferiority testing functions for the FTB design.
#
# Contents:
#   1. analyze_c1_zstat_ordinal()
#      Computes the one-sided Wald-type Z-statistic for a
#      candidate non-inferiority threshold c1 using logistic
#      regression adjusted for ordinal Decipher bin.
#
#   2. lrt_NI_pval()
#      Computes the one-degree-of-freedom likelihood ratio test
#      p-value for non-inferiority at a fixed threshold c1.
#
# Notes:
#   These functions are used in:
#   - Stage 1 bootstrap-based threshold selection
#   - Stage 2 fixed-order sequential expansion
# ============================================================

analyze_c1_zstat_ordinal <- function(data, 
                                     c1, 
                                     non_inferiority_margin) {
  
  
  
  # 2) Subset for decipher_bin <= c1 and only treatments D or C
  #    Here we assume c1 is an integer bin index (1..10).
  subset_data <- data[data$decipher_bin <= c1 & data$treatment %in% c("D", "C"), ]
  
  # 3) Convert decipher_bin to an ordered factor
  #    so that R sees it as an ordinal variable
  subset_data$decipher_bin <- factor(subset_data$decipher_bin, 
                                     levels = sort(unique(subset_data$decipher_bin)),
                                     ordered = TRUE)
  
  # 4) Ensure that treatment is a factor with levels C then D
  subset_data$treatment <- factor(subset_data$treatment, levels = c("C","D"))
  
  # 5) Check we have at least 2 events in each group for stable logistic regression
  events_in_D <- sum(subset_data$treatment == "D" & subset_data$status_5yr == 1, na.rm=TRUE)
  events_in_C <- sum(subset_data$treatment == "C" & subset_data$status_5yr == 1, na.rm=TRUE)
  
  if (events_in_D < 2 || events_in_C < 2) {
    # Not enough events for a stable model
    return(list(
      c1       = c1,
      z_stat   = NA,
      p_value  = NA,
      or       = NA,
      log_or   = NA,
      se_log_or= NA
    ))
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
  
  if (is.null(glm_model)) {
    return(list(
      c1       = c1,
      z_stat   = NA,
      p_value  = NA,
      or       = NA,
      log_or   = NA,
      se_log_or= NA
    ))
  }
  
  coefs <- summary(glm_model)$coefficients
  if (!"treatmentD" %in% rownames(coefs)) {
    return(list(
      c1       = c1,
      z_stat   = NA,
      p_value  = NA,
      or       = NA,
      log_or   = NA,
      se_log_or= NA
    ))
  }
  
  log_or    <- coefs["treatmentD", "Estimate"]
  se_log_or <- coefs["treatmentD", "Std. Error"]
  
  if (is.na(log_or) || is.na(se_log_or) || is.infinite(log_or) || is.infinite(se_log_or)) {
    return(list(
      c1       = c1,
      z_stat   = NA,
      p_value  = NA,
      or       = NA,
      log_or   = NA,
      se_log_or= NA
    ))
  }
  
  or <- exp(log_or)
  
  # 8) Non-inferiority Z-statistic:
  #    H0: OR >= margin   vs.   H1: OR < margin
  #    equivalently  H0: log(OR) >= log(margin).
  #    so  z = [log_or - log(margin)] / SE
  z_stat <- (log_or - log(non_inferiority_margin)) / se_log_or
  
  # One-sided p-value
  one_sided_p_value <- pnorm(z_stat)
  
  return(list(
    c1       = c1,
    z_stat   = z_stat,
    p_value  = one_sided_p_value,
    or       = or,
    log_or   = log_or,
    se_log_or= se_log_or
  ))
}



## ------------------------------------------------------------------ ##
##  one-df likelihood-ratio NI test at a single cutoff   ##
## ------------------------------------------------------------------ ##
lrt_NI_pval <- function(data, c1, non_inferiority_margin) {
  sub <- data[data$decipher_bin <= c1 & data$treatment %in% c("C","D"), ]
  
  
  
  sub$decipher_bin <- factor(sub$decipher_bin,
                             levels = sort(unique(sub$decipher_bin)),
                             ordered = TRUE)
  sub$treat_num  <- as.numeric(sub$treatment == "D")     # 1 = D, 0 = C
  log_m <- log(non_inferiority_margin)
  
  
  ## full vs. null (offset = log(margin) * treat_num)
  mod_full <- tryCatch(
    glm(status_5yr ~ decipher_bin + treat_num,
        data = sub, family = binomial()),
    error = function(e) NULL)
  
  mod_null <- tryCatch(
    glm(status_5yr ~ decipher_bin +
          offset(log_m * treat_num),
        data = sub, family = binomial()),
    error = function(e) NULL)
  
  if (is.null(mod_full) || is.null(mod_null))
    return(NA)
  
  lrt_stat <- 2 * (logLik(mod_full) - logLik(mod_null))
  pchisq(lrt_stat, df = 1, lower.tail = FALSE)   # one-df χ² p-value
  
}










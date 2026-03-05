# ============================================================
# File: 80_simulation_functions.R
# Purpose:
#   Simulation functions for evaluating the flexible threshold
#   biomarker-based (FTB) design under manuscript Scenarios C-F.
#
# Notes:
#   1. This file defines shared helper functions and
#      scenario-specific simulation functions.
#   2. It does NOT run parallel simulations directly.
#      Parallel execution should be handled in scenario-specific
#      runner scripts.
#   3. These functions assume that `prob_table` has already been
#      defined in the calling script, e.g.:
#         prob_table <- scenario_C$prob_table
#   4. These functions also assume that the testing / smoothing /
#      bootstrap functions have already been sourced.
# ============================================================


# ============================================================
# Common bin setup
# ============================================================
bin_breaks <- seq(0, 1, by = 0.1)


# ============================================================
# Dynamic treatment assignment
# ============================================================
assign_treatment_dynamic <- function(decipher, sup_done) {
  if (!sup_done) {
    if (decipher <= 0.3) {
      available_treatments <- c("D", "C")
    } else if (decipher <= 0.7) {
      available_treatments <- c("D", "I", "C")
    } else {
      available_treatments <- c("I", "C")
    }
  } else {
    if (decipher <= 0.7) {
      available_treatments <- c("D", "C")
    } else {
      return(NA)
    }
  }
  
  sample(available_treatments, size = 1)
}


# ============================================================
# Simulate binary endpoint from scenario-specific probability table
# ============================================================
simulate_ordinal <- function(decipher, treatment, prob_table) {
  
  idx <- findInterval(decipher, bin_breaks, rightmost.closed = TRUE)
  
  if (treatment == "C") {
    p_event <- prob_table$Prob_Control[idx]
  } else if (treatment == "I") {
    p_event <- prob_table$Prob_Treat_I[idx]
  } else if (treatment == "D") {
    p_event <- prob_table$Prob_Treat_D[idx]
  } else {
    stop("Unknown treatment!")
  }
  
  rbinom(1, 1, p_event)
}


# ============================================================
# Internal helper:
# Run one simulated trial given target sample sizes and max total
# ============================================================
.run_ftb_trial <- function(sim_id,
                           prob_table,
                           target_Sup,
                           target_NI,
                           max_patients,
                           non_inferiority_margin = 1.5,
                           n_bootstrap = 1000,
                           alpha_sup = 0.025,
                           alpha_ni = 0.05) {
  
  # Storage vectors
  decipher_scores   <- numeric(max_patients)
  treatment_assigns <- character(max_patients)
  status_5yr_all    <- numeric(max_patients)
  used_in_sup       <- logical(max_patients)
  used_in_ni        <- logical(max_patients)
  
  # Counters
  N_total    <- 0
  N_SupCount <- 0
  N_NI       <- 0
  sup_done   <- FALSE
  
  while (N_total < max_patients) {
    
    # 1) Generate Decipher
    decipher <- runif(1, 0, 1)
    
    # 2) Assign treatment or skip
    trt <- assign_treatment_dynamic(decipher, sup_done)
    if (is.na(trt)) {
      next
    }
    
    # 3) Valid patient
    N_total <- N_total + 1
    
    # 4) Simulate 5-year event
    event_5yr <- simulate_ordinal(
      decipher   = decipher,
      treatment  = trt,
      prob_table = prob_table
    )
    
    # 5) Store
    decipher_scores[N_total]   <- decipher
    treatment_assigns[N_total] <- trt
    status_5yr_all[N_total]    <- event_5yr
    
    # 6A) NI inclusion: Decipher <= 0.7 and D/C
    ni_incl <- (decipher <= 0.7 && trt %in% c("D", "C"))
    used_in_ni[N_total] <- ni_incl
    if (ni_incl) {
      N_NI <- N_NI + 1
    }
    
    # 6B) Sup inclusion: Decipher >= 0.3 and I/C while Sup still open
    sup_incl <- (!sup_done && decipher >= 0.3 && trt %in% c("I", "C"))
    used_in_sup[N_total] <- sup_incl
    if (sup_incl) {
      N_SupCount <- N_SupCount + 1
      if (N_SupCount >= target_Sup) {
        sup_done <- TRUE
      }
    }
  }
  
  # Final simulated dataset
  simulated_data <- data.frame(
    patient_ID   = 1:max_patients,
    decipher     = decipher_scores,
    treatment    = treatment_assigns,
    status_5yr   = status_5yr_all,
    used_in_sup  = used_in_sup,
    used_in_ni   = used_in_ni
  )
  
  simulated_data$decipher_bin <- cut(
    simulated_data$decipher,
    breaks = bin_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE
  )
  
  # Counts
  n_sup <- sum(simulated_data$used_in_sup)
  n_ni  <- sum(simulated_data$used_in_ni)
  
  control_used_in_both <- sum(
    simulated_data$treatment == "C" &
      simulated_data$decipher <= 0.7 &
      simulated_data$decipher >= 0.3 &
      simulated_data$used_in_sup
  )
  
  control_shared <- sum(
    simulated_data$treatment == "C" &
      simulated_data$decipher >= 0.3 &
      simulated_data$decipher <= 0.7
  )
  
  # Initialize outputs
  c2_hat         <- 1
  c2_significant <- FALSE
  c2_decision    <- "Null"
  c2_max         <- NA
  
  c1_hat         <- 0
  c1_significant <- FALSE
  c1_decision    <- "Null"
  c1_min         <- NA
  
  # Candidate cutoffs
  c2_values_ord <- 4:7
  c1_values_ord <- 4:7
  
  # ------------------------------------------------------------
  # Superiority analysis
  # ------------------------------------------------------------
  data_Sup <- subset(simulated_data, used_in_sup)
  
  c2_bootstrap_results <- c2_bootstrap_LRT_Smooth(
    data         = data_Sup,
    c2_values    = c2_values_ord,
    alpha        = alpha_sup,
    n_bootstrap  = n_bootstrap
  )
  
  Sup_or_values <- c2_bootstrap_results$OR
  Sup_OR <- rep(NA, 7)
  Sup_OR <- Sup_or_values[1:7]
  
  if (!is.null(c2_bootstrap_results$rejected_c2_values) &&
      length(c2_bootstrap_results$rejected_c2_values) > 0) {
    c2_hat <- min(unlist(c2_bootstrap_results$rejected_c2_values), na.rm = TRUE)
    c2_max <- max(unlist(c2_bootstrap_results$rejected_c2_values), na.rm = TRUE)
    c2_significant <- TRUE
    c2_decision <- "Alternative"
  } else {
    c2_hat <- 1
    c2_significant <- FALSE
    c2_decision <- "Null"
  }
  
  # ------------------------------------------------------------
  # Non-inferiority analysis
  # ------------------------------------------------------------
  data_NI <- subset(simulated_data, used_in_ni)
  
  c1_bootstrap_results <- c1_bootstrap_LRT_Smooth(
    data                    = data_NI,
    c1_values               = c1_values_ord,
    non_inferiority_margin  = non_inferiority_margin,
    alpha                   = alpha_ni,
    n_bootstrap             = n_bootstrap
  )
  
  NI_or_values <- c1_bootstrap_results$OR
  NI_OR <- rep(NA, 7)
  NI_OR <- NI_or_values[1:7]
  
  if (!is.null(c1_bootstrap_results$rejected_c1_values) &&
      length(c1_bootstrap_results$rejected_c1_values) > 0) {
    c1_hat <- max(unlist(c1_bootstrap_results$rejected_c1_values), na.rm = TRUE)
    c1_min <- min(unlist(c1_bootstrap_results$rejected_c1_values), na.rm = TRUE)
    c1_significant <- TRUE
    c1_decision <- "Alternative"
  } else {
    c1_hat <- 0
    c1_significant <- FALSE
    c1_decision <- "Null"
  }
  
  # ------------------------------------------------------------
  # Return trial-level result
  # ------------------------------------------------------------
  return(list(
    sim_id = sim_id,
    n_sup = n_sup,
    n_ni = n_ni,
    control_used_in_both = control_used_in_both,
    control_shared = control_shared,
    
    c2_hat = c2_hat,
    c2_max = c2_max,
    c2_significant = c2_significant,
    c2_decision = c2_decision,
    
    Sup_0.4 = Sup_OR[1],
    Sup_0.5 = Sup_OR[2],
    Sup_0.6 = Sup_OR[3],
    Sup_0.7 = Sup_OR[4],
    Sup_0.8 = Sup_OR[5],
    Sup_0.9 = Sup_OR[6],
    Sup_1.0 = Sup_OR[7],
    
    c1_hat = c1_hat,
    c1_min = c1_min,
    c1_significant = c1_significant,
    c1_decision = c1_decision,
    
    NI_0.1 = NI_OR[1],
    NI_0.2 = NI_OR[2],
    NI_0.3 = NI_OR[3],
    NI_0.4 = NI_OR[4],
    NI_0.5 = NI_OR[5],
    NI_0.6 = NI_OR[6],
    NI_0.7 = NI_OR[7]
  ))
}


# ============================================================
# Scenario C:
# Global null, used for Type I error assessment
# Same sample-size benchmarked on conventional fixed-threshold design
# ============================================================
simulate_Bootstrap_LRT_Smooth <- function(sim_id) {
  
  n_bootstrap <- 1000
  non_inferiority_margin <- 1.5
  
  c_val <- 1.095
  target_Sup <- c_val * 1360
  target_NI  <- c_val * 2340
  max_patients <- 3700
  
  .run_ftb_trial(
    sim_id                  = sim_id,
    prob_table              = prob_table,
    target_Sup              = target_Sup,
    target_NI               = target_NI,
    max_patients            = max_patients,
    non_inferiority_margin  = non_inferiority_margin,
    n_bootstrap             = n_bootstrap,
    alpha_sup               = 0.025,
    alpha_ni                = 0.05
  )
}


# ============================================================
# Scenario D:
# Overall power using same sample size benchmarked on the
# conventional fixed-threshold design
# ============================================================
ha_benchmark_simulate_LRT_Smooth <- function(sim_id) {
  
  n_bootstrap <- 1000
  non_inferiority_margin <- 1.5
  
  c_val <- 1.095
  target_Sup <- c_val * 1360
  target_NI  <- c_val * 2340
  max_patients <- 3700
  
  .run_ftb_trial(
    sim_id                  = sim_id,
    prob_table              = prob_table,
    target_Sup              = target_Sup,
    target_NI               = target_NI,
    max_patients            = max_patients,
    non_inferiority_margin  = non_inferiority_margin,
    n_bootstrap             = n_bootstrap,
    alpha_sup               = 0.025,
    alpha_ni                = 0.05
  )
}


# ============================================================
# Scenario D:
# Find sample size achieving 85% targeted power for both
# superiority and non-inferiority
# ============================================================
ha_simulate_LRT_Smooth <- function(sim_id) {
  
  n_bootstrap <- 1000
  non_inferiority_margin <- 1.5
  
  target_Sup <- 1280
  target_NI  <- 3000
  max_patients <- 2860
  
  .run_ftb_trial(
    sim_id                  = sim_id,
    prob_table              = prob_table,
    target_Sup              = target_Sup,
    target_NI               = target_NI,
    max_patients            = max_patients,
    non_inferiority_margin  = non_inferiority_margin,
    n_bootstrap             = n_bootstrap,
    alpha_sup               = 0.025,
    alpha_ni                = 0.05
  )
}


# ============================================================
# Scenario E:
# Find sample size achieving 85% targeted power for both
# superiority and non-inferiority
# ============================================================
case1_simulate_LRT_Smooth <- function(sim_id) {
  
  n_bootstrap <- 1000
  non_inferiority_margin <- 1.5
  
  target_Sup <- 1650
  target_NI  <- 4000
  max_patients <- 4000
  
  .run_ftb_trial(
    sim_id                  = sim_id,
    prob_table              = prob_table,
    target_Sup              = target_Sup,
    target_NI               = target_NI,
    max_patients            = max_patients,
    non_inferiority_margin  = non_inferiority_margin,
    n_bootstrap             = n_bootstrap,
    alpha_sup               = 0.025,
    alpha_ni                = 0.05
  )
}


# ============================================================
# Scenario F:
# Find sample size achieving 85% targeted power for both
# superiority and non-inferiority
# ============================================================
case2_simulate_LRT_Smooth <- function(sim_id) {
  
  n_bootstrap <- 1000
  non_inferiority_margin <- 1.5
  
  target_Sup <- 2700
  target_NI  <- 5000
  max_patients <- 5500
  
  .run_ftb_trial(
    sim_id                  = sim_id,
    prob_table              = prob_table,
    target_Sup              = target_Sup,
    target_NI               = target_NI,
    max_patients            = max_patients,
    non_inferiority_margin  = non_inferiority_margin,
    n_bootstrap             = n_bootstrap,
    alpha_sup               = 0.025,
    alpha_ni                = 0.05
  )
}
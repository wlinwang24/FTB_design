# ============================================================
# Script: run_scenarioC_type1.R
# Purpose:
#   Run Scenario C (global null) to assess type I error of the
#   flexible threshold biomarker-based design.
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# Source project files
# ------------------------------------------------------------
source("R/00_packages.R")
source("R/10_NI_tests.R")
source("R/20_Sup_tests.R")
source("R/30_NI_smoothing.R")
source("R/40_Sup_smoothing.R")
source("R/50_NI_bootstrap_stage1_stage2.R")
source("R/60_Sup_bootstrap_stage1_stage2.R")
source("R/70_scenario_definitions.R")
source("R/80_simulation_functions.R")

# ------------------------------------------------------------
# Scenario-specific truth
# ------------------------------------------------------------
prob_table <- scenario_C$prob_table

# ------------------------------------------------------------
# Parallel setup
# ------------------------------------------------------------
library(doParallel)
library(doRNG)

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ------------------------------------------------------------
# Run simulations
# ------------------------------------------------------------
n_simulations <- 2000
set.seed(10012024)

results_list <- foreach(
  sim_id = 1:n_simulations,
  .packages = c("survival")
) %dorng% {
  simulate_Bootstrap_LRT_Smooth(sim_id)
}

stopCluster(cl)
rm(cl)
gc()

# ------------------------------------------------------------
# Convert and save
# ------------------------------------------------------------
h0_sim2000 <- do.call(rbind, lapply(results_list, as.data.frame))


write.csv(h0_sim2000, "results/h0_sim2000.csv", row.names = FALSE)

# ------------------------------------------------------------
# Quick checks
# ------------------------------------------------------------
sum(h0_sim2000$c2_significant == TRUE, na.rm = TRUE)
sum(h0_sim2000$c1_significant == TRUE, na.rm = TRUE)

table(h0_sim2000$c2_hat, useNA = "ifany")
table(h0_sim2000$c1_hat, useNA = "ifany")
# ============================================================
# File: 90_scenario_results.R
# Purpose:
#   Summarize simulation results for the FTB design and prepare
#   data objects used in manuscript tables and figures.
#
# Contents:
#   1. Type I error / liberal power summaries
#   2. Threshold-selection frequency data preparation
#   3. Cumulative inclusion data preparation
#
# Notes:
#   Plotting itself is handled in 95_figures.R.
# ============================================================


# ============================================================
# Common mapping objects
# ============================================================
all_bins <- seq(0.05, 0.95, by = 0.10)

sup_map <- c(
  `1` = NA,
  `4` = 0.35,
  `5` = 0.45,
  `6` = 0.55,
  `7` = 0.65,
  `8` = 0.75
)

ni_map <- c(
  `0` = NA,
  `4` = 0.35,
  `5` = 0.45,
  `6` = 0.55,
  `7` = 0.65
)


# ============================================================
# Part 1: Simple operating characteristics
# ============================================================

# ------------------------------------------------------------
# Scenario C: Type I error (global null)
# ------------------------------------------------------------
type1_error_scenarioC <- list(
  sup_count = sum(h0_sim2000$c2_significant == TRUE, na.rm = TRUE),
  ni_count  = sum(h0_sim2000$c1_significant == TRUE, na.rm = TRUE),
  sup_rate  = mean(h0_sim2000$c2_significant == TRUE, na.rm = TRUE),
  ni_rate   = mean(h0_sim2000$c1_significant == TRUE, na.rm = TRUE)
)


# ------------------------------------------------------------
# Liberal power helper
# ------------------------------------------------------------
compute_liberal_power <- function(file) {
  list(
    sup_count = sum(file$c2_significant == TRUE, na.rm = TRUE),
    ni_count  = sum(file$c1_significant == TRUE, na.rm = TRUE),
    sup_rate  = mean(file$c2_significant == TRUE, na.rm = TRUE),
    ni_rate   = mean(file$c1_significant == TRUE, na.rm = TRUE)
  )
}


# ------------------------------------------------------------
# Scenario D benchmarked at standard fixed-threshold sample size
# ------------------------------------------------------------
liberal_power_scenarioD_benchmark <- compute_liberal_power(ha_benched_sim2000)

# ------------------------------------------------------------
# Scenario D targeted-power sample size
# ------------------------------------------------------------
liberal_power_scenarioD_targeted <- compute_liberal_power(ha_sup1280tot2860_sim2000)

# ------------------------------------------------------------
# Scenario E targeted-power sample size
# ------------------------------------------------------------
liberal_power_scenarioE <- compute_liberal_power(case1_sup1650tot4000_sim2000)

# ------------------------------------------------------------
# Scenario F targeted-power sample size
# ------------------------------------------------------------
liberal_power_scenarioF <- compute_liberal_power(case2_sup2700tot5500_sim2000)


# ============================================================
# Part 2: Helper functions for threshold-selection frequencies
# ============================================================

# ------------------------------------------------------------
# Prepare threshold-selection data from one simulation dataset
# ------------------------------------------------------------
prepare_threshold_selection_data <- function(file) {
  file %>%
    transmute(
      sim_id = row_number(),
      th_sup = unname(sup_map[as.character(c2_hat)]),
      th_ni  = unname(ni_map[as.character(c1_hat)])
    )
}


# ------------------------------------------------------------
# Convert threshold selections to plotting data
# ------------------------------------------------------------
prepare_plot_dat_selection <- function(df_threshold, n_sim = 2000) {
  df_threshold %>%
    pivot_longer(c(th_sup, th_ni), names_to = "trial", values_to = "threshold") %>%
    mutate(trial = recode(trial, th_sup = "Superiority", th_ni = "Non-inferiority")) %>%
    filter(!is.na(threshold)) %>%
    count(trial, threshold, name = "n") %>%
    complete(trial, threshold = all_bins, fill = list(n = 0)) %>%
    group_by(trial) %>%
    mutate(pct = 100 * n / n_sim) %>%
    ungroup()
}


# ------------------------------------------------------------
# Prepare cumulative inclusion data from selection data
# ------------------------------------------------------------
prepare_plot_dat_cumulative <- function(plot_dat,
                                        sup_domain = c(0.35, 0.95),
                                        ni_domain  = c(0.05, 0.65)) {
  plot_dat %>%
    mutate(thr = as.numeric(as.character(threshold))) %>%
    arrange(trial, thr) %>%
    group_by(trial) %>%
    mutate(
      pct_in_domain = case_when(
        trial == "Superiority" ~ if_else(thr >= sup_domain[1] & thr <= sup_domain[2], pct, 0),
        trial == "Non-inferiority" ~ if_else(thr >= ni_domain[1] & thr <= ni_domain[2], pct, 0)
      ),
      cum_calc = if (unique(trial) == "Superiority") {
        cumsum(pct_in_domain)
      } else {
        rev(cumsum(rev(pct_in_domain)))
      }
    ) %>%
    ungroup() %>%
    group_by(trial, thr) %>%
    summarise(cum_pct = max(cum_calc), .groups = "drop") %>%
    mutate(threshold = factor(thr, levels = all_bins))
}


# ============================================================
# Part 3: Scenario D benchmarked at standard fixed-threshold size
# ============================================================

df_ha_benchmark <- prepare_threshold_selection_data(ha_benched_sim2000)

plot_dat_ha_benchmark <- prepare_plot_dat_selection(
  df_threshold = df_ha_benchmark,
  n_sim = nrow(ha_benched_sim2000)
)

ha_benchmark_true_sup_bins <- c(0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
ha_benchmark_true_ni_bins  <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)

plot_dat_cum_ha_benchmark <- prepare_plot_dat_cumulative(
  plot_dat = plot_dat_ha_benchmark,
  sup_domain = c(0.35, 0.95),
  ni_domain  = c(0.05, 0.65)
)


# ============================================================
# Part 4: Scenario D targeted-power sample size
# ============================================================

df_ha_targeted <- prepare_threshold_selection_data(ha_sup1280tot2860_sim2000)

plot_dat_ha_targeted <- prepare_plot_dat_selection(
  df_threshold = df_ha_targeted,
  n_sim = nrow(ha_sup1280tot2860_sim2000)
)

ha_targeted_true_sup_bins <- c(0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
ha_targeted_true_ni_bins  <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)

plot_dat_cum_ha_targeted <- prepare_plot_dat_cumulative(
  plot_dat = plot_dat_ha_targeted,
  sup_domain = c(0.35, 0.95),
  ni_domain  = c(0.05, 0.65)
)


# ============================================================
# Part 5: Scenario E targeted-power sample size
# ============================================================

df_case1 <- prepare_threshold_selection_data(case1_sup1650tot4000_sim2000)

plot_dat_case1 <- prepare_plot_dat_selection(
  df_threshold = df_case1,
  n_sim = nrow(case1_sup1650tot4000_sim2000)
)

case1_true_sup_bins <- c(0.65, 0.75, 0.85, 0.95)
case1_true_ni_bins  <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65)

plot_dat_cum_case1 <- prepare_plot_dat_cumulative(
  plot_dat = plot_dat_case1,
  sup_domain = c(0.35, 0.95),
  ni_domain  = c(0.05, 0.65)
)


# ============================================================
# Part 6: Scenario F targeted-power sample size
# ============================================================

df_case2 <- prepare_threshold_selection_data(case2_sup2700tot5500_sim2000)

plot_dat_case2 <- prepare_plot_dat_selection(
  df_threshold = df_case2,
  n_sim = nrow(case2_sup2700tot5500_sim2000)
)

case2_true_sup_bins <- c(0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
case2_true_ni_bins  <- c(0.05, 0.15, 0.25, 0.35, 0.45)

plot_dat_cum_case2 <- prepare_plot_dat_cumulative(
  plot_dat = plot_dat_case2,
  sup_domain = c(0.35, 0.95),
  ni_domain  = c(0.05, 0.65)
)


# ============================================================
# Part 7: Optional summary table across scenarios
# ============================================================

operating_characteristics_summary <- data.frame(
  Scenario = c(
    "C: Global null",
    "D: Benchmark sample size",
    "D: Targeted-power sample size",
    "E: Targeted-power sample size",
    "F: Targeted-power sample size"
  ),
  Sup_Count = c(
    type1_error_scenarioC$sup_count,
    liberal_power_scenarioD_benchmark$sup_count,
    liberal_power_scenarioD_targeted$sup_count,
    liberal_power_scenarioE$sup_count,
    liberal_power_scenarioF$sup_count
  ),
  Sup_Rate = c(
    type1_error_scenarioC$sup_rate,
    liberal_power_scenarioD_benchmark$sup_rate,
    liberal_power_scenarioD_targeted$sup_rate,
    liberal_power_scenarioE$sup_rate,
    liberal_power_scenarioF$sup_rate
  ),
  NI_Count = c(
    type1_error_scenarioC$ni_count,
    liberal_power_scenarioD_benchmark$ni_count,
    liberal_power_scenarioD_targeted$ni_count,
    liberal_power_scenarioE$ni_count,
    liberal_power_scenarioF$ni_count
  ),
  NI_Rate = c(
    type1_error_scenarioC$ni_rate,
    liberal_power_scenarioD_benchmark$ni_rate,
    liberal_power_scenarioD_targeted$ni_rate,
    liberal_power_scenarioE$ni_rate,
    liberal_power_scenarioF$ni_rate
  )
)
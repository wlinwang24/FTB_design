# ============================================================
# File: 70_scenario_definitions.R
# Purpose:
#   Define manuscript simulation scenarios for evaluation of
#   the flexible threshold biomarker-based (FTB) design.
#
# Notes:
#   Scenario C is the global null scenario used to assess
#   type I error. Scenarios D-F are alternative scenarios used
#   to assess operating characteristics under broad and narrow
#   regions of non-inferiority and superiority.
# ============================================================

# ============================================================
# Common Biomarker level setup
# ============================================================
decipher_breaks <- seq(0, 1, by = 0.1)

decipher_bins <- cut(
  x = decipher_breaks,
  breaks = decipher_breaks,
  include.lowest = TRUE,
  right = TRUE
)

bin_midpoints <- seq(0.05, 0.95, by = 0.1)



# ============================================================
# Scenario III: global null
# Purpose:
#   Used to assess type I error of the flexible design.
#   No non-inferiority and no superiority anywhere.
# ============================================================
prob_treat_D_III <- c(0.1014, 0.1223, 0.1469, 0.1751, 0.2080,
                    0.2464, 0.2900, 0.3392, 0.3958, 0.4598)

prob_control_III <- c(0.070, 0.085, 0.103, 0.124, 0.149,
                    0.179, 0.214, 0.255, 0.304, 0.362)

prob_treat_I_III <- c(0.070, 0.085, 0.103, 0.124, 0.149,
                    0.179, 0.214, 0.255, 0.304, 0.362)

scenario_III <- list(
  scenario_id = "III",
  manuscript_label = "Scenario III",
  scenario_type = "global_null",
  description = "Global null: no non-inferiority and no superiority anywhere.",
  prob_treat_D = prob_treat_D_III,
  prob_control = prob_control_III,
  prob_treat_I = prob_treat_I_III,
  prob_table = data.frame(
    Decipher_Midpoint = bin_midpoints,
    Prob_Treat_D = prob_treat_D_III,
    Prob_Control = prob_control_III,
    Prob_Treat_I = prob_treat_I_III
  )
)


# ============================================================
# Scenario IV: alternative scenario
# Purpose:
#   Used to assess performance of the flexible design under an
#   alternative with a relatively broad region structure.
# ============================================================
prob_treat_D_IV <- c(0.070, 0.085, 0.103, 0.124, 0.149,
                    0.179, 0.214, 0.255, 0.304, 0.362)

prob_control_IV <- c(0.070, 0.085, 0.103, 0.124, 0.149,
                    0.179, 0.214, 0.255, 0.304, 0.362)

prob_treat_I_IV <- c(0.0466, 0.0569, 0.0695, 0.0843, 0.1022,
                    0.1241, 0.1504, 0.1820, 0.2211, 0.2694)

scenario_IV <- list(
  scenario_id = "IV",
  manuscript_label = "Scenario IV",
  scenario_type = "alternative",
  description = "Alternative scenario with broad region structure where a fixed cutoff may perform reasonably well.",
  prob_treat_D = prob_treat_D_IV,
  prob_control = prob_control_IV,
  prob_treat_I = prob_treat_I_IV,
  prob_table = data.frame(
    Decipher_Midpoint = bin_midpoints,
    Prob_Treat_D = prob_treat_D_IV,
    Prob_Control = prob_control_IV,
    Prob_Treat_I = prob_treat_I_IV
  )
)



# ============================================================
# Scenario V: alternative scenario
# Purpose:
#   Used to assess performance of the flexible design under an
#   alternative with a more localized or asymmetric region
#   pattern.
# ============================================================
prob_treat_D_V <- c(0.05, 0.067, 0.083, 0.103, 0.13,
                    0.15, 0.226, 0.37, 0.39, 0.43)

prob_control_V <- c(0.05, 0.067, 0.083, 0.103, 0.13,
                    0.15, 0.226, 0.281, 0.333, 0.4)

prob_treat_I_V <- c(0.05, 0.067, 0.083, 0.103, 0.13,
                    0.15, 0.16, 0.17, 0.18, 0.19)

scenario_V <- list(
  scenario_id = "V",
  manuscript_label = "Scenario V",
  scenario_type = "alternative",
  description = "Alternative scenario with narrower or more localized benefit regions, where a fixed cutoff may be inadequate.",
  prob_treat_D = prob_treat_D_V,
  prob_control = prob_control_V,
  prob_treat_I = prob_treat_I_V,
  prob_table = data.frame(
    Decipher_Midpoint = bin_midpoints,
    Prob_Treat_D = prob_treat_D_V,
    Prob_Control = prob_control_V,
    Prob_Treat_I = prob_treat_I_V
  )
)





# ============================================================
# Scenario VI: alternative scenario
# Purpose:
#   Used to assess performance of the flexible design under an
#   alternative with strong departure from a single fixed-cutoff
#   structure.
# ============================================================
prob_treat_D_VI <- c(0.070, 0.085, 0.103, 0.124, 0.162,
                    0.26, 0.341, 0.43, 0.43, 0.43)

prob_control_VI <- c(0.070, 0.085, 0.103, 0.124, 0.149,
                    0.179, 0.214, 0.255, 0.304, 0.362)

prob_treat_I_VI <- c(0.070, 0.085, 0.103, 0.0842, 0.1023,
                    0.1242, 0.1503, 0.1819, 0.2217, 0.2696)

scenario_VI <- list(
  scenario_id = "VI",
  manuscript_label = "Scenario VI",
  scenario_type = "alternative",
  description = "Alternative scenario representing a setting where a fixed cutoff is clearly inadequate.",
  prob_treat_D = prob_treat_D_VI,
  prob_control = prob_control_VI,
  prob_treat_I = prob_treat_I_VI,
  prob_table = data.frame(
    Decipher_Midpoint = bin_midpoints,
    Prob_Treat_D = prob_treat_D_VI,
    Prob_Control = prob_control_VI,
    Prob_Treat_I = prob_treat_I_VI
  )
)



# ============================================================
# Combined list of manuscript scenarios
# ============================================================
scenario_list <- list(
  III = scenario_III,
  IV = scenario_IV,
  V = scenario_V,
  VI = scenario_VI
)














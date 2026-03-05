# ============================================================
# File: 95_figures.R
# Purpose:
#   Generate manuscript figures for threshold-selection
#   frequencies and cumulative inclusion percentages for the
#   FTB design simulation study.
#
# Notes:
#   This file assumes that 90_scenario_results.R has already
#   been sourced, so that the plot data objects already exist.
# ============================================================


# ============================================================
# Helper: common legend dummy data
# ============================================================
star_legend_df <- data.frame(
  x = c(-Inf, -Inf),
  y = c(-Inf, -Inf),
  star_label = c("True NI Level", "True Sup Level")
)


# ============================================================
# Helper function: threshold-selection frequency plot
# ============================================================
make_selection_plot <- function(plot_dat,
                                true_sup_bins,
                                true_ni_bins,
                                title_label = "",
                                x_lab = "Biomarker Level Midpoint",
                                y_lab = "Percent of simulations selecting threshold") {
  
  star_y_sup <- max(plot_dat$pct, na.rm = TRUE) + 7
  star_y_ni  <- max(plot_dat$pct, na.rm = TRUE) + 4
  
  ggplot(
    plot_dat,
    aes(x = factor(threshold, levels = all_bins),
        y = pct,
        fill = trial)
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(
      aes(label = ifelse(pct == 0, "", sprintf("%.1f%%", pct))),
      position = position_dodge(width = 0.8),
      vjust = -0.25,
      size = 2
    ) +
    scale_fill_manual(
      values = c("Superiority" = "#619CFF",
                 "Non-inferiority" = "#7CAE00")
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      title = title_label
    ) +
    scale_x_discrete(labels = paste("Level", 1:10)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      plot.title = element_text(size = 9),
      plot.caption = element_text(size = 6),
      legend.title = element_blank(),
      legend.box = "vertical"
    ) +
    geom_point(
      data = data.frame(threshold = factor(true_sup_bins, levels = all_bins)),
      aes(x = threshold, y = star_y_sup),
      shape = 8, color = "#619CFF", size = 1.5, inherit.aes = FALSE
    ) +
    geom_point(
      data = data.frame(threshold = factor(true_ni_bins, levels = all_bins)),
      aes(x = threshold, y = star_y_ni),
      shape = 8, color = "#7CAE00", size = 1.5, inherit.aes = FALSE
    ) +
    geom_point(
      data = star_legend_df,
      aes(x = x, y = y, shape = star_label, color = star_label),
      size = 1,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      name = NULL,
      values = c("True NI Level" = "#7CAE00",
                 "True Sup Level" = "#619CFF")
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("True NI Level" = 8,
                 "True Sup Level" = 8),
      guide = "none"
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(shape = 8, size = 2)
      ),
      fill = guide_legend(
        order = 2,
        override.aes = list(size = 2)
      )
    )
}


# ============================================================
# Helper function: cumulative inclusion plot
# ============================================================
make_cumulative_plot <- function(plot_dat_cum,
                                 true_sup_bins,
                                 true_ni_bins,
                                 title_label = "",
                                 x_lab = "Biomarker Level Midpoint",
                                 y_lab = "Cumulative % of simulations including level") {
  
  star_y_sup <- max(plot_dat_cum$cum_pct, na.rm = TRUE) + 8
  star_y_ni  <- max(plot_dat_cum$cum_pct, na.rm = TRUE) + 5
  
  ggplot(
    plot_dat_cum,
    aes(x = factor(threshold, levels = all_bins),
        y = cum_pct,
        fill = trial)
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, na.rm = TRUE) +
    geom_text(
      aes(label = ifelse(is.na(cum_pct) | cum_pct == 0, "", sprintf("%.1f%%", cum_pct))),
      position = position_dodge(width = 0.8),
      vjust = -0.25,
      size = 2,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = c("Superiority" = "#619CFF",
                 "Non-inferiority" = "#7CAE00")
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      title = title_label
    ) +
    scale_x_discrete(labels = paste("Level", 1:10)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      plot.title = element_text(size = 9),
      plot.caption = element_text(size = 6),
      legend.title = element_blank(),
      legend.box = "vertical"
    ) +
    geom_point(
      data = data.frame(threshold = factor(true_sup_bins, levels = all_bins)),
      aes(x = threshold, y = star_y_sup),
      shape = 8, color = "#619CFF", size = 1.5, inherit.aes = FALSE
    ) +
    geom_point(
      data = data.frame(threshold = factor(true_ni_bins, levels = all_bins)),
      aes(x = threshold, y = star_y_ni),
      shape = 8, color = "#7CAE00", size = 1.5, inherit.aes = FALSE
    ) +
    geom_point(
      data = star_legend_df,
      aes(x = x, y = y, shape = star_label, color = star_label),
      size = 1,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      name = NULL,
      values = c("True NI Level" = "#7CAE00",
                 "True Sup Level" = "#619CFF")
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("True NI Level" = 8,
                 "True Sup Level" = 8),
      guide = "none"
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(shape = 8, size = 2)
      ),
      fill = guide_legend(
        order = 2,
        override.aes = list(size = 2)
      )
    )
}


# ============================================================
# Scenario D: benchmark sample size
# ============================================================
selection_ha_benchmark <- make_selection_plot(
  plot_dat      = plot_dat_ha_benchmark,
  true_sup_bins = ha_benchmark_true_sup_bins,
  true_ni_bins  = ha_benchmark_true_ni_bins,
  title_label   = " "
)

selection_cum_ha_benchmark <- make_cumulative_plot(
  plot_dat_cum  = plot_dat_cum_ha_benchmark,
  true_sup_bins = ha_benchmark_true_sup_bins,
  true_ni_bins  = ha_benchmark_true_ni_bins,
  title_label   = "A"
)


# ============================================================
# Scenario D: targeted-power sample size
# ============================================================
selection_ha_targeted <- make_selection_plot(
  plot_dat      = plot_dat_ha_targeted,
  true_sup_bins = ha_targeted_true_sup_bins,
  true_ni_bins  = ha_targeted_true_ni_bins,
  title_label   = "A"
)

selection_cum_ha_targeted <- make_cumulative_plot(
  plot_dat_cum  = plot_dat_cum_ha_targeted,
  true_sup_bins = ha_targeted_true_sup_bins,
  true_ni_bins  = ha_targeted_true_ni_bins,
  title_label   = "A"
)


# ============================================================
# Scenario E
# ============================================================
selection_case1 <- make_selection_plot(
  plot_dat      = plot_dat_case1,
  true_sup_bins = case1_true_sup_bins,
  true_ni_bins  = case1_true_ni_bins,
  title_label   = "B",
  x_lab         = "Biomarker Level Midpoint"
)

selection_cum_case1 <- make_cumulative_plot(
  plot_dat_cum  = plot_dat_cum_case1,
  true_sup_bins = case1_true_sup_bins,
  true_ni_bins  = case1_true_ni_bins,
  title_label   = "B",
  x_lab         = "Biomarker Level Midpoint"
)


# ============================================================
# Scenario F
# ============================================================
selection_case2 <- make_selection_plot(
  plot_dat      = plot_dat_case2,
  true_sup_bins = case2_true_sup_bins,
  true_ni_bins  = case2_true_ni_bins,
  title_label   = "C",
  x_lab         = "Biomarker Level Midpoint"
)

selection_cum_case2 <- make_cumulative_plot(
  plot_dat_cum  = plot_dat_cum_case2,
  true_sup_bins = case2_true_sup_bins,
  true_ni_bins  = case2_true_ni_bins,
  title_label   = "C",
  x_lab         = "Biomarker Level Midpoint",
  y_lab         = "Cumulative % of simulations including bin"
)


# ============================================================
# Optional ggsave examples
# Uncomment if you want automatic file export
# ============================================================

# ggsave("figures/selection_ha_benchmark.png",
#        selection_ha_benchmark, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_cum_ha_benchmark.png",
#        selection_cum_ha_benchmark, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_ha_targeted.png",
#        selection_ha_targeted, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_cum_ha_targeted.png",
#        selection_cum_ha_targeted, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_case1.png",
#        selection_case1, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_cum_case1.png",
#        selection_cum_case1, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_case2.png",
#        selection_case2, width = 8, height = 3.5, dpi = 300)

# ggsave("figures/selection_cum_case2.png",
#        selection_cum_case2, width = 8, height = 3.5, dpi = 300)
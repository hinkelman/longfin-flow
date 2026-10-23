# Summarize the two extension studies:
#  dlm_variants  why does the paper's DLM shrink the flow effect? (truth = growth)
#  level_truth   how do the methods do when the DLM's own model is the truth?

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dir.create(file.path("output", "figures"), showWarnings = FALSE)
scenario_labels <- c(null = "Null (no effect)", constant = "Constant",
                     decline = "Linear decline", step_1987 = "Step at 1987",
                     step_2002 = "Step at 2002")
variant_labels <- c(A_paper = "A. Paper DLM (smoothed states, level, RW intercept)",
                    B_true_states = "B. True states instead of smoothed",
                    C_growth = "C. Flow drives growth, not level",
                    D_growth_no_int = "D. Growth, no RW intercept",
                    E_one_stage = "E. One-stage MARSS, time-varying C")
variant_cols <- c(A_paper = "#2a78d6", B_true_states = "#eb6834", C_growth = "#1baf7a",
                  D_growth_no_int = "#eda100", E_one_stage = "#e87ba4")

summarize_variants <- function(v) {
  v |>
    mutate(tv = q_slope > 1e-8) |>
    group_by(truth, scenario, b, variant) |>
    summarize(n_fits = n(),
              true_mean = mean(true_mean),
              mean_slope = mean(mean_slope),
              ratio_to_truth = ifelse(true_mean == 0, NA, mean_slope / true_mean),
              true_change = mean(true_first - true_last),
              est_change = mean(first - last),
              flags_time_varying = mean(tv),
              ci_coverage = mean(cover),
              rmse = mean(rmse),
              converged = mean(conv == 0),
              .groups = "drop") |>
    mutate(scenario = factor(scenario, names(scenario_labels))) |>
    arrange(truth, b, scenario, variant)
}

# Study 1: decomposing the DLM's shrinkage ---------------------------------------------
results_dir <- file.path("output", "results")
dv <- read.csv(file.path(results_dir, "dlm_variants_variants.csv"))
tab_variants <- summarize_variants(dv)
write.csv(tab_variants, file.path("output", "ext_dlm_variants.csv"), row.names = FALSE)
cat("DLM variants, truth = growth (the paper's MARSS model):\n")
print(as.data.frame(tab_variants |> select(scenario, b, variant, ratio_to_truth, true_change,
                                           est_change, flags_time_varying, ci_coverage, converged)),
      digits = 2)

paths_fig <- function(paths, b_val, title, subtitle) {
  med <- paths |>
    filter(b == b_val) |>
    mutate(scenario = factor(scenario, names(scenario_labels))) |>
    rename(slope = median)
  truth <- distinct(med, scenario, year, true_beta)
  ggplot(med, aes(year)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_line(aes(y = slope, colour = variant), linewidth = 0.8) +
    geom_step(data = truth, aes(y = true_beta), colour = "#0b0b0b", linewidth = 0.8,
              direction = "mid") +
    facet_wrap(~ scenario, nrow = 1, labeller = as_labeller(scenario_labels)) +
    scale_colour_manual(NULL, values = variant_cols, labels = variant_labels) +
    scale_x_continuous(breaks = c(1985, 2000, 2015)) +
    guides(colour = guide_legend(ncol = 2)) +
    labs(x = NULL, y = "Flow effect (per SD of outflow)", title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title.position = "plot")
}
ggsave(file.path("output", "figures", "ext_dlm_variants_b024.png"),
       paths_fig(filter(read.csv(file.path(results_dir, "dlm_variants_variant_path_quantiles.csv")),
                        scenario != "null"), 0.24,
                 "Which modeling choice shrinks the DLM's flow effect? (starting effect b = 0.24)",
                 "Median slope over 100 replicates x 4 states; black = true effect. B-D use the simulated true states (diagnostic only)"),
       width = 13, height = 5, dpi = 200, bg = "white")

ratio_plot <- tab_variants |>
  filter(scenario != "null") |>
  mutate(b_lab = sprintf("b = %.2f", b)) |>
  ggplot(aes(ratio_to_truth, scenario, colour = variant)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "grey40") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  facet_wrap(~ b_lab) +
  scale_y_discrete(labels = scenario_labels, limits = rev) +
  scale_colour_manual(NULL, values = variant_cols, labels = variant_labels) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(x = "Mean estimated effect / mean true effect (1 = unbiased)", y = NULL,
       title = "How much of the true flow effect each DLM variant recovers") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title.position = "plot")
ggsave(file.path("output", "figures", "ext_dlm_variants_ratio.png"), ratio_plot,
       width = 11, height = 5.5, dpi = 200, bg = "white")

# Study 2: the DLM's own model as truth -------------------------------------------------
s <- read.csv(file.path(results_dir, "level_truth_summary.csv")) |> mutate(scenario = factor(scenario, names(scenario_labels)),
                          dlm_tv = dlm_q_slope > 1e-8,
                          marss_lo = marss_C - 1.96 * marss_se,
                          marss_hi = marss_C + 1.96 * marss_se)
reps <- distinct(s, scenario, b, rep, dAICc_era)
# pooled-null threshold computed in R/05_followup_checks.R
threshold <- read.csv(file.path("output", "followup_era_calibrated_power.csv"))$calibrated_threshold[1]
tab_level <- s |>
  group_by(scenario, b) |>
  summarize(true_mean = mean(true_mean),
            marss_C_mean = mean(marss_C),
            marss_cover_true_mean = mean(marss_lo <= true_mean & true_mean <= marss_hi, na.rm = TRUE),
            dlm_mean_slope = mean(dlm_mean),
            dlm_ratio_to_truth = ifelse(true_mean == 0, NA, mean(dlm_mean) / true_mean),
            dlm_flags_time_varying = mean(dlm_tv),
            dlm_pointwise_cover = mean(dlm_cover),
            .groups = "drop") |>
  left_join(reps |> group_by(scenario, b) |>
              summarize(era_flags_naive = mean(dAICc_era < -2),
                        era_flags_calibrated = mean(dAICc_era < threshold), .groups = "drop"),
            by = c("scenario", "b")) |>
  mutate(era_calibrated_threshold = threshold)
write.csv(tab_level, file.path("output", "ext_level_truth.csv"), row.names = FALSE)
cat("\nTruth = level (the paper's DLM model): main methods\n")
cat("Pooled null-calibrated era threshold: dAICc <", round(threshold, 2), "\n")
print(as.data.frame(tab_level), digits = 2)

tab_level_variants <- summarize_variants(read.csv(file.path(results_dir, "level_truth_variants.csv")))
write.csv(tab_level_variants, file.path("output", "ext_level_truth_variants.csv"), row.names = FALSE)
cat("\nTruth = level: DLM variants\n")
print(as.data.frame(tab_level_variants |> select(scenario, b, variant, ratio_to_truth, true_change,
                                                 est_change, flags_time_varying, ci_coverage)),
      digits = 2)

ggsave(file.path("output", "figures", "ext_level_truth_variants_b024.png"),
       paths_fig(filter(read.csv(file.path(results_dir, "level_truth_variant_path_quantiles.csv")),
                        scenario != "null"), 0.24,
                 "When flow drives the level (the DLM's own model), starting effect b = 0.24",
                 "Median slope over 100 replicates x 4 states; black = true effect. B-D use the simulated true states (diagnostic only)"),
       width = 13, height = 5, dpi = 200, bg = "white")

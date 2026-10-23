# Summarize the two extension studies:
#  dlm_variants  why does the paper's DLM shrink the flow effect? (truth = growth)
#  level_truth   how do the methods do when the DLM's own model is the truth?

source("R/plot_helpers.R")

VARIANT_LABELS <- c(A_paper = "A. Paper DLM (smoothed states, level, RW intercept)",
                    B_true_states = "B. True states instead of smoothed",
                    C_growth = "C. Flow drives growth, not level",
                    D_growth_no_int = "D. Growth, no RW intercept",
                    E_one_stage = "E. One-stage MARSS, time-varying C")
VARIANT_COLS <- c(A_paper = COL[["blue"]], B_true_states = COL[["orange"]],
                  C_growth = COL[["aqua"]], D_growth_no_int = COL[["yellow"]],
                  E_one_stage = COL[["magenta"]])
PATHS_SUBTITLE <- paste("Median slope over 100 replicates x 4 states; black = true effect.",
                        "B-D use the simulated true states (diagnostic only)")

summarize_variants <- function(v) {
  v |>
    mutate(tv = q_slope > Q_TV_THRESHOLD) |>
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
    mutate(scenario = factor_scenario(scenario)) |>
    arrange(truth, b, scenario, variant)
}

# Median slope paths of each variant against the true effect, at one effect size
paths_fig <- function(study, b_val, title) {
  med <- read_result(study, "variant_path_quantiles") |>
    filter(scenario != "null", b == b_val) |>
    mutate(scenario = factor_scenario(scenario))
  truth <- distinct(med, scenario, year, true_beta)
  ggplot(med, aes(year)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_line(aes(y = median, colour = variant), linewidth = 0.8) +
    geom_step(data = truth, aes(y = true_beta), colour = COL[["truth"]], linewidth = 0.8,
              direction = "mid") +
    facet_wrap(~ scenario, nrow = 1, labeller = as_labeller(SCENARIO_LABELS)) +
    scale_colour_manual(NULL, values = VARIANT_COLS, labels = VARIANT_LABELS) +
    scale_x_continuous(breaks = c(1985, 2000, 2015)) +
    guides(colour = guide_legend(ncol = 2)) +
    labs(x = NULL, y = "Flow effect (per SD of outflow)", title = title, subtitle = PATHS_SUBTITLE) +
    theme_talk() +
    theme(panel.grid.minor = element_blank())
}

# Study 1: decomposing the DLM's shrinkage ---------------------------------------------
tab_variants <- summarize_variants(read_result("dlm_variants", "variants"))
write.csv(tab_variants, file.path("output", "ext_dlm_variants.csv"), row.names = FALSE)
cat("DLM variants, truth = growth (the paper's MARSS model):\n")
print(as.data.frame(tab_variants |> select(scenario, b, variant, ratio_to_truth, true_change,
                                           est_change, flags_time_varying, ci_coverage, converged)),
      digits = 2)

save_figure("ext_dlm_variants_b024.png",
            paths_fig("dlm_variants", 0.24,
                      "Which modeling choice shrinks the DLM's flow effect? (starting effect b = 0.24)"),
            width = 13, height = 5)

ratio_plot <- tab_variants |>
  filter(scenario != "null") |>
  mutate(b_lab = b_label(scenario, b)) |>
  ggplot(aes(ratio_to_truth, scenario, colour = variant)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "grey40") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  facet_wrap(~ b_lab) +
  scale_y_discrete(labels = SCENARIO_LABELS, limits = rev) +
  scale_colour_manual(NULL, values = VARIANT_COLS, labels = VARIANT_LABELS) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(x = "Mean estimated effect / mean true effect (1 = unbiased)", y = NULL,
       title = "How much of the true flow effect each DLM variant recovers") +
  theme_talk()
save_figure("ext_dlm_variants_ratio.png", ratio_plot, width = 11, height = 5.5)

# Study 2: the DLM's own model as truth -------------------------------------------------
s <- read_result("level_truth", "summary") |>
  mutate(scenario = factor_scenario(scenario),
         dlm_tv = dlm_q_slope > Q_TV_THRESHOLD,
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

tab_level_variants <- summarize_variants(read_result("level_truth", "variants"))
write.csv(tab_level_variants, file.path("output", "ext_level_truth_variants.csv"), row.names = FALSE)
cat("\nTruth = level: DLM variants\n")
print(as.data.frame(tab_level_variants |> select(scenario, b, variant, ratio_to_truth, true_change,
                                                 est_change, flags_time_varying, ci_coverage)),
      digits = 2)

save_figure("ext_level_truth_variants_b024.png",
            paths_fig("level_truth", 0.24,
                      "When flow drives the level (the DLM's own model), starting effect b = 0.24"),
            width = 13, height = 5)

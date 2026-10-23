# Summarize the structure and interval checks (studies extra_growth, extra_level).
#  structure: how often does AICc pick the true structure (flow on growth vs level)?
#  intervals: do one-stage time-varying slopes with a prior on x0 give honest 95%
#             intervals, compared with x0 estimated, and the paper's DLM with a prior?

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

results_dir <- file.path("output", "results")
scenario_labels <- c(null = "Null (no effect)", constant = "Constant",
                     decline = "Linear decline", step_1987 = "Step at 1987",
                     step_2002 = "Step at 2002")
truth_labels <- c(growth = "Truth: flow drives growth", level = "Truth: flow drives level")

read_study <- function(name) {
  bind_rows(lapply(c("extra_growth", "extra_level"), function(s)
    read.csv(file.path(results_dir, paste0(s, "_", name, ".csv")))))
}

# Structure ---------------------------------------------------------------------------
st <- read_study("structure") |>
  mutate(scenario = factor(scenario, names(scenario_labels)))
reps <- distinct(st, truth, scenario, b, rep, AICc_stage1, AICc_C, AICc_D, AICc_era_C, AICc_era_D)

# null-calibrated thresholds for "era beats static" within each form, pooling the
# 200 null replicates (identical process under both truths)
nulls <- filter(reps, scenario == "null")
thr_era_C <- unname(quantile(nulls$AICc_era_C - nulls$AICc_C, 0.05))
thr_era_D <- unname(quantile(nulls$AICc_era_D - nulls$AICc_D, 0.05))

tab_structure <- reps |>
  group_by(truth, scenario, b) |>
  summarize(n_reps = n(),
            picks_level = mean(AICc_D < AICc_C),
            picks_level_by_2 = mean(AICc_D < AICc_C - 2),
            picks_growth_by_2 = mean(AICc_C < AICc_D - 2),
            level_beats_none_by_2 = mean(AICc_D < AICc_stage1 - 2),
            growth_beats_none_by_2 = mean(AICc_C < AICc_stage1 - 2),
            era_C_flags_calibrated = mean(AICc_era_C - AICc_C < thr_era_C),
            era_D_flags_calibrated = mean(AICc_era_D - AICc_D < thr_era_D),
            .groups = "drop") |>
  mutate(threshold_era_C = thr_era_C, threshold_era_D = thr_era_D)

# true means come from the interval metrics (same replicates, scored 1981-2020)
v <- read_study("variants") |> mutate(scenario = factor(scenario, names(scenario_labels)))
true_means <- distinct(v, truth, scenario, b, rep, state, true_mean)
tab_static <- st |>
  left_join(true_means, by = c("truth", "scenario", "b", "rep", "state")) |>
  group_by(truth, scenario, b) |>
  summarize(true_mean = mean(true_mean),
            C_mean = mean(C_est), D_mean = mean(D_est),
            C_cover = mean(abs(C_est - true_mean) <= 1.96 * C_se, na.rm = TRUE),
            D_cover = mean(abs(D_est - true_mean) <= 1.96 * D_se, na.rm = TRUE),
            .groups = "drop")

structure_out <- left_join(tab_structure, tab_static, by = c("truth", "scenario", "b"))
write.csv(structure_out, file.path("output", "check_structure.csv"), row.names = FALSE)
cat("Null-calibrated era thresholds: growth form dAICc <", round(thr_era_C, 2),
    "; level form dAICc <", round(thr_era_D, 2), "\n\n")
cat("Structure choice (share of replicates):\n")
print(as.data.frame(structure_out |> select(truth, scenario, b, picks_level, picks_level_by_2,
                                            picks_growth_by_2, level_beats_none_by_2,
                                            growth_beats_none_by_2)), digits = 2, row.names = FALSE)
cat("\nStatic effect estimates and 95% CI coverage of the true mean:\n")
print(as.data.frame(structure_out |> select(truth, scenario, b, true_mean, C_mean, C_cover,
                                            D_mean, D_cover, era_C_flags_calibrated,
                                            era_D_flags_calibrated)), digits = 2, row.names = FALSE)

p_structure <- structure_out |>
  mutate(b_lab = ifelse(scenario == "null", "no effect", sprintf("b = %.2f", b)),
         truth = truth_labels[truth]) |>
  select(truth, scenario, b_lab, `Level preferred (dAICc > 2)` = picks_level_by_2,
         `Growth preferred (dAICc > 2)` = picks_growth_by_2) |>
  pivot_longer(-c(truth, scenario, b_lab), names_to = "choice", values_to = "rate") |>
  ggplot(aes(scenario, rate, fill = choice)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_grid(truth ~ b_lab, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = scenario_labels) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(NULL, values = c("Level preferred (dAICc > 2)" = "#eb6834",
                                     "Growth preferred (dAICc > 2)" = "#2a78d6")) +
  labs(x = NULL, y = "Share of replicates",
       title = "Can AICc tell whether flow drives growth or level?",
       subtitle = "Static flow effect fit both ways to the same simulated station data") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1), plot.title.position = "plot")
ggsave(file.path("output", "figures", "check_structure.png"), p_structure,
       width = 11, height = 6, dpi = 200, bg = "white")

# Intervals ---------------------------------------------------------------------------
interval_labels <- c(F_tv_matched = "One-stage, x0 estimated",
                     G_tv_matched_prior = "One-stage, prior on x0",
                     H_paper_dlm_prior = "Paper DLM, prior on x0")
tab_intervals <- v |>
  mutate(tv = q_slope > 1e-8) |>
  group_by(truth, scenario, b, variant) |>
  summarize(ratio_to_truth = ifelse(mean(true_mean) == 0, NA, mean(mean_slope) / mean(true_mean)),
            flags_time_varying = mean(tv),
            ci_coverage = mean(cover),
            mean_se = mean(mean_se),
            rmse = mean(rmse),
            converged = mean(conv == 0),
            .groups = "drop")
write.csv(tab_intervals, file.path("output", "check_intervals.csv"), row.names = FALSE)
cat("\nInterval check:\n")
print(as.data.frame(tab_intervals |> select(truth, scenario, b, variant, ratio_to_truth,
                                            flags_time_varying, ci_coverage, mean_se, converged)),
      digits = 2, row.names = FALSE)

p_intervals <- tab_intervals |>
  mutate(b_lab = ifelse(scenario == "null", "no effect", sprintf("b = %.2f", b)),
         truth = truth_labels[truth], variant = interval_labels[variant]) |>
  ggplot(aes(ci_coverage, scenario, colour = variant)) +
  geom_vline(xintercept = 0.95, linetype = "22", colour = "grey40") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  facet_grid(truth ~ b_lab, scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = scenario_labels, limits = rev) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(NULL, values = c("#6250d6", "#1baf7a", "#2a78d6")) +
  labs(x = "Share of years where the 95% interval covers the true effect (target 95%)", y = NULL,
       title = "Does a prior on x0 fix the slope intervals?") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title.position = "plot")
ggsave(file.path("output", "figures", "check_intervals.png"), p_intervals,
       width = 11, height = 6, dpi = 200, bg = "white")

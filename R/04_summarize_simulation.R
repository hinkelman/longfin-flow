# Summarize the simulation study: recovery tables and figures for the talk.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

res <- readRDS(file.path("output", "sim_results.rds"))
dir.create(file.path("output", "figures"), showWarnings = FALSE)

scenario_labels <- c(null = "Null (no effect)", constant = "Constant",
                     decline = "Linear decline", step_1987 = "Step at 1987",
                     step_2002 = "Step at 2002")
col_truth <- "#0b0b0b"
col_marss <- "#2a78d6"
col_dlm <- "#eb6834"
col_era <- "#1baf7a"

s <- res$summary |>
  mutate(scenario = factor(scenario, names(scenario_labels)),
         marss_lo = marss_C - 1.96 * marss_se,
         marss_hi = marss_C + 1.96 * marss_se,
         dlm_tv = dlm_q_slope > 1e-8)

# Recovery table --------------------------------------------------------------
# State-level metrics pool the 4 habitat states; the era AICc test is per replicate
tab_state <- s |>
  group_by(scenario, b) |>
  summarize(
    n_fits = n(),
    true_mean = mean(true_mean),
    marss_C_mean = mean(marss_C),
    marss_cover_true_mean = mean(marss_lo <= true_mean & true_mean <= marss_hi, na.rm = TRUE),
    marss_sig_pos = mean(marss_lo > 0, na.rm = TRUE),
    dlm_mean_slope = mean(dlm_mean),
    dlm_ratio_to_truth = mean(dlm_mean) / mean(true_mean),
    dlm_flags_time_varying = mean(dlm_tv),
    dlm_first_minus_last = mean(dlm_first - dlm_last),
    true_first_minus_last = mean(true_first - true_last),
    dlm_pointwise_cover = mean(dlm_cover),
    dlm_static_cover = mean(dlm_cover[!dlm_tv]),
    sign_agree_marss_dlm = mean(sign(marss_C) == sign(dlm_mean)),
    .groups = "drop")

tab_rep <- s |>
  distinct(scenario, b, rep, dAICc_era) |>
  group_by(scenario, b) |>
  summarize(era_preferred_dAICc_lt_minus2 = mean(dAICc_era < -2), .groups = "drop")

tab_era <- s |>
  group_by(scenario, b) |>
  summarize(across(c(era_pre1987, era_1987_2002, era_post2002,
                     true_pre1987, true_1987_2002, true_post2002), mean),
            .groups = "drop")

recovery <- tab_state |> left_join(tab_rep, by = c("scenario", "b")) |>
  left_join(tab_era, by = c("scenario", "b"))
write.csv(recovery, file.path("output", "sim_recovery_table.csv"), row.names = FALSE)

convergence <- s |> summarize(across(starts_with("conv_"), ~ mean(.x == 0)))
cat("Fraction of fits converged:\n"); print(convergence)
cat("\nRecovery summary:\n")
print(as.data.frame(recovery |> select(scenario, b, true_mean, marss_C_mean, marss_cover_true_mean,
                                       dlm_mean_slope, dlm_ratio_to_truth, dlm_flags_time_varying,
                                       dlm_pointwise_cover, era_preferred_dAICc_lt_minus2,
                                       sign_agree_marss_dlm)), digits = 2)
cat("\nEra estimates vs truth:\n"); print(as.data.frame(tab_era), digits = 2)

# Figure 1: what each method recovers over time ----------------------------------
paths <- res$paths |> mutate(scenario = factor(scenario, names(scenario_labels)))
dlm_band <- paths |>
  group_by(scenario, b, year) |>
  summarize(true_beta = first(true_beta),
            lo = quantile(dlm_slope, 0.1), mid = median(dlm_slope), hi = quantile(dlm_slope, 0.9),
            .groups = "drop")

marss_line <- s |> group_by(scenario, b) |>
  summarize(marss_mid = median(marss_C), .groups = "drop")

era_steps <- s |>
  group_by(scenario, b) |>
  summarize(pre1987 = median(era_pre1987), `1987_2002` = median(era_1987_2002),
            post2002 = median(era_post2002), .groups = "drop") |>
  pivot_longer(c(pre1987, `1987_2002`, post2002), names_to = "era", values_to = "est") |>
  mutate(start = c(pre1987 = 1980, `1987_2002` = 1987, post2002 = 2003)[era],
         end = c(pre1987 = 1986, `1987_2002` = 2002, post2002 = 2020)[era])

recovery_plot <- function(b_val) {
  ggplot(filter(dlm_band, b == b_val), aes(year)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col_dlm, alpha = 0.2) +
    geom_line(aes(y = mid, colour = "DLM slope (paper)"), linewidth = 0.8) +
    geom_hline(data = filter(marss_line, b == b_val),
               aes(yintercept = marss_mid, colour = "MARSS C (paper)"), linewidth = 0.8,
               linetype = "22") +
    geom_segment(data = filter(era_steps, b == b_val),
                 aes(x = start, xend = end, y = est, yend = est, colour = "Era MARSS"),
                 linewidth = 1.1) +
    geom_step(aes(y = true_beta, colour = "True effect"), linewidth = 0.8, direction = "mid") +
    facet_wrap(~ scenario, nrow = 1, labeller = as_labeller(scenario_labels)) +
    scale_colour_manual(NULL, values = c("True effect" = col_truth, "MARSS C (paper)" = col_marss,
                                         "DLM slope (paper)" = col_dlm, "Era MARSS" = col_era),
                        breaks = c("True effect", "MARSS C (paper)", "Era MARSS", "DLM slope (paper)")) +
    labs(x = NULL, y = "Flow effect (per SD of outflow)",
         title = sprintf("Recovered flow effects, starting effect b = %.2f", b_val),
         subtitle = "Medians over 100 replicates x 4 habitat states; orange band = 10th-90th percentile of DLM slopes") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title.position = "plot")
}
ggsave(file.path("output", "figures", "sim_recovery_b008.png"), recovery_plot(0.08),
       width = 13, height = 4.5, dpi = 200, bg = "white")
ggsave(file.path("output", "figures", "sim_recovery_b024.png"), recovery_plot(0.24),
       width = 13, height = 4.5, dpi = 200, bg = "white")

# Figure 2: does each method say "the effect changed"? --------------------------
detect <- recovery |>
  select(scenario, b, `DLM: slope variance > 0` = dlm_flags_time_varying,
         `Era MARSS: dAICc < -2` = era_preferred_dAICc_lt_minus2) |>
  pivot_longer(-c(scenario, b), names_to = "method", values_to = "rate") |>
  mutate(b_lab = ifelse(scenario == "null", "no effect", sprintf("b = %.2f", b)))

p_detect <- ggplot(detect, aes(scenario, rate, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = scales::percent(rate, accuracy = 1)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3.2, colour = "grey20") +
  facet_wrap(~ b_lab, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = scenario_labels) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08)) +
  scale_fill_manual(NULL, values = c(col_dlm, col_era)) +
  labs(x = NULL, y = "Share of fits flagging a time-varying effect",
       title = "How often each method flags a change in the flow effect") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1), plot.title.position = "plot")
ggsave(file.path("output", "figures", "sim_detection.png"), p_detect,
       width = 11, height = 4.8, dpi = 200, bg = "white")

# Figure 3: DLM interval coverage --------------------------------------------------
p_cover <- s |>
  mutate(type = ifelse(dlm_tv, "Slope allowed to vary (q > 0)", "Slope collapsed to static (q ~ 0)")) |>
  ggplot(aes(dlm_mean_se, dlm_cover, colour = type)) +
  geom_point(alpha = 0.35, size = 1.4) +
  geom_hline(yintercept = 0.95, linetype = "22", colour = "grey40") +
  scale_x_log10() +
  scale_colour_manual(NULL, values = c(col_dlm, "#6250d6")) +
  labs(x = "Mean DLM slope SE (log scale)", y = "Share of years where 95% CI covers true effect",
       title = "DLM 95% intervals: nominal vs actual coverage") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title.position = "plot")
ggsave(file.path("output", "figures", "sim_dlm_coverage.png"), p_cover,
       width = 7, height = 5, dpi = 200, bg = "white")

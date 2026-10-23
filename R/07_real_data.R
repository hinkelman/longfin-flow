# Fit every flow-model form to the real Longfin Smelt station data (age-0 and
# age-1+, contemporary Jan-Jun outflow), using the authors' station selection,
# Z6 habitat grouping, and gear-specific observation variances:
#   none           no flow effect (the paper's stage-1 model)
#   static         one flow effect per state, as growth (C, the paper's MARSS) or level (D)
#   era            separate effects for 1980-86, 1987-2002, 2003-20, as growth or level
#   time-varying   random-walk flow effect, as growth or level (one-stage models)
# Time-varying models are fit twice: with x0 estimated (likelihood comparable with
# the other models, used for AICc) and with a prior on x0 (used for intervals).
# Run after R/00_run_original.R 1 and R/02_sim_setup.R.

suppressMessages({
  library(dplyr)
  library(ggplot2)
})
source("R/sim_functions.R")

age0 <- readRDS(file.path("output", "sim_setup.rds"))
# script 1 ends with the age-1+ section, so its saved environment holds age-1+ objects
e1 <- readRDS(file.path("output", "original_run", "script1", "env.rds"))
age1 <- list(years = age0$years, flow_z = age0$flow_z, y_obs = e1$MWT_OT_matrix2,
             Z = e1$Z6, R_model = e1$R.model, gear = e1$z_matrix$Gear,
             x0 = setNames(coef(e1$z6)$x0[, 1], rownames(e1$z6$states)))
setups <- list(`age-0` = age0, `age-1+` = age1)

state_labels <- c(`1 MWT` = "Channel, midwater", `1 OT` = "Channel, otter trawl",
                  `2 MWT` = "Shoal, midwater", `2 OT` = "Shoal, otter trawl")

fit_age <- function(setup) {
  y <- setup$y_obs
  k <- length(setup$x0)
  s1 <- fit_stage1(y, setup)
  s1_inits <- list(Q = pmax(coef(s1)$Q, 1e-4), R = coef(s1)$R, x0 = coef(s1)$x0)
  refit <- function(f, form) {
    fit <- f(y, setup, inits = s1_inits, form = form)
    if (is.null(fit$par)) fit <- f(y, setup, form = form)
    fit
  }
  x0_mean <- c(coef(s1)$x0[, 1], rep(0, k))
  list(
    none = s1,
    static_growth = refit(fit_marss_c, "growth"),
    static_level = refit(fit_marss_c, "level"),
    era_growth = refit(fit_marss_era, "growth"),
    era_level = refit(fit_marss_era, "level"),
    tv_growth = fit_tv_marss(y, setup, inits = list(x0 = matrix(x0_mean, ncol = 1)), form = "growth"),
    tv_level = fit_tv_marss(y, setup, inits = list(x0 = matrix(x0_mean, ncol = 1)), form = "level"),
    tv_growth_prior = fit_tv_marss(y, setup, form = "growth", prior_var = 1, x0_mean = x0_mean),
    tv_level_prior = fit_tv_marss(y, setup, form = "level", prior_var = 1, x0_mean = x0_mean)
  )
}
fits <- lapply(setups, fit_age)

# Model comparison ---------------------------------------------------------------------
comparison <- bind_rows(lapply(names(fits), function(age) {
  f <- fits[[age]]
  ms <- c("none", "static_growth", "static_level", "era_growth", "era_level", "tv_growth", "tv_level")
  data.frame(age = age, model = ms,
             logLik = sapply(f[ms], `[[`, "logLik"),
             K = sapply(f[ms], function(x) x$num.params),
             AICc = sapply(f[ms], `[[`, "AICc"),
             converged = sapply(f[ms], function(x) x$convergence == 0))
})) |>
  group_by(age) |>
  mutate(dAICc = AICc - min(AICc)) |>
  ungroup()
write.csv(comparison, file.path("output", "real_data_model_comparison.csv"), row.names = FALSE)
cat("Model comparison (dAICc relative to the best model for each age):\n")
print(as.data.frame(comparison), digits = 4, row.names = FALSE)

# Static and era flow effects -----------------------------------------------------------
effects <- bind_rows(lapply(names(fits), function(age) {
  setup <- setups[[age]]
  bind_rows(lapply(c("static_growth", "static_level", "era_growth", "era_level"), function(m) {
    form <- sub(".*_", "", m)
    prefix <- if (form == "growth") "C." else "D."
    est <- flow_coefs(fits[[age]][[m]], setup, form)
    ct <- coef_table(fits[[age]][[m]], prefix)
    covs <- if (ncol(est) == 1) "all years" else c("1980-86", "1987-2002", "2003-20")
    terms <- if (ncol(est) == 1) paste0(prefix, rownames(est)) else
      paste0(prefix, rep(rownames(est), 3), " ", rep(c("pre1987", "1987_2002", "post2002"), each = 4))
    data.frame(age = age, model = m, state = rep(rownames(est), ncol(est)),
               period = rep(covs, each = nrow(est)), estimate = as.vector(est),
               se = if (is.null(ct)) NA else ct$std.error[match(terms, ct$term)])
  }))
})) |>
  mutate(lo95 = estimate - 1.96 * se, hi95 = estimate + 1.96 * se)
write.csv(effects, file.path("output", "real_data_flow_effects.csv"), row.names = FALSE)
cat("\nStatic flow effects (per SD of Jan-Jun outflow):\n")
print(as.data.frame(filter(effects, grepl("static", model))), digits = 3, row.names = FALSE)

# Time-varying slopes (prior versions, for intervals) ---------------------------------------
tv <- bind_rows(lapply(names(fits), function(age) {
  bind_rows(lapply(c("growth", "level"), function(form) {
    f <- fits[[age]][[paste0("tv_", form, "_prior")]]
    p <- tv_slopes(f, setups[[age]], form)
    q <- coef(f)$Q[, 1]
    data.frame(age = age, form = form, state = rep(rownames(p$slope), ncol(p$slope)),
               year = rep(setups[[age]]$years, each = nrow(p$slope)),
               slope = as.vector(p$slope), se = as.vector(p$se),
               q_slope = rep(q[paste0("q_slo_", rownames(p$slope))], ncol(p$slope)))
  }))
}))
write.csv(tv, file.path("output", "real_data_tv_slopes.csv"), row.names = FALSE)
cat("\nTime-varying slope variances (q_slope; ~0 means the effect is estimated as constant):\n")
print(as.data.frame(distinct(tv, age, form, state, q_slope)), digits = 3, row.names = FALSE)

p_tv <- tv |>
  filter(!is.na(slope)) |>
  mutate(state = factor(state_labels[state], state_labels),
         form = ifelse(form == "growth", "Flow drives growth (process eq.)",
                       "Flow drives level (observation eq.)")) |>
  ggplot(aes(year, slope, colour = form, fill = form)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_ribbon(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(age ~ state) +
  scale_colour_manual(NULL, values = c("#2a78d6", "#eb6834"), aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = c(1985, 2000, 2015)) +
  labs(x = NULL, y = "Flow effect (per SD of outflow)",
       title = "Time-varying flow effects in the real Longfin Smelt data",
       subtitle = "One-stage models with a prior on x0; bands are 95% intervals") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path("output", "figures", "real_data_tv_slopes.png"), p_tv,
       width = 12, height = 6, dpi = 200, bg = "white")

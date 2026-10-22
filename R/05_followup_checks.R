# Follow-up checks on the simulation results.
#  1. What does the static MARSS C average over? Compare it with the plain mean
#     of the true effect, a flow-weighted mean, and a flow x sampling-weighted mean.
#  2. Why does the era model flag change in 20% of null replicates? Refit the
#     worst null replicates with a BFGS polish (and a nested start) to rule out
#     incomplete optimization, then calibrate the dAICc threshold on the nulls.

source("R/sim_functions.R")
setup <- readRDS(file.path("output", "sim_setup.rds"))
s <- readRDS(file.path("output", "sim_results.rds"))$summary
yrs <- setup$years

# 1. MARSS C as a weighted average ---------------------------------------------------
w_flow <- setup$flow_z^2
w_flow_obs <- setup$flow_z^2 * colMeans(setup$mask)
weighted <- do.call(rbind, lapply(split(s, list(s$scenario, s$b), drop = TRUE), function(d) {
  beta <- beta_path(d$scenario[1], yrs, d$b[1])
  data.frame(scenario = d$scenario[1], b = d$b[1],
             true_mean = mean(beta),
             flow_weighted = sum(beta * w_flow) / sum(w_flow),
             flow_x_sampling_weighted = sum(beta * w_flow_obs) / sum(w_flow_obs),
             marss_C_mean = mean(d$marss_C),
             marss_C_mcse = sd(d$marss_C) / sqrt(nrow(d)))
}))
rownames(weighted) <- NULL
write.csv(weighted, file.path("output", "followup_marss_weighting.csv"), row.names = FALSE)
cat("Static MARSS C vs candidate averages of the true effect:\n")
print(weighted, digits = 3)
cat("\nFraction of series observed: 1980-86", round(mean(setup$mask[, yrs < 1987]), 2),
    " 1987-2002", round(mean(setup$mask[, yrs >= 1987 & yrs <= 2002]), 2),
    " 2003-20", round(mean(setup$mask[, yrs > 2002]), 2), "\n")

# 2a. Are the null false positives an optimization artifact? -------------------------
reps_dAICc <- unique(s[, c("scenario", "b", "rep", "dAICc_era")])
null_reps <- reps_dAICc[reps_dAICc$scenario == "null", ]
worst <- head(null_reps$rep[order(null_reps$dAICc_era)], 8)

polish <- function(fit, y, model) {
  b <- quiet_marss(y, model = model, inits = fit, method = "BFGS", control = list(maxit = 5000))
  if (!is.null(b$par) && b$logLik > fit$logLik) b else fit
}
refits <- parallel::mclapply(worst, function(r) {
  # null replicates are grid rows 1-100, so the seed is the one 03_run_simulation.R used
  set.seed(20261112 + r)
  sim <- simulate_data(setup, beta_path("null", yrs, 0.08))
  base <- list(Z = setup$Z, Q = "diagonal and unequal", R = setup$R_model, B = "identity", A = "zero")
  m_c <- c(base, list(c = matrix(setup$flow_z, nrow = 1), C = "unequal"))
  m_e <- c(base, list(c = era_covariates(setup), C = "unconstrained"))
  fc0 <- quiet_marss(sim$y, model = m_c, control = list(maxit = 10000))
  fe0 <- quiet_marss(sim$y, model = m_e, control = list(maxit = 10000))
  fc <- polish(fc0, sim$y, m_c)
  fe <- polish(fe0, sim$y, m_e)
  # the era model nests the static one: also start it at the static solution
  nested <- list(Q = coef(fc)$Q, R = coef(fc)$R, x0 = coef(fc)$x0,
                 C = matrix(rep(coef(fc)$C[, 1], 3), ncol = 1))
  fe2 <- polish(quiet_marss(sim$y, model = m_e, inits = nested, control = list(maxit = 10000)),
                sim$y, m_e)
  if (fe2$logLik > fe$logLik) fe <- fe2
  data.frame(rep = r,
             sim_dAICc = null_reps$dAICc_era[null_reps$rep == r],
             logLik_gain_static = fc$logLik - fc0$logLik,
             logLik_gain_era = fe$logLik - fe0$logLik,
             polished_dAICc = fe$AICc - fc$AICc,
             LR = 2 * (fe$logLik - fc$logLik))
}, mc.cores = min(8, parallel::detectCores()))
refits <- do.call(rbind, refits)
refits$p_chi2_8 <- 1 - pchisq(refits$LR, df = 8)
write.csv(refits, file.path("output", "followup_null_refits.csv"), row.names = FALSE)
cat("\nWorst null replicates refit with BFGS polish:\n")
print(refits, digits = 3)
cat("Expected share of nulls with LR >", round(min(refits$LR), 1), "under chi2_8:",
    signif(1 - pchisq(min(refits$LR), 8), 2), "; observed:", nrow(refits), "of", nrow(null_reps), "\n")

# 2b. Null-calibrated threshold --------------------------------------------------------
# With beta = 0 the growth and level truths are the same process, so the null
# replicates of the main and level_truth studies are 200 draws from one null
# distribution. A 5th percentile from 100 draws alone is too noisy (-7.45 vs -3.21).
level <- readRDS(file.path("output", "sim_results_level_truth.rds"))$summary
null_level <- unique(level[level$scenario == "null", c("rep", "dAICc_era")])
null_pooled <- c(null_reps$dAICc_era, null_level$dAICc_era)
threshold <- unname(quantile(null_pooled, 0.05))
set.seed(1)
boot <- replicate(2000, quantile(sample(null_pooled, replace = TRUE), 0.05))
threshold_ci <- unname(quantile(boot, c(0.05, 0.95)))

power <- aggregate(dAICc_era ~ scenario + b, reps_dAICc,
                   function(x) c(naive = mean(x < -2), calibrated = mean(x < threshold)))
power <- do.call(data.frame, power)
names(power) <- c("scenario", "b", "flag_rate_dAICc_lt_minus2", "flag_rate_calibrated")
power$calibrated_threshold <- threshold
power$threshold_lo90 <- threshold_ci[1]
power$threshold_hi90 <- threshold_ci[2]
write.csv(power, file.path("output", "followup_era_calibrated_power.csv"), row.names = FALSE)
cat("\nNaive (dAICc < -2) false-positive rate, pooled nulls:", mean(null_pooled < -2), "\n")
cat("Null-calibrated threshold (5% false positives, 200 pooled nulls): dAICc <", round(threshold, 2),
    " bootstrap 90% interval", round(threshold_ci, 2), "\n")
print(power[, 1:4], digits = 2)

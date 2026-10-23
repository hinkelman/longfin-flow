# Functions for the simulation study and real-data fits: generate data with a known
# (possibly time-varying) flow effect, fit the paper's MARSS and DLM models and the
# alternatives to it, and score the fitted flow effects against the truth.

suppressMessages(library(MARSS))
source("R/design.R")

# Simulate station-level log CPUE ------------------------------------------------
# truth = "growth" (the paper's MARSS model; flow changes the growth of log CPUE):
#   x_t = x_{t-1} + beta_t * flow_t + w_t,  w_t ~ N(0, Q)
# truth = "level" (the paper's DLM model; flow shifts the level of log CPUE):
#   x_t = a_t + beta_t * flow_t,  a_t = a_{t-1} + w_t
# Observation: y_it = x_{state(i),t} + v_it,  v_it ~ N(0, R_gear(i))
simulate_data <- function(setup, beta, truth = "growth") {
  states <- names(setup$x0)
  TT <- length(setup$years)
  x <- matrix(NA_real_, length(states), TT, dimnames = list(states, setup$years))
  prev <- setup$x0
  for (t in seq_len(TT)) {
    w <- rnorm(length(states), 0, sqrt(setup$Q))
    if (truth == "growth") {
      prev <- prev + beta[t] * setup$flow_z[t] + w
      x[, t] <- prev
    } else {
      prev <- prev + w
      x[, t] <- prev + beta[t] * setup$flow_z[t]
    }
  }
  series_state <- as.character(setup$Z)
  obs_sd <- sqrt(setup$R[as.character(setup$gear)])
  y <- x[series_state, , drop = FALSE] +
    matrix(rnorm(length(series_state) * TT, 0, obs_sd), ncol = TT)
  y[!setup$mask] <- NA
  dimnames(y) <- dimnames(setup$y_obs)
  list(x = x, y = y)
}

# Fitting helpers -------------------------------------------------------------------
quiet_marss <- function(...) {
  fit <- NULL
  invisible(capture.output(fit <- suppressWarnings(MARSS(..., silent = TRUE))))
  fit
}

# Stage 1 of the paper: Z6 MARSS without covariates, used to produce DLM inputs
fit_stage1 <- function(y, setup) {
  quiet_marss(y, model = list(Z = setup$Z, Q = "diagonal and unequal", R = setup$R_model,
                              B = "identity", A = "zero"),
              control = list(maxit = 10000))
}

# Fit a flow model warm-started at the stage-1 variance estimates (a variance of
# exactly 0 is not a valid start, so floor it), falling back to a cold start
fit_from_stage1 <- function(fit_fun, y, setup, s1, form = "growth") {
  inits <- list(Q = pmax(coef(s1)$Q, 1e-4), R = coef(s1)$R, x0 = coef(s1)$x0)
  fit <- fit_fun(y, setup, inits = inits, form = form)
  if (is.null(fit$par)) fit <- fit_fun(y, setup, form = form)
  fit
}

# Refit with BFGS starting from a fitted model, keeping whichever fit is better
polish_bfgs <- function(fit, y) {
  b <- quiet_marss(y, model = fit$call$model, inits = fit, method = "BFGS",
                   control = list(maxit = 5000))
  if (!is.null(b$par) && b$logLik > fit$logLik) b else fit
}

# Static and era flow models ----------------------------------------------------------
# Flow can enter in two places (form):
#   "growth": in the process equation, x_t = x_{t-1} + C c_t + w_t (the paper's MARSS)
#   "level":  in the observation equation, y_t = Z x_t + D d_t + v_t, so flow shifts
#             the observed level of log CPUE (the structure the paper's DLM assumes)
# Each habitat state gets its own coefficient per covariate row, named "<state>" for
# one covariate and "<state> <covariate>" for several.
flow_coef_matrix <- function(setup, covariates, form) {
  states <- names(setup$x0)
  # in the level form D has one row per observed series, shared within a state
  rows <- if (form == "growth") states else as.character(setup$Z)
  cols <- rownames(covariates)
  M <- matrix(list(0), length(rows), length(cols))
  for (j in seq_along(cols)) {
    M[, j] <- if (length(cols) == 1) rows else paste(rows, cols[j])
  }
  M
}

fit_flow_model <- function(y, setup, covariates, form = "growth", inits = NULL) {
  model <- list(Z = setup$Z, Q = "diagonal and unequal", R = setup$R_model,
                B = "identity", A = "zero")
  if (form == "growth") {
    model$c <- covariates
    model$C <- flow_coef_matrix(setup, covariates, form)
  } else {
    model$d <- covariates
    model$D <- flow_coef_matrix(setup, covariates, form)
  }
  quiet_marss(y, model = model, inits = inits, control = list(maxit = 10000))
}

flow_covariate <- function(setup) matrix(setup$flow_z, nrow = 1, dimnames = list("flow", NULL))

# Era covariates: flow in one era, 0 in the others (one row per era)
era_covariates <- function(setup) {
  era <- era_of(setup$years)
  cc <- t(sapply(levels(era), function(e) setup$flow_z * (era == e)))
  dimnames(cc) <- list(levels(era), setup$years)
  cc
}

# One flow effect per state (the paper's MARSS in the growth form)
fit_static_flow <- function(y, setup, inits = NULL, form = "growth") {
  fit_flow_model(y, setup, flow_covariate(setup), form, inits)
}

# Separate flow effects per state for each era
fit_era_flow <- function(y, setup, inits = NULL, form = "growth") {
  fit_flow_model(y, setup, era_covariates(setup), form, inits)
}

# Flow effects of a fitted static or era model: one row per state and covariate,
# with Hessian-based SEs (NA if the Hessian fails)
flow_coef_table <- function(fit, setup, form = "growth") {
  covariates <- if (form == "growth") fit$call$model$c else fit$call$model$d
  M <- flow_coef_matrix(setup, covariates, form)
  # first row of each state (in the level form, rows are series within states)
  rows <- if (form == "growth") seq_len(nrow(M)) else match(names(setup$x0), as.character(setup$Z))
  par_names <- unlist(M[rows, , drop = FALSE])
  prefix <- if (form == "growth") "C" else "D"
  td <- tryCatch(MARSS::tidy(fit), error = function(e) NULL)
  data.frame(state = rep(names(setup$x0), ncol(M)),
             covariate = rep(rownames(covariates), each = length(rows)),
             estimate = coef(fit)[[prefix]][par_names, 1],
             se = if (is.null(td)) NA_real_ else
               td$std.error[match(paste0(prefix, ".", par_names), td$term)],
             row.names = NULL)
}

# Time-varying flow models --------------------------------------------------------------
# DLM with random-walk coefficients, fitted by BFGS as in the paper:
#   y_t = alpha_t + beta_t * flow_t + v_t  (intercept = TRUE, the paper's 3_DLM_multivariate_v4.R)
#   y_t = beta_t * flow_t + v_t            (intercept = FALSE)
# y has one row per habitat state; the slopes are the last nrow(y) DLM states.
# By default x0 is estimated (as in the paper). With prior_var, x0 is instead a
# prior with mean x0_mean (default 0) and variance prior_var, so the reported state
# SEs include the uncertainty about the starting values (see fit_tv_marss).
fit_dlm <- function(y, flow, intercept = TRUE, prior_var = NULL, x0_mean = NULL) {
  k <- nrow(y)
  TT <- ncol(y)
  m <- if (intercept) 2 * k else k
  q_names <- paste0("q_slo_", rownames(y))
  if (intercept) q_names <- c(paste0("q_int_", rownames(y)), q_names)
  Q <- matrix(list(0), m, m)
  diag(Q) <- q_names
  Z <- array(0, c(k, m, TT))
  for (i in seq_len(k)) {
    if (intercept) Z[i, i, ] <- 1
    Z[i, m - k + i, ] <- flow
  }
  mod <- list(B = diag(m), U = matrix(0, m, 1), Q = Q, Z = Z,
              A = matrix(0, k, 1), R = "diagonal and unequal")
  inits <- list(x0 = matrix(0, nrow = m))
  if (!is.null(prior_var)) {
    mod$x0 <- matrix(if (is.null(x0_mean)) 0 else x0_mean, m, 1)
    mod$V0 <- diag(prior_var, m)
    mod$tinitx <- 0
    inits <- NULL
  }
  quiet_marss(y, inits = inits, model = mod, method = "BFGS", control = list(maxit = 10000))
}

# One-stage time-varying flow effect, fitted to the station data (not in the paper).
# The state vector is (x, beta), with beta_t = beta_{t-1} + e_t, and
#   form = "growth": x_t = x_{t-1} + beta_{t-1} * flow_t + w_t   (B_t carries flow_t)
#   form = "level":  x_t = x_{t-1} + w_t,  y_t = x_t + beta_t * flow_t + v_t   (Z_t carries flow_t)
# EM alone crawls here (10,000 iterations without converging), so run 100 EM
# iterations and finish with BFGS, which reaches a higher likelihood in seconds.
# By default x0 is estimated, and MARSS then reports state SEs as if x0 were
# known: when a slope's variance goes to 0 its SE collapses toward 0. With
# prior_var, x0 is a prior with mean x0_mean and variance prior_var instead.
fit_tv_marss <- function(y, setup, inits = NULL, form = "growth", prior_var = NULL,
                         x0_mean = NULL) {
  states <- names(setup$x0)
  k <- length(states)
  TT <- ncol(y)
  series_state <- match(as.character(setup$Z), states)
  Q <- matrix(list(0), 2 * k, 2 * k)
  diag(Q) <- c(paste0("q_", states), paste0("q_slo_", states))
  if (form == "growth") {
    B <- array(0, c(2 * k, 2 * k, TT))
    for (t in seq_len(TT)) {
      B[, , t] <- diag(2 * k)
      B[cbind(1:k, k + 1:k, t)] <- setup$flow_z[t]
    }
    Z <- matrix(0, nrow(y), 2 * k)
    Z[cbind(seq_len(nrow(y)), series_state)] <- 1
  } else {
    B <- diag(2 * k)
    Z <- array(0, c(nrow(y), 2 * k, TT))
    for (t in seq_len(TT)) {
      Z[cbind(seq_len(nrow(y)), series_state, t)] <- 1
      Z[cbind(seq_len(nrow(y)), k + series_state, t)] <- setup$flow_z[t]
    }
  }
  model <- list(B = B, U = "zero", Q = Q, Z = Z, A = "zero", R = setup$R_model,
                x0 = "unequal", tinitx = 0)
  if (!is.null(prior_var)) {
    # x0 is fixed (the prior mean), so it can't also have starting values
    model$x0 <- matrix(x0_mean, ncol = 1)
    model$V0 <- diag(prior_var, 2 * k)
    inits$x0 <- NULL
  }
  em <- quiet_marss(y, model = model, inits = inits, control = list(maxit = 100))
  if (is.null(em$par)) return(em)
  # BFGS cannot start a variance at exactly 0, so floor the EM variances
  start <- list(Q = pmax(coef(em)$Q, 1e-6), R = coef(em)$R)
  if (is.null(prior_var)) start$x0 <- coef(em)$x0
  bf <- quiet_marss(y, model = model, inits = start, method = "BFGS", control = list(maxit = 5000))
  if (!is.null(bf$par) && bf$logLik >= em$logLik) bf else em
}

# Slope paths (states x years) and SEs ----------------------------------------------------
# From fit_dlm: the slopes are the last length(states) DLM states, one column per year
dlm_slopes <- function(fit, states, years) {
  idx <- nrow(fit$states) - length(states) + seq_along(states)
  sl <- fit$states[idx, , drop = FALSE]
  se <- fit$states.se[idx, , drop = FALSE]
  dimnames(sl) <- dimnames(se) <- list(states, years)
  list(slope = sl, se = se)
}

# From fit_tv_marss, aligned to setup$years. In the growth form beta_{t-1}
# multiplies flow_t, so year t uses state column t-1 and the first year has no
# estimate (NA).
tv_slopes <- function(fit, setup, form = "growth") {
  k <- length(setup$x0)
  idx <- k + seq_len(k)
  sl <- fit$states[idx, , drop = FALSE]
  se <- fit$states.se[idx, , drop = FALSE]
  if (form == "growth") {
    sl <- cbind(NA, sl[, -ncol(sl), drop = FALSE])
    se <- cbind(NA, se[, -ncol(se), drop = FALSE])
  }
  dimnames(sl) <- dimnames(se) <- list(names(setup$x0), setup$years)
  list(slope = sl, se = se)
}

# Scoring against the truth ----------------------------------------------------------------
# Slope summaries against the true effect, for a matrix of slope paths (states x years)
slope_metrics <- function(slope, slope_se, beta, q_slope) {
  true <- matrix(beta, nrow(slope), length(beta), byrow = TRUE)
  data.frame(state = rownames(slope),
             mean_slope = rowMeans(slope),
             first = slope[, 1],
             last = slope[, ncol(slope)],
             q_slope = unname(q_slope),
             mean_se = rowMeans(slope_se),
             cover = rowMeans(abs(slope - true) <= 1.96 * slope_se),
             rmse = sqrt(rowMeans((slope - true)^2)),
             true_mean = mean(beta), true_first = beta[1], true_last = beta[length(beta)])
}

# Score a named list of slope fits on the years in score_years.
#   slopes: function(name, fit) returning list(slope, se) with states x years dimnames
# Failed fits are dropped. Returns per-state metrics and the slope paths.
score_slope_fits <- function(fits, slopes, beta, years, score_years) {
  out <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    if (is.null(f$par)) return(NULL)
    p <- slopes(v, f)
    cols <- colnames(p$slope) %in% score_years
    sl <- p$slope[, cols, drop = FALSE]
    se <- p$se[, cols, drop = FALSE]
    states <- rownames(sl)
    true <- beta[years %in% score_years]
    m <- slope_metrics(sl, se, true, coef(f)$Q[paste0("q_slo_", states), 1])
    m$variant <- v
    m$conv <- f$convergence
    list(metrics = m,
         paths = data.frame(variant = v, state = rep(states, times = ncol(sl)),
                            year = rep(as.numeric(colnames(sl)), each = length(states)),
                            true_beta = rep(true, each = length(states)),
                            slope = as.vector(sl), se = as.vector(se)))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  list(metrics = do.call(rbind, lapply(out, `[[`, "metrics")),
       paths = do.call(rbind, lapply(out, `[[`, "paths")))
}

# Studies ---------------------------------------------------------------------------------
# The DLM variants and interval checks are scored on 1981-2020, the years with a
# growth increment. (The main study scores the paper's DLM on all of 1980-2020; it
# predates the variants and is kept as run.)
score_years <- function(setup) setup$years[-1]

# Decomposing the paper's DLM: change one modeling choice at a time.
#   A_paper          level ~ flow, stage-1 smoothed states, random-walk intercept (paper)
#   B_true_states    level ~ flow, true states, random-walk intercept
#   C_growth         change in true states ~ flow, random-walk intercept
#   D_growth_no_int  change in true states ~ flow, no intercept
#   E_one_stage      time-varying flow effect in the MARSS process equation, station data
# B-D use the simulated true states, which are never observed; they are diagnostics.
fit_dlm_variants <- function(setup, sim, s1, beta) {
  states <- names(setup$x0)
  years <- setup$years
  flow <- setup$flow_z
  st1 <- s1$states
  rownames(st1) <- states
  growth <- sim$x[, -1, drop = FALSE] - sim$x[, -ncol(sim$x), drop = FALSE]
  x0_mean <- c(coef(s1)$x0[, 1], rep(0, length(states)))

  fits <- list(
    A_paper = fit_dlm(st1, flow),
    B_true_states = fit_dlm(sim$x, flow),
    C_growth = fit_dlm(growth, flow[-1]),
    D_growth_no_int = fit_dlm(growth, flow[-1], intercept = FALSE),
    E_one_stage = fit_tv_marss(sim$y, setup, inits = list(x0 = matrix(x0_mean, ncol = 1)))
  )
  slopes <- function(v, f) {
    switch(v,
           A_paper = , B_true_states = dlm_slopes(f, states, years),
           C_growth = , D_growth_no_int = dlm_slopes(f, states, years[-1]),
           E_one_stage = tv_slopes(f, setup, "growth"))
  }
  score_slope_fits(fits, slopes, beta, years, score_years(setup))
}

# Structure and interval checks:
#  structure: static flow effect as growth (C) or level (D), and the era versions of
#             both; AICc says which structure the data favor
#  intervals: the one-stage time-varying model in the form matching the truth, with
#             x0 estimated or given a prior, and the paper's DLM with a prior
fit_checks <- function(setup, sim, s1, beta, truth) {
  states <- names(setup$x0)
  k <- length(states)
  fc <- fit_from_stage1(fit_static_flow, sim$y, setup, s1, "growth")
  fd <- fit_from_stage1(fit_static_flow, sim$y, setup, s1, "level")
  ec <- fit_from_stage1(fit_era_flow, sim$y, setup, s1, "growth")
  ed <- fit_from_stage1(fit_era_flow, sim$y, setup, s1, "level")
  C <- flow_coef_table(fc, setup, "growth")
  D <- flow_coef_table(fd, setup, "level")
  structure <- data.frame(state = states,
                          C_est = C$estimate, C_se = C$se, D_est = D$estimate, D_se = D$se,
                          AICc_stage1 = s1$AICc, AICc_C = fc$AICc, AICc_D = fd$AICc,
                          AICc_era_C = ec$AICc, AICc_era_D = ed$AICc,
                          conv_C = fc$convergence, conv_D = fd$convergence,
                          conv_era_C = ec$convergence, conv_era_D = ed$convergence)

  # x0 prior: abundance states centred on the stage-1 estimates, slopes on 0;
  # variance 1 is wide on both scales (log CPUE units; flow effect per SD of flow)
  x0_mean <- c(coef(s1)$x0[, 1], rep(0, k))
  st1 <- s1$states
  rownames(st1) <- states
  fits <- list(
    F_tv_matched = fit_tv_marss(sim$y, setup, inits = list(x0 = matrix(x0_mean, ncol = 1)),
                                form = truth),
    G_tv_matched_prior = fit_tv_marss(sim$y, setup, form = truth, prior_var = 1, x0_mean = x0_mean),
    H_paper_dlm_prior = fit_dlm(st1, setup$flow_z, prior_var = 1, x0_mean = c(st1[, 1], rep(0, k)))
  )
  slopes <- function(v, f) {
    if (v == "H_paper_dlm_prior") dlm_slopes(f, states, setup$years) else tv_slopes(f, setup, truth)
  }
  scored <- score_slope_fits(fits, slopes, beta, setup$years, score_years(setup))
  c(list(structure = structure), scored)
}

# The paper's three models on one simulated data set: static MARSS C, era MARSS
# (growth form), and the DLM on the stage-1 states
fit_main <- function(setup, sim, s1, beta) {
  states <- names(setup$x0)
  years <- setup$years
  mc <- fit_from_stage1(fit_static_flow, sim$y, setup, s1)
  me <- fit_from_stage1(fit_era_flow, sim$y, setup, s1)
  dlm_states <- s1$states
  rownames(dlm_states) <- states
  dl <- fit_dlm(dlm_states, setup$flow_z)

  C <- flow_coef_table(mc, setup)
  marss <- data.frame(state = states, marss_C = C$estimate, marss_se = C$se)

  era_C <- flow_coef_table(me, setup)
  era <- data.frame(state = states,
                    era_pre1987 = era_C$estimate[era_C$covariate == "pre1987"],
                    era_1987_2002 = era_C$estimate[era_C$covariate == "1987_2002"],
                    era_post2002 = era_C$estimate[era_C$covariate == "post2002"],
                    dAICc_era = me$AICc - mc$AICc)

  p <- dlm_slopes(dl, states, years)
  q <- coef(dl)$Q[, 1]
  m <- slope_metrics(p$slope, p$se, beta, q[paste0("q_slo_", states)])
  dlm <- data.frame(state = states,
                    dlm_mean = m$mean_slope, dlm_first = m$first, dlm_last = m$last,
                    dlm_q_slope = m$q_slope, dlm_q_int = unname(q[paste0("q_int_", states)]),
                    dlm_mean_se = m$mean_se, dlm_cover = m$cover, dlm_rmse = m$rmse)

  era_of_year <- era_of(years)
  summary <- Reduce(function(a, b) merge(a, b, by = "state"), list(marss, era, dlm))
  summary <- cbind(summary,
                   true_mean = mean(beta), true_first = beta[1], true_last = beta[length(beta)],
                   true_pre1987 = mean(beta[era_of_year == "pre1987"]),
                   true_1987_2002 = mean(beta[era_of_year == "1987_2002"]),
                   true_post2002 = mean(beta[era_of_year == "post2002"]),
                   conv_s1 = s1$convergence, conv_c = mc$convergence,
                   conv_era = me$convergence, conv_dlm = dl$convergence)
  paths <- data.frame(state = rep(states, times = length(years)),
                      year = rep(years, each = length(states)),
                      true_beta = rep(beta, each = length(states)),
                      dlm_slope = as.vector(p$slope), dlm_se = as.vector(p$se))
  list(summary = summary, paths = paths)
}

# One replicate: simulate, fit, return tidy summaries ----------------------------------------
#   truth:    "growth" (paper's MARSS model) or "level" (paper's DLM model)
#   main:     fit the paper's MARSS C, the paper's DLM, and the era model (fit_main)
#   variants: also fit the DLM variants (fit_dlm_variants)
#   checks:   also fit the structure and interval checks (fit_checks)
run_replicate <- function(setup, scenario, b, rep_id, truth = "growth", main = TRUE,
                          variants = FALSE, checks = FALSE) {
  beta <- beta_path(scenario, setup$years, b)
  sim <- simulate_data(setup, beta, truth)
  ids <- data.frame(scenario = scenario, b = b, rep = rep_id, truth = truth)
  s1 <- fit_stage1(sim$y, setup)
  result <- list()

  if (main) {
    f <- fit_main(setup, sim, s1, beta)
    summary <- cbind(ids[c("scenario", "b", "rep")], f$summary)
    # the main study's results have no truth column (all growth); level_truth's do
    if (truth != "growth") summary$truth <- truth
    result$summary <- summary
    result$paths <- cbind(ids[c("scenario", "b", "rep")], f$paths)
  }
  if (variants) {
    v <- fit_dlm_variants(setup, sim, s1, beta)
    result$variants <- cbind(ids, v$metrics)
    result$variant_paths <- cbind(ids, v$paths)
  }
  if (checks) {
    e <- fit_checks(setup, sim, s1, beta, truth)
    result$structure <- cbind(ids, e$structure)
    result$variants <- rbind(result$variants, cbind(ids, e$metrics))
    result$variant_paths <- rbind(result$variant_paths, cbind(ids, e$paths))
  }
  result
}

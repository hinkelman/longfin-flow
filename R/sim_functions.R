# Functions for the simulation study: generate data with a known (possibly
# time-varying) flow effect, then fit the paper's MARSS and DLM models to it.

suppressMessages(library(MARSS))

# True flow-effect trajectories -------------------------------------------------
# b is the flow effect at the start of the series (log-CPUE growth per SD of flow)
beta_path <- function(scenario, years, b) {
  switch(scenario,
         null      = rep(0, length(years)),
         constant  = rep(b, length(years)),
         decline   = b * (max(years) - years) / (max(years) - min(years)),
         step_1987 = ifelse(years < 1987, b, 0),
         step_2002 = ifelse(years < 2002, b, 0),
         stop("unknown scenario"))
}

# Simulate station-level log CPUE ------------------------------------------------
# State model: x_t = x_{t-1} + beta_t * flow_t + w_t,  w_t ~ N(0, Q)
# Observation: y_it = x_{state(i),t} + v_it,          v_it ~ N(0, R_gear(i))
simulate_data <- function(setup, beta) {
  states <- names(setup$x0)
  TT <- length(setup$years)
  x <- matrix(NA_real_, length(states), TT, dimnames = list(states, setup$years))
  prev <- setup$x0
  for (t in seq_len(TT)) {
    prev <- prev + beta[t] * setup$flow_z[t] + rnorm(length(states), 0, sqrt(setup$Q))
    x[, t] <- prev
  }
  series_state <- as.character(setup$Z)
  obs_sd <- sqrt(setup$R[as.character(setup$gear)])
  y <- x[series_state, , drop = FALSE] +
    matrix(rnorm(length(series_state) * TT, 0, obs_sd), ncol = TT)
  y[!setup$mask] <- NA
  dimnames(y) <- dimnames(setup$y_obs)
  list(x = x, y = y)
}

# Model fits ----------------------------------------------------------------------
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

# Paper's time-averaged model: one flow coefficient (C) per state
fit_marss_c <- function(y, setup, inits = NULL) {
  quiet_marss(y, model = list(Z = setup$Z, Q = "diagonal and unequal", R = setup$R_model,
                              B = "identity", A = "zero",
                              c = matrix(setup$flow_z, nrow = 1), C = "unequal"),
              inits = inits, control = list(maxit = 10000))
}

# Era model: separate flow coefficient per state for 1980-86, 1987-2002, 2003-20
era_covariates <- function(setup) {
  era <- cut(setup$years, c(-Inf, 1986, 2002, Inf), labels = c("pre1987", "1987_2002", "post2002"))
  cc <- t(sapply(levels(era), function(e) setup$flow_z * (era == e)))
  dimnames(cc) <- list(levels(era), setup$years)
  cc
}

fit_marss_era <- function(y, setup, inits = NULL) {
  quiet_marss(y, model = list(Z = setup$Z, Q = "diagonal and unequal", R = setup$R_model,
                              B = "identity", A = "zero",
                              c = era_covariates(setup), C = "unconstrained"),
              inits = inits, control = list(maxit = 10000))
}

# Paper's DLM (3_DLM_multivariate_v4.R): states_t = alpha_t + beta_t * flow_t + v_t,
# random-walk alpha and beta for each of the 4 groups, fitted by BFGS
fit_dlm <- function(states, setup) {
  m <- 8
  TT <- ncol(states)
  Q <- matrix(list(0), m, m)
  diag(Q) <- c(paste0("q_int_", rownames(states)), paste0("q_slo_", rownames(states)))
  Z <- array(0, c(4, m, TT))
  for (i in 1:4) {
    Z[i, i, ] <- 1
    Z[i, i + 4, ] <- setup$flow_z
  }
  mod <- list(B = diag(m), U = matrix(0, m, 1), Q = Q, Z = Z,
              A = matrix(0, 4, 1), R = "diagonal and unequal")
  quiet_marss(states, inits = list(x0 = matrix(0, nrow = m)), model = mod,
              method = "BFGS", control = list(maxit = 10000))
}

# Coefficient table (estimate, SE) from a MARSS fit via the Hessian
coef_table <- function(fit, prefix = "C.") {
  td <- tryCatch(MARSS::tidy(fit), error = function(e) NULL)
  if (is.null(td)) return(NULL)
  td <- td[startsWith(td$term, prefix), c("term", "estimate", "std.error")]
  rownames(td) <- NULL
  td
}

# One replicate: simulate, fit everything, return tidy summaries -----------------
run_replicate <- function(setup, scenario, b, rep_id) {
  beta <- beta_path(scenario, setup$years, b)
  sim <- simulate_data(setup, beta)
  states <- names(setup$x0)

  s1 <- fit_stage1(sim$y, setup)
  # warm-start the covariate models at stage-1 variance estimates (a variance
  # of exactly 0 is not a valid start, so floor it), falling back to a cold start
  s1_inits <- list(Q = pmax(coef(s1)$Q, 1e-4), R = coef(s1)$R, x0 = coef(s1)$x0)
  mc <- fit_marss_c(sim$y, setup, inits = s1_inits)
  if (is.null(mc$par)) mc <- fit_marss_c(sim$y, setup)
  me <- fit_marss_era(sim$y, setup, inits = s1_inits)
  if (is.null(me$par)) me <- fit_marss_era(sim$y, setup)
  dlm_states <- s1$states
  rownames(dlm_states) <- states
  dl <- fit_dlm(dlm_states, setup)

  ct <- coef_table(mc)
  marss <- data.frame(state = states,
                      marss_C = coef(mc)$C[paste0(states), 1] %||% NA,
                      marss_se = NA_real_)
  if (!is.null(ct)) marss$marss_se <- ct$std.error[match(paste0("C.", states), ct$term)]

  era_C <- coef(me, type = "matrix")$C
  era <- data.frame(state = states,
                    era_pre1987 = era_C[states, 1],
                    era_1987_2002 = era_C[states, 2],
                    era_post2002 = era_C[states, 3],
                    dAICc_era = me$AICc - mc$AICc)

  q <- coef(dl)$Q[, 1]
  slope <- dl$states[5:8, , drop = FALSE]
  slope_se <- dl$states.se[5:8, , drop = FALSE]
  dlm <- data.frame(state = states,
                    dlm_mean = rowMeans(slope),
                    dlm_first = slope[, 1],
                    dlm_last = slope[, ncol(slope)],
                    dlm_q_slope = q[paste0("q_slo_", states)],
                    dlm_q_int = q[paste0("q_int_", states)],
                    dlm_mean_se = rowMeans(slope_se),
                    dlm_cover = rowMeans(abs(slope - matrix(beta, 4, length(beta), byrow = TRUE)) <=
                                           1.96 * slope_se),
                    dlm_rmse = sqrt(rowMeans((slope - matrix(beta, 4, length(beta), byrow = TRUE))^2)))

  summary <- Reduce(function(a, b) merge(a, b, by = "state"), list(marss, era, dlm))
  summary <- cbind(scenario = scenario, b = b, rep = rep_id, summary,
                   true_mean = mean(beta), true_first = beta[1], true_last = beta[length(beta)],
                   true_pre1987 = mean(beta[setup$years < 1987]),
                   true_1987_2002 = mean(beta[setup$years >= 1987 & setup$years <= 2002]),
                   true_post2002 = mean(beta[setup$years > 2002]),
                   conv_s1 = s1$convergence, conv_c = mc$convergence,
                   conv_era = me$convergence, conv_dlm = dl$convergence)

  paths <- data.frame(scenario = scenario, b = b, rep = rep_id,
                      state = rep(states, times = length(setup$years)),
                      year = rep(setup$years, each = 4),
                      true_beta = rep(beta, each = 4),
                      dlm_slope = as.vector(slope),
                      dlm_se = as.vector(slope_se))
  list(summary = summary, paths = paths)
}

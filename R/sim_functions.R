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

# Flow can enter in two places (form):
#   "growth": in the process equation, x_t = x_{t-1} + C c_t + w_t (the paper's MARSS)
#   "level":  in the observation equation, y_t = Z x_t + D d_t + v_t, so flow shifts
#             the observed level of log CPUE (the structure the paper's DLM assumes)
# Each habitat state gets its own coefficient per covariate row, named "<state>" for
# one covariate and "<state> <covariate>" for several.
flow_coef_matrix <- function(setup, covariates, form) {
  states <- names(setup$x0)
  if (form == "growth") {
    rows <- states
  } else {
    rows <- as.character(setup$Z)  # one row per observed series, shared within a state
  }
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

# Paper's time-averaged model: one flow coefficient per state
fit_marss_c <- function(y, setup, inits = NULL, form = "growth") {
  fit_flow_model(y, setup, matrix(setup$flow_z, nrow = 1, dimnames = list("flow", NULL)),
                 form, inits)
}

# Era model: separate flow coefficient per state for 1980-86, 1987-2002, 2003-20
era_covariates <- function(setup) {
  era <- cut(setup$years, c(-Inf, 1986, 2002, Inf), labels = c("pre1987", "1987_2002", "post2002"))
  cc <- t(sapply(levels(era), function(e) setup$flow_z * (era == e)))
  dimnames(cc) <- list(levels(era), setup$years)
  cc
}

fit_marss_era <- function(y, setup, inits = NULL, form = "growth") {
  fit_flow_model(y, setup, era_covariates(setup), form, inits)
}

# Flow coefficients of a fitted flow model as a states x covariates matrix
flow_coefs <- function(fit, setup, form = "growth") {
  M <- coef(fit, type = "matrix")[[if (form == "growth") "C" else "D"]]
  # in the level form D has one row per series; take the first series of each state
  M <- M[match(names(setup$x0), if (form == "growth") rownames(M) else as.character(setup$Z)), ,
         drop = FALSE]
  rownames(M) <- names(setup$x0)
  M
}

# DLM with random-walk coefficients, fitted by BFGS as in the paper:
#   y_t = alpha_t + beta_t * flow_t + v_t  (intercept = TRUE, the paper's 3_DLM_multivariate_v4.R)
#   y_t = beta_t * flow_t + v_t            (intercept = FALSE)
# y has one row per habitat state; the slopes are the last nrow(y) DLM states.
# By default x0 is estimated (as in the paper). With prior_var, x0 is instead a
# prior with mean x0_mean and variance prior_var, so the reported state SEs include the
# uncertainty about the starting values (see fit_tv_marss).
fit_dlm <- function(y, flow, intercept = TRUE, prior_var = NULL, x0_mean = rep(0, m)) {
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
    mod$x0 <- matrix(x0_mean, m, 1)
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

# Slope paths (states x years) and SEs from fit_tv_marss, aligned to setup$years.
# In the growth form beta_{t-1} multiplies flow_t, so year t uses state column t-1
# and the first year has no estimate (NA).
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

# Coefficient table (estimate, SE) from a MARSS fit via the Hessian
coef_table <- function(fit, prefix = "C.") {
  td <- tryCatch(MARSS::tidy(fit), error = function(e) NULL)
  if (is.null(td)) return(NULL)
  td <- td[startsWith(td$term, prefix), c("term", "estimate", "std.error")]
  rownames(td) <- NULL
  td
}

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

# Decomposing the paper's DLM: change one modeling choice at a time.
#   A_paper          level ~ flow, stage-1 smoothed states, random-walk intercept (paper)
#   B_true_states    level ~ flow, true states, random-walk intercept
#   C_growth         change in true states ~ flow, random-walk intercept
#   D_growth_no_int  change in true states ~ flow, no intercept
#   E_one_stage      time-varying flow effect in the MARSS process equation, station data
# B-D use the simulated true states, which are never observed; they are diagnostics.
# All variants are scored on 1981-2020, the years with a growth increment.
fit_dlm_variants <- function(setup, sim, s1, beta) {
  states <- names(setup$x0)
  k <- length(states)
  keep <- -1
  flow <- setup$flow_z
  st1 <- s1$states
  rownames(st1) <- states
  growth <- sim$x[, -1, drop = FALSE] - sim$x[, -ncol(sim$x), drop = FALSE]

  fits <- list(
    A_paper = fit_dlm(st1, flow),
    B_true_states = fit_dlm(sim$x, flow),
    C_growth = fit_dlm(growth, flow[-1]),
    D_growth_no_int = fit_dlm(growth, flow[-1], intercept = FALSE),
    E_one_stage = fit_tv_marss(sim$y, setup,
                               inits = list(x0 = matrix(c(coef(s1)$x0[, 1], rep(0, k)), ncol = 1)))
  )

  out <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    if (is.null(f$par)) return(NULL)
    q <- coef(f)$Q[, 1]
    q_slope <- q[paste0("q_slo_", states)]
    sl <- f$states[nrow(f$states) - k + 1:k, , drop = FALSE]
    se <- f$states.se[nrow(f$states) - k + 1:k, , drop = FALSE]
    if (v %in% c("A_paper", "B_true_states")) {
      sl <- sl[, keep, drop = FALSE]
      se <- se[, keep, drop = FALSE]
    }
    if (v == "E_one_stage") {
      # beta_{t-1} multiplies flow_t, so the effect in year t is state column t-1
      sl <- sl[, -ncol(sl), drop = FALSE]
      se <- se[, -ncol(se), drop = FALSE]
    }
    rownames(sl) <- rownames(se) <- states
    m <- slope_metrics(sl, se, beta[keep], q_slope)
    m$variant <- v
    m$conv <- f$convergence
    list(metrics = m,
         paths = data.frame(variant = v, state = rep(states, times = ncol(sl)),
                            year = rep(setup$years[keep], each = k),
                            true_beta = rep(beta[keep], each = k),
                            slope = as.vector(sl), se = as.vector(se)))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  list(metrics = do.call(rbind, lapply(out, `[[`, "metrics")),
       paths = do.call(rbind, lapply(out, `[[`, "paths")))
}

# Extra fits (Item 1 and Item 3 of the follow-up list):
#  structure: static flow effect as growth (C) or level (D), and the era versions of
#             both; AICc says which structure the data favor
#  intervals: the one-stage time-varying model in the form matching the truth, with
#             x0 estimated (as before) or given a prior, and the paper's DLM with a prior
# Interval variants are scored on 1981-2020, like fit_dlm_variants.
fit_extra <- function(setup, sim, s1, beta, truth) {
  states <- names(setup$x0)
  k <- length(states)
  s1_inits <- list(Q = pmax(coef(s1)$Q, 1e-4), R = coef(s1)$R, x0 = coef(s1)$x0)
  refit <- function(f, form) {
    fit <- f(sim$y, setup, inits = s1_inits, form = form)
    if (is.null(fit$par)) fit <- f(sim$y, setup, form = form)
    fit
  }
  fc <- refit(fit_marss_c, "growth")
  fd <- refit(fit_marss_c, "level")
  ec <- refit(fit_marss_era, "growth")
  ed <- refit(fit_marss_era, "level")
  se_of <- function(fit, prefix) {
    ct <- coef_table(fit, prefix)
    if (is.null(ct)) return(rep(NA_real_, k))
    ct$std.error[match(paste0(prefix, states), ct$term)]
  }
  structure <- data.frame(state = states,
                          C_est = flow_coefs(fc, setup, "growth")[, 1], C_se = se_of(fc, "C."),
                          D_est = flow_coefs(fd, setup, "level")[, 1], D_se = se_of(fd, "D."),
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
  keep <- -1
  out <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    if (is.null(f$par)) return(NULL)
    q <- coef(f)$Q[, 1]
    if (v == "H_paper_dlm_prior") {
      sl <- f$states[k + 1:k, keep, drop = FALSE]
      se <- f$states.se[k + 1:k, keep, drop = FALSE]
      rownames(sl) <- rownames(se) <- states
    } else {
      p <- tv_slopes(f, setup, truth)
      sl <- p$slope[, keep, drop = FALSE]
      se <- p$se[, keep, drop = FALSE]
    }
    m <- slope_metrics(sl, se, beta[keep], q[paste0("q_slo_", states)])
    m$variant <- v
    m$conv <- f$convergence
    list(metrics = m,
         paths = data.frame(variant = v, state = rep(states, times = ncol(sl)),
                            year = rep(setup$years[keep], each = k),
                            true_beta = rep(beta[keep], each = k),
                            slope = as.vector(sl), se = as.vector(se)))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  list(structure = structure,
       metrics = do.call(rbind, lapply(out, `[[`, "metrics")),
       paths = do.call(rbind, lapply(out, `[[`, "paths")))
}

# One replicate: simulate, fit everything, return tidy summaries -----------------
#   truth:    "growth" (paper's MARSS model) or "level" (paper's DLM model)
#   main:     fit the paper's MARSS C, the paper's DLM, and the era model
#   variants: also fit the DLM variants (fit_dlm_variants)
#   extra:    also fit the structure and interval checks (fit_extra)
run_replicate <- function(setup, scenario, b, rep_id, truth = "growth", main = TRUE,
                          variants = FALSE, extra = FALSE) {
  beta <- beta_path(scenario, setup$years, b)
  sim <- simulate_data(setup, beta, truth)
  states <- names(setup$x0)
  ids <- data.frame(scenario = scenario, b = b, rep = rep_id, truth = truth)

  s1 <- fit_stage1(sim$y, setup)
  result <- list()

  if (main) {
    # warm-start the covariate models at stage-1 variance estimates (a variance
    # of exactly 0 is not a valid start, so floor it), falling back to a cold start
    s1_inits <- list(Q = pmax(coef(s1)$Q, 1e-4), R = coef(s1)$R, x0 = coef(s1)$x0)
    mc <- fit_marss_c(sim$y, setup, inits = s1_inits)
    if (is.null(mc$par)) mc <- fit_marss_c(sim$y, setup)
    me <- fit_marss_era(sim$y, setup, inits = s1_inits)
    if (is.null(me$par)) me <- fit_marss_era(sim$y, setup)
    dlm_states <- s1$states
    rownames(dlm_states) <- states
    dl <- fit_dlm(dlm_states, setup$flow_z)

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
    if (truth != "growth") summary$truth <- truth

    result$summary <- summary
    result$paths <- data.frame(scenario = scenario, b = b, rep = rep_id,
                               state = rep(states, times = length(setup$years)),
                               year = rep(setup$years, each = 4),
                               true_beta = rep(beta, each = 4),
                               dlm_slope = as.vector(slope),
                               dlm_se = as.vector(slope_se))
  }

  if (variants) {
    v <- fit_dlm_variants(setup, sim, s1, beta)
    result$variants <- cbind(ids, v$metrics)
    result$variant_paths <- cbind(ids, v$paths)
  }

  if (extra) {
    e <- fit_extra(setup, sim, s1, beta, truth)
    result$structure <- cbind(ids, e$structure)
    result$variants <- rbind(result$variants, cbind(ids, e$metrics))
    result$variant_paths <- rbind(result$variant_paths, cbind(ids, e$paths))
  }
  result
}

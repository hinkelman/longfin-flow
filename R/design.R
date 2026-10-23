# Study design shared by the simulation, real-data, and summary scripts:
# years, true flow-effect scenarios, eras, the replicate grid and seeds, and the
# threshold used to call a fitted slope "time-varying".

YEARS <- 1980:2020

# True flow-effect trajectories ----------------------------------------------------
# b is the flow effect at the start of the series (per SD of flow)
SCENARIO_LABELS <- c(null = "Null (no effect)", constant = "Constant",
                     decline = "Linear decline", step_1987 = "Step at 1987",
                     step_2002 = "Step at 2002")

beta_path <- function(scenario, years, b) {
  switch(scenario,
         null      = rep(0, length(years)),
         constant  = rep(b, length(years)),
         decline   = b * (max(years) - years) / (max(years) - min(years)),
         step_1987 = ifelse(years < 1987, b, 0),
         step_2002 = ifelse(years < 2002, b, 0),
         stop("unknown scenario"))
}

# Eras for the era model and era summaries: before the 1987 clam invasion, before
# the Pelagic Organism Decline (through 2002), and after
ERA_LABELS <- c(pre1987 = "1980-86", `1987_2002` = "1987-2002", post2002 = "2003-20")

era_of <- function(years) {
  cut(years, c(-Inf, 1986, 2002, Inf), labels = names(ERA_LABELS))
}

# Replicate grid and seeds ------------------------------------------------------------
# b = 0.08 is roughly the largest time-averaged MARSS flow effect in the paper
# (age-0 channel midwater); b = 0.24 is a 3x "strong signal" check. The null does
# not depend on b, so it is run once.
sim_grid <- function(n_reps) {
  grid <- expand.grid(rep = seq_len(n_reps), scenario = names(SCENARIO_LABELS),
                      b = c(0.08, 0.24), stringsAsFactors = FALSE)
  grid <- grid[!(grid$scenario == "null" & grid$b == 0.24), ]
  rownames(grid) <- NULL
  grid
}

# Replicate-specific seeds, so results don't depend on scheduling or restarts.
# Studies with the same offset fit models to identical simulated data:
# offset 0 for growth truth, 10000 for level truth.
SEED_BASE <- 20261112
grid_seed <- function(grid_row, seed_offset = 0) SEED_BASE + seed_offset + grid_row

# Model summaries -------------------------------------------------------------------
# A fitted slope variance above this is counted as "time-varying"; below it the
# random walk has collapsed to a constant slope (fits give ~1e-11 or smaller)
Q_TV_THRESHOLD <- 1e-8

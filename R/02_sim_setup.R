# Build the inputs for the simulation study from the authors' own age-0 pipeline.
# Sources the age-0 section (lines 1-370) of the authors' MARSS script so the
# station set, missing-data pattern, Z6 grouping, and R structure are exactly
# theirs, then saves those pieces plus the fitted Z6 parameters and the
# z-scored Jan-Jun outflow covariate.

suppressMessages(library(MARSS))

root <- normalizePath(".")
run_dir <- file.path(root, "output", "original_run", "script1")
script1 <- list.files(file.path(root, "original", "code"), pattern = "^1_", full.names = TRUE)

age0_code <- readLines(script1)[1:370]
tmp <- tempfile(fileext = ".R")
writeLines(age0_code, tmp)

old_wd <- setwd(run_dir)
pdf(NULL)  # the age-0 section draws plots; discard them
env <- new.env()
source(tmp, local = env)
dev.off()
setwd(old_wd)

dayflow <- read.csv(file.path(root, "original", "data", "annual_average_deltawide_flow_jan_june.csv"))
years <- 1980:2020
flow <- dayflow$Mean_outflow_CY[match(years, dayflow$Calendar_Year)]
flow_z <- as.numeric(zscore(flow))

fit <- env$z6
sim_setup <- list(
  years = years,
  flow_z = flow_z,
  y_obs = env$MWT_OT_matrix2,          # observed age-0 log CPUE (series x years)
  mask = !is.na(env$MWT_OT_matrix2),   # real missing-data pattern
  Z = env$Z6,                          # series -> state grouping (factor)
  R_model = env$R.model,               # gear-specific observation variances
  gear = env$z_matrix$Gear,
  Q = setNames(coef(fit)$Q[, 1], rownames(fit$states)),
  R = coef(fit)$R[, 1],
  x0 = setNames(coef(fit)$x0[, 1], rownames(fit$states)),
  fit = fit
)
saveRDS(sim_setup, file.path(root, "output", "sim_setup.rds"))

cat("Series:", nrow(sim_setup$y_obs), " Years:", ncol(sim_setup$y_obs),
    " Fraction observed:", round(mean(sim_setup$mask), 3), "\n")
cat("Series per state:\n"); print(table(sim_setup$Z))
cat("Q:\n"); print(sim_setup$Q)
cat("R:\n"); print(sim_setup$R)
cat("AICc of refit age-0 Z6 (paper: 803.29):", round(fit$AICc, 2), "\n")
cat("Flow lag-1 autocorrelation:", round(acf(flow_z, plot = FALSE)$acf[2], 3), "\n")

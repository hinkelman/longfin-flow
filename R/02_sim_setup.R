# Build the model inputs ("setups") for age-0 and age-1+ from the authors' own
# pipeline. The authors' MARSS script has an age-0 section followed by an age-1+
# section; each is sourced in turn so the station set, missing-data pattern, Z6
# habitat grouping, and R structure are exactly theirs. Each setup also holds the
# fitted Z6 parameters and the z-scored Jan-Jun outflow covariate.
#   output/setups.rds: list(age0, age1); the simulations use age0
# Run after R/00_run_original.R 1 (which creates output/original_run/script1/).

suppressMessages(library(MARSS))
source("R/design.R")

root <- normalizePath(".")
run_dir <- file.path(root, "output", "original_run", "script1")
script1 <- readLines(list.files(file.path(root, "original", "code"), pattern = "^1_",
                                full.names = TRUE))
age1_start <- grep("^#+AGE-1\\+#+$", script1)
stopifnot(length(age1_start) == 1)

dayflow <- read.csv(file.path(root, "original", "data", "annual_average_deltawide_flow_jan_june.csv"))
flow <- dayflow$Mean_outflow_CY[match(YEARS, dayflow$Calendar_Year)]
flow_z <- as.numeric(zscore(flow))

# Source lines of the authors' script into env (in its run directory) and collect
# the pieces the models need
source_section <- function(lines, env) {
  tmp <- tempfile(fileext = ".R")
  writeLines(lines, tmp)
  old_wd <- setwd(run_dir)
  on.exit(setwd(old_wd))
  pdf(NULL)  # the sections draw plots; discard them
  source(tmp, local = env)
  dev.off()
  fit <- env$z6
  list(years = YEARS,
       flow_z = flow_z,
       y_obs = env$MWT_OT_matrix2,          # observed log CPUE (series x years)
       mask = !is.na(env$MWT_OT_matrix2),   # real missing-data pattern
       Z = env$Z6,                          # series -> state grouping (factor)
       R_model = env$R.model,               # gear-specific observation variances
       gear = env$z_matrix$Gear,
       Q = setNames(coef(fit)$Q[, 1], rownames(fit$states)),
       R = coef(fit)$R[, 1],
       x0 = setNames(coef(fit)$x0[, 1], rownames(fit$states)),
       AICc = fit$AICc)
}

# the age-1+ section reuses packages and helper objects from the age-0 section,
# so both run in the same environment
env <- new.env()
setups <- list(age0 = source_section(script1[seq_len(age1_start - 1)], env),
               age1 = source_section(script1[age1_start:length(script1)], env))
saveRDS(setups, file.path(root, "output", "setups.rds"))

for (age in names(setups)) {
  s <- setups[[age]]
  cat("\n", age, ": series", nrow(s$y_obs), " years", ncol(s$y_obs),
      " fraction observed", round(mean(s$mask), 3), " Z6 AICc", round(s$AICc, 2), "\n")
  print(table(s$Z))
}
cat("\nPaper's Table 1 AICc for Z6: age-0 803.29, age-1+ 248.84\n")
cat("Flow lag-1 autocorrelation:", round(acf(flow_z, plot = FALSE)$acf[2], 3), "\n")

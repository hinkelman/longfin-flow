# Run a simulation study in parallel and save replicate-level results.
# Each replicate is checkpointed to output/sim_reps/<study>/, so a rerun skips
# finished replicates; progress is appended to output/sim_progress_<study>.log.
# Usage: Rscript R/03_run_simulation.R [n_reps] [n_workers] [study]
#
# Studies:
#   main          truth = growth (the paper's MARSS model); fits the paper's MARSS C,
#                 the paper's DLM, and the era model
#   dlm_variants  same data as main (same seeds); fits the DLM variants that
#                 change one modeling choice at a time (fit_dlm_variants)
#   level_truth   truth = level (the paper's DLM model); fits everything

suppressMessages({
  library(future.apply)
})
source("R/sim_functions.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100
n_workers <- if (length(args) >= 2) as.integer(args[2]) else max(1, parallelly::availableCores() - 1)
study <- if (length(args) >= 3) args[3] else "main"

study_opts <- switch(study,
  main = list(truth = "growth", main = TRUE, variants = FALSE, seed_offset = 0,
              out = "sim_results.rds"),
  dlm_variants = list(truth = "growth", main = FALSE, variants = TRUE, seed_offset = 0,
                      out = "sim_results_dlm_variants.rds"),
  level_truth = list(truth = "level", main = TRUE, variants = TRUE, seed_offset = 10000,
                     out = "sim_results_level_truth.rds"),
  stop("unknown study: ", study))

setup <- readRDS(file.path("output", "sim_setup.rds"))
setup$fit <- NULL  # not needed on workers

rep_dir <- file.path("output", "sim_reps", study)
dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
progress_log <- file.path("output", paste0("sim_progress_", study, ".log"))

# b = 0.08 is roughly the largest time-averaged MARSS flow effect in the paper
# (age-0 channel midwater); b = 0.24 is a 3x "strong signal" check
grid <- expand.grid(rep = seq_len(n_reps),
                    scenario = c("null", "constant", "decline", "step_1987", "step_2002"),
                    b = c(0.08, 0.24),
                    stringsAsFactors = FALSE)
grid <- grid[!(grid$scenario == "null" & grid$b == 0.24), ]  # null does not depend on b
grid$file <- file.path(rep_dir, sprintf("%s_b%.2f_rep%03d.rds", grid$scenario, grid$b, grid$rep))
# replicate-specific seeds so results don't depend on scheduling or restarts;
# dlm_variants shares seeds (and so simulated data) with main
grid$seed <- 20261112 + study_opts$seed_offset + seq_len(nrow(grid))
todo <- which(!file.exists(grid$file))
cat("Study:", study, " replicates:", nrow(grid), " remaining:", length(todo), "\n")

plan(multisession, workers = n_workers)
t0 <- Sys.time()
invisible(future_lapply(todo, function(i) {
  g <- grid[i, ]
  set.seed(g$seed)
  t_start <- Sys.time()
  r <- tryCatch(run_replicate(setup, g$scenario, g$b, g$rep, truth = study_opts$truth,
                              main = study_opts$main, variants = study_opts$variants),
                error = function(e) list(error = conditionMessage(e)))
  secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (!is.null(r$summary)) r$summary$secs <- secs
  saveRDS(r, g$file)
  cat(sprintf("%s %s b=%.2f rep=%d %.0fs %s\n", format(Sys.time(), "%H:%M:%S"), g$scenario,
              g$b, g$rep, secs, if (is.null(r$error)) "ok" else paste("ERROR:", r$error)),
      file = progress_log, append = TRUE)
  NULL
}, future.seed = NULL, future.chunk.size = 1))
plan(sequential)

results <- lapply(grid$file, readRDS)
ok <- vapply(results, function(r) is.null(r$error), logical(1))
cat("Failed replicates:", sum(!ok), "\n")
bind <- function(part) {
  parts <- lapply(results[ok], `[[`, part)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts)) do.call(rbind, parts) else NULL
}
saveRDS(list(summary = bind("summary"), paths = bind("paths"),
             variants = bind("variants"), variant_paths = bind("variant_paths"), grid = grid),
        file.path("output", study_opts$out))
cat("Elapsed:", format(Sys.time() - t0), "\n")

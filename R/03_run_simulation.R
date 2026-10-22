# Run the simulation study in parallel and save replicate-level results.
# Each replicate is checkpointed to output/sim_reps/, so a rerun skips finished
# replicates; progress is appended to output/sim_progress.log.
# Usage: Rscript R/03_run_simulation.R [n_reps] [n_workers]

suppressMessages({
  library(future.apply)
})
source("R/sim_functions.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100
n_workers <- if (length(args) >= 2) as.integer(args[2]) else max(1, parallelly::availableCores() - 1)

setup <- readRDS(file.path("output", "sim_setup.rds"))
setup$fit <- NULL  # not needed on workers

rep_dir <- file.path("output", "sim_reps")
dir.create(rep_dir, showWarnings = FALSE)
progress_log <- file.path("output", "sim_progress.log")

# b = 0.08 is roughly the largest time-averaged MARSS flow effect in the paper
# (age-0 channel midwater); b = 0.24 is a 3x "strong signal" check
grid <- expand.grid(rep = seq_len(n_reps),
                    scenario = c("null", "constant", "decline", "step_1987", "step_2002"),
                    b = c(0.08, 0.24),
                    stringsAsFactors = FALSE)
grid <- grid[!(grid$scenario == "null" & grid$b == 0.24), ]  # null does not depend on b
grid$file <- file.path(rep_dir, sprintf("%s_b%.2f_rep%03d.rds", grid$scenario, grid$b, grid$rep))
# replicate-specific seeds so results don't depend on scheduling or restarts
grid$seed <- 20261112 + seq_len(nrow(grid))
todo <- which(!file.exists(grid$file))
cat("Replicates:", nrow(grid), " remaining:", length(todo), "\n")

plan(multisession, workers = n_workers)
t0 <- Sys.time()
invisible(future_lapply(todo, function(i) {
  g <- grid[i, ]
  set.seed(g$seed)
  t_start <- Sys.time()
  r <- tryCatch(run_replicate(setup, g$scenario, g$b, g$rep),
                error = function(e) list(error = conditionMessage(e)))
  secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (is.null(r$error)) r$summary$secs <- secs
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
summary_df <- do.call(rbind, lapply(results[ok], `[[`, "summary"))
paths_df <- do.call(rbind, lapply(results[ok], `[[`, "paths"))
saveRDS(list(summary = summary_df, paths = paths_df, grid = grid),
        file.path("output", "sim_results.rds"))
cat("Elapsed:", format(Sys.time() - t0), "\n")

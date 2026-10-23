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
#   checks_growth same data as main; structure (growth vs level) and interval checks
#   checks_level  same data as level_truth; structure and interval checks

suppressMessages({
  library(future.apply)
})
source("R/sim_functions.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100
n_workers <- if (length(args) >= 2) as.integer(args[2]) else max(1, parallelly::availableCores() - 1)
study <- if (length(args) >= 3) args[3] else "main"

# seed_offset: studies with the same truth share seeds, so fit identical data
study_opts <- switch(study,
  main = list(truth = "growth", main = TRUE, variants = FALSE, checks = FALSE, seed_offset = 0),
  dlm_variants = list(truth = "growth", main = FALSE, variants = TRUE, checks = FALSE, seed_offset = 0),
  level_truth = list(truth = "level", main = TRUE, variants = TRUE, checks = FALSE, seed_offset = 10000),
  checks_growth = list(truth = "growth", main = FALSE, variants = FALSE, checks = TRUE, seed_offset = 0),
  checks_level = list(truth = "level", main = FALSE, variants = FALSE, checks = TRUE, seed_offset = 10000),
  stop("unknown study: ", study))

setup <- readRDS(file.path("output", "setups.rds"))$age0

rep_dir <- file.path("output", "sim_reps", study)
dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
progress_log <- file.path("output", paste0("sim_progress_", study, ".log"))

grid <- sim_grid(n_reps)
grid$file <- file.path(rep_dir, sprintf("%s_b%.2f_rep%03d.rds", grid$scenario, grid$b, grid$rep))
grid$seed <- grid_seed(seq_len(nrow(grid)), study_opts$seed_offset)
todo <- which(!file.exists(grid$file))
cat("Study:", study, " replicates:", nrow(grid), " remaining:", length(todo), "\n")

plan(multisession, workers = n_workers)
t0 <- Sys.time()
invisible(future_lapply(todo, function(i) {
  g <- grid[i, ]
  set.seed(g$seed)
  t_start <- Sys.time()
  r <- tryCatch(run_replicate(setup, g$scenario, g$b, g$rep, truth = study_opts$truth,
                              main = study_opts$main, variants = study_opts$variants,
                              checks = study_opts$checks),
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
# Compact, committed results (output/results/<study>_*.csv): replicate-level
# summaries, plus slope paths reduced to quantiles across replicates. The
# summary scripts read only these, so they run without the checkpoints.
results_dir <- file.path("output", "results")
dir.create(results_dir, showWarnings = FALSE)
write_result <- function(df, name) {
  if (is.null(df)) return(invisible())
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], signif, digits = 6)
  write.csv(df, file.path(results_dir, paste0(study, "_", name, ".csv")), row.names = FALSE)
}
path_quantiles <- function(paths, keys, value) {
  if (is.null(paths)) return(NULL)
  groups <- split(paths, paths[keys], drop = TRUE)
  do.call(rbind, lapply(groups, function(d) {
    q <- quantile(d[[value]], c(0.1, 0.5, 0.9), names = FALSE)
    cbind(d[1, keys, drop = FALSE], true_beta = d$true_beta[1], q10 = q[1], median = q[2], q90 = q[3])
  }))
}
write_result(bind("summary"), "summary")
write_result(bind("variants"), "variants")
write_result(bind("structure"), "structure")
write_result(path_quantiles(bind("paths"), c("scenario", "b", "year"), "dlm_slope"),
             "dlm_path_quantiles")
write_result(path_quantiles(bind("variant_paths"), c("scenario", "b", "variant", "year"), "slope"),
             "variant_path_quantiles")
cat("Elapsed:", format(Sys.time() - t0), "\n")

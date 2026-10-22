# Run the authors' Zenodo scripts unmodified against the Dryad data.
# Each script runs in its own scratch directory (output/original_run/script<N>/)
# holding copies of the Dryad files, so the original folders stay untouched and
# scripts can run in parallel without sharing files.
# Usage: Rscript R/00_run_original.R <script-number>

args <- commandArgs(trailingOnly = TRUE)
script_num <- if (length(args)) args[1] else "1"

root <- normalizePath(".")
run_dir <- file.path(root, "output", "original_run", paste0("script", script_num))
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
file.copy(list.files(file.path(root, "original", "data"), full.names = TRUE),
          run_dir, overwrite = TRUE)
# the scripts read "z6_states_age1+.csv"; Dryad names the file "z6_states_age1_.csv"
file.copy(file.path(run_dir, "z6_states_age1_.csv"),
          file.path(run_dir, "z6_states_age1+.csv"), overwrite = TRUE)

script <- list.files(file.path(root, "original", "code"),
                     pattern = paste0("^", script_num, "_"), full.names = TRUE)
setwd(run_dir)
pdf("plots.pdf")
env <- new.env()
t0 <- Sys.time()
source(script, local = env)
dev.off()
cat("Elapsed:", format(Sys.time() - t0), "\n")
saveRDS(as.list(env), "env.rds")

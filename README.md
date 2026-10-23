# longfin-flow

**Static or shifting? Stress-testing time-varying flow–ecology models.**
Material for the IEP Data Science PWT meeting, 12 November 2026.

This project (1) reruns the analysis in Saffarinia et al. (2026, *Ecological Applications* 36:e70178,
<https://doi.org/10.1002/eap.70178>) using the authors' published code and data, and (2) simulates
Longfin Smelt–like data where the true effect of outflow is known (constant, slowly declining, or a
step change at 1987 or 2002) to see what the paper's MARSS and DLM approaches recover.

## Source code and data

This repository does not redistribute the authors' files. Get them from the archives the authors
chose, and cite those archives:

- **Code:** Saffarinia et al. (2025b), Zenodo, <https://doi.org/10.5281/zenodo.14503433> (CC-BY 4.0)
- **Data:** Saffarinia et al. (2025a), Dryad, <https://doi.org/10.5061/dryad.x3ffbg7w0> (CC0; version
  2.0.0, 2026-01-09). The authors ask re-users to contact them about possible collaboration.

`original-manifest.csv` records the md5 checksum of every file this analysis used. The Zenodo checksums
match the ones Zenodo publishes. `R/fetch_original.R` downloads the code from Zenodo and verifies
everything against the manifest. Dryad blocks scripted downloads, so download the dataset zip from the
DOI page in a browser and unzip it into `original/data/` first. If a later Dryad or Zenodo version
changes a file, the script flags the mismatch.

## Layout

| Path | Contents |
|---|---|
| `original-manifest.csv` | Checksums of the authors' files used here |
| `R/fetch_original.R` | Downloads the code from Zenodo; verifies code and data against the manifest |
| `R/00_run_original.R` | Runs one of the authors' scripts unmodified, in `output/original_run/script<N>/` |
| `R/01_replication_check.R` | Compares rerun to published coefficients; DLM variance diagnostics |
| `R/design.R` | Study design shared by all scripts: years, scenarios, eras, replicate grid and seeds, thresholds |
| `R/02_sim_setup.R` | Builds the age-0 and age-1+ model inputs (stations, missing-data pattern, habitat grouping, fitted variances) from the authors' pipeline |
| `R/sim_functions.R` | Data simulation, model fitting, and scoring functions |
| `R/plot_helpers.R` | Shared labels, colors, theme, and result readers for the summary scripts |
| `R/03_run_simulation.R` | Runs one simulation study in parallel (checkpointed, resumable); see the study list below |
| `R/04_summarize_simulation.R` | Recovery tables and figures (study `main`) |
| `R/05_followup_checks.R` | What static MARSS C averages over; era-model false positives and calibrated threshold |
| `R/06_summarize_extensions.R` | Why the DLM shrinks the flow effect; results when the DLM's own model is the truth |
| `R/07_real_data.R` | Fits every model form (growth vs level, static, era, time-varying) to the real age-0 and age-1+ data |
| `R/08_summarize_checks.R` | Can AICc tell growth from level? Does a prior on x0 fix the slope intervals? |
| `output/results/` | Replicate-level simulation results and slope-path quantiles (the summary scripts read only these) |
| `output/` | Summary tables (`*.csv`) and figures (`figures/`); model objects and checkpoints are not tracked |

Simulation studies (`Rscript R/03_run_simulation.R <n_reps> <n_workers> <study>`):

| Study | Truth | What it fits | Hours on 7 cores |
|---|---|---|---|
| `main` | flow drives growth | paper's MARSS C, paper's DLM, era model | 1.5 |
| `dlm_variants` | growth (same data as `main`) | DLM variants changing one choice at a time | 1.7 |
| `level_truth` | flow drives level | everything in `main` and `dlm_variants` | 2.5 |
| `checks_growth` | growth (same data as `main`) | growth vs level structure; interval checks | 1.8 |
| `checks_level` | level (same data as `level_truth`) | growth vs level structure; interval checks | 1.1 |

## Reproduce

Requires R 4.6.1; package versions are pinned in `renv.lock`. Run `renv::restore()` once to install them.

To regenerate the tables and figures from the committed simulation results (minutes):

```bash
Rscript R/fetch_original.R
for i in 1 2 3 4; do Rscript R/00_run_original.R $i; done
Rscript R/01_replication_check.R
Rscript R/02_sim_setup.R
Rscript R/04_summarize_simulation.R
Rscript R/05_followup_checks.R
Rscript R/06_summarize_extensions.R
Rscript R/07_real_data.R
Rscript R/08_summarize_checks.R
```

To rerun the simulations themselves (about 9.5 hours on 7 cores), run the studies before the summary
scripts:

```bash
for s in main dlm_variants level_truth checks_growth checks_level; do Rscript R/03_run_simulation.R 100 7 $s; done
```

Each replicate is seeded on its own and checkpointed to `output/sim_reps/<study>/`, so an interrupted run
picks up where it stopped. Studies on the same truth share seeds, so they fit models to identical
simulated data. Each study writes its compact results to `output/results/<study>_*.csv`.
`02_sim_setup.R` writes `output/setups.rds` (not tracked), which the simulations, `05_followup_checks.R`
(it refits a few null replicates), and `07_real_data.R` read.

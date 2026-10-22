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
| `R/02_sim_setup.R` | Pulls the age-0 station layout, missing-data pattern, and fitted variances from the authors' pipeline |
| `R/sim_functions.R` | Data simulation and model-fitting functions |
| `R/03_run_simulation.R` | Runs the simulation in parallel (checkpointed, resumable) |
| `R/04_summarize_simulation.R` | Recovery tables and figures |
| `output/` | Result tables (`*.csv`) and figures (`figures/`); model objects and checkpoints are not tracked |

## Reproduce

Requires R (tested on 4.6.1) with `MARSS` (3.11.10), `tidyverse`, `viridis`, `zoo`, `future.apply`.

```bash
Rscript R/fetch_original.R
for i in 1 2 3 4; do Rscript R/00_run_original.R $i; done
Rscript R/01_replication_check.R
Rscript R/02_sim_setup.R
Rscript R/03_run_simulation.R 100 7
Rscript R/04_summarize_simulation.R
```

The simulation (900 replicates) takes about 1.5 hours on 7 cores. Each replicate is seeded on its own
and checkpointed to `output/sim_reps/`, so an interrupted run picks up where it stopped.

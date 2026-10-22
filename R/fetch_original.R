# Get the authors' code and data from their archives and verify them against
# original-manifest.csv (md5 checksums of the versions this analysis used).
# The files are not redistributed in this repository (original/ is gitignored).
#
#   Code: Zenodo, https://doi.org/10.5281/zenodo.14503433 (CC-BY 4.0)
#         downloaded automatically through the Zenodo API
#   Data: Dryad, https://doi.org/10.5061/dryad.x3ffbg7w0 (CC0; version 2.0.0, 2026-01-09)
#         Dryad blocks scripted downloads, so download the zip from the DOI page in a
#         browser and unzip it into original/data/; this script then verifies it
#
# Usage: Rscript R/fetch_original.R

manifest <- read.csv("original-manifest.csv")
code_dir <- file.path("original", "code")
data_dir <- file.path("original", "data")
dir.create(code_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# Code from Zenodo ---------------------------------------------------------------
zenodo_files <- manifest$file[manifest$source == "zenodo"]
for (f in zenodo_files) {
  dest <- file.path(code_dir, f)
  if (!file.exists(dest)) {
    url <- sprintf("https://zenodo.org/api/records/14503433/files/%s/content", utils::URLencode(f))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
}

# Verify everything ----------------------------------------------------------------
manifest$path <- file.path(ifelse(manifest$source == "zenodo", code_dir, data_dir), manifest$file)
manifest$present <- file.exists(manifest$path)
manifest$md5_ok <- FALSE
manifest$md5_ok[manifest$present] <-
  unname(tools::md5sum(manifest$path[manifest$present])) == manifest$md5[manifest$present]

missing_data <- manifest$source == "dryad" & !manifest$present
if (any(missing_data)) {
  cat("Dryad data not found in original/data/. Download the dataset zip from\n",
      "  https://doi.org/10.5061/dryad.x3ffbg7w0\n",
      "unzip it into original/data/, and rerun this script.\n", sep = "")
}
bad <- manifest$present & !manifest$md5_ok
if (any(bad)) {
  cat("Checksum mismatch (file differs from the version this analysis used):\n")
  print(manifest[bad, c("source", "file")], row.names = FALSE)
}
cat(sprintf("Verified %d of %d files.\n", sum(manifest$md5_ok), nrow(manifest)))
if (!all(manifest$md5_ok)) quit(status = 1)

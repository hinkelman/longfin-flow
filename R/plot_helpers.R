# Shared labels, colors, and theme for the summary figures.

suppressMessages({
  library(dplyr)
  library(ggplot2)
})
source("R/design.R")

FIG_DIR <- file.path("output", "figures")
RESULTS_DIR <- file.path("output", "results")
dir.create(FIG_DIR, showWarnings = FALSE)

# Categorical palette (validated default order) and the true-effect ink
COL <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
         magenta = "#e87ba4", violet = "#6250d6", truth = "#0b0b0b")

TRUTH_LABELS <- c(growth = "Truth: flow drives growth", level = "Truth: flow drives level")

# Read one study's compact results (written by R/03_run_simulation.R)
read_result <- function(study, part) {
  read.csv(file.path(RESULTS_DIR, paste0(study, "_", part, ".csv")))
}

# Order scenarios for tables and facets
factor_scenario <- function(x) factor(x, names(SCENARIO_LABELS))

# Facet label for the effect size; the null has no b
b_label <- function(scenario, b) ifelse(scenario == "null", "no effect", sprintf("b = %.2f", b))

theme_talk <- function() {
  theme_minimal(base_size = 12) +
    theme(legend.position = "top", plot.title.position = "plot")
}

save_figure <- function(name, plot, width, height) {
  ggsave(file.path(FIG_DIR, name), plot, width = width, height = height, dpi = 200, bg = "white")
}

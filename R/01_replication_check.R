# Compare our rerun of the authors' scripts (R/00_run_original.R 1-4) with the
# published outputs on Dryad, and pull out DLM diagnostics used in the talk.
# Run from the project root after R/00_run_original.R has been run for 1-4.

suppressMessages({
  library(MARSS)
  library(dplyr)
  library(tidyr)
})

env <- lapply(setNames(1:4, paste0("s", 1:4)), function(i)
  readRDS(file.path("output", "original_run", paste0("script", i), "env.rds")))

# MARSS flow coefficients (C) ---------------------------------------------
# ageclass labels follow the authors' combined file: "1_cont" is the lagged model
marss_rerun <- bind_rows(
  `0` = env$s2$age_0_marss_covariate_parameters,
  `1` = env$s2$age_1_marss_covariate_parameters,
  `1_cont` = env$s2$age_1_contemporary_marss_covariate_parameters,
  .id = "ageclass"
) |>
  mutate(ChanShoal = as.integer(substr(Parameter, 3, 3)),
         gear = sub("C\\.[12] ", "", Parameter)) |>
  select(ageclass, ChanShoal, gear, marss_rerun = Coefficient,
         marss_rerun_lo = LowerCI, marss_rerun_hi = UpperCI)

# DLM time-averaged slopes --------------------------------------------------
dlm_rerun <- bind_rows(
  `0` = env$s3$DLM_Join_averaged_age0,
  `1` = env$s3$DLM_Join_averaged_age1,
  `1_cont` = env$s4$DLM_Join_averaged_age1_lag,
  .id = "ageclass"
) |>
  mutate(ChanShoal = as.integer(ChanShoal)) |>
  select(ageclass, ChanShoal, gear, dlm_rerun = DLM_slope, dlm_rerun_se = DLM_SE)

published <- read.csv(file.path("original", "data", "marss_DLM_covariate_coefficients.csv")) |>
  mutate(ageclass = as.character(ageclass)) |>
  select(ageclass, ChanShoal, gear, dlm_pub = DLM_slope, dlm_pub_se = DLM_SE,
         marss_pub = marss_coefficient, marss_pub_lo = LowerCI_marss,
         marss_pub_hi = UpperCI_marss)

comparison <- published |>
  left_join(marss_rerun, by = c("ageclass", "ChanShoal", "gear")) |>
  left_join(dlm_rerun, by = c("ageclass", "ChanShoal", "gear")) |>
  mutate(marss_abs_diff = abs(marss_rerun - marss_pub),
         dlm_abs_diff = abs(dlm_rerun - dlm_pub),
         same_sign = sign(marss_pub) == sign(dlm_pub),
         dlm_to_marss_ratio = dlm_pub / marss_pub)

write.csv(comparison, file.path("output", "replication_comparison.csv"), row.names = FALSE)

cat("\nMax |MARSS C rerun - published|:", signif(max(comparison$marss_abs_diff), 3), "\n")
cat("Max |DLM mean slope rerun - published|:", signif(max(comparison$dlm_abs_diff), 3), "\n")
cat("Sign agreement MARSS vs DLM (published):", sum(comparison$same_sign), "of",
    nrow(comparison), "\n\n")
print(comparison |> select(ageclass, ChanShoal, gear, marss_pub, dlm_pub, same_sign,
                           dlm_to_marss_ratio), digits = 3)

# DLM variance diagnostics ----------------------------------------------------
dlm_fits <- list(age0 = env$s3$dlmage0, age1 = env$s3$dlmage1, age1_lag = env$s4$dlmage1lag)
dlm_diag <- bind_rows(lapply(names(dlm_fits), function(nm) {
  f <- dlm_fits[[nm]]
  q <- coef(f)$Q[, 1]
  r <- coef(f)$R[, 1]
  data.frame(model = nm,
             group = sub("q_slo_", "", grep("q_slo", names(q), value = TRUE)),
             q_slope = q[grep("q_slo", names(q))],
             r_obs = unname(r),
             mean_slope_se = rowMeans(f$states.se[5:8, ]))
}))
rownames(dlm_diag) <- NULL
write.csv(dlm_diag, file.path("output", "replication_dlm_diagnostics.csv"), row.names = FALSE)
cat("\nDLM variance estimates (slope process variance, observation variance, mean slope SE):\n")
print(dlm_diag, digits = 3)

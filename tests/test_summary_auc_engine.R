source("R/summary_auc_engine.R")

stopifnot(abs(pilot_summary$residual_sd - 0.394) < 0.01)
stopifnot(abs(pilot_summary$baseline_slope - 0.976) < 0.02)

set.seed(1)
dat_null <- simulate_summary_trial(
  n_per_arm = 50000,
  reduction = 0,
  allocation_kenya = 0.5
)
set.seed(1)
dat_effect <- simulate_summary_trial(
  n_per_arm = 50000,
  reduction = 0.25,
  allocation_kenya = 0.5
)

# Identical seeds guarantee identical baseline draws. Only the randomized
# outcome is shifted, so treatment cannot alter baseline.
stopifnot(identical(dat_null$baseline_c, dat_effect$baseline_c))

observed_shift <- with(
  dat_effect,
  mean(log_sauc[treatment == 1]) - mean(log_sauc[treatment == 0])
)
stopifnot(abs(observed_shift - log(0.75)) < 0.02)

fit_null <- fit_summary_bayes(
  dat_null,
  prior_method = "robust",
  historical_ess = 10,
  robust_weight = 0.5
)
stopifnot(fit_null$prob_benefit > 0, fit_null$prob_benefit < 1)
stopifnot(
  fit_null$historical_posterior_weight >= 0,
  fit_null$historical_posterior_weight <= 1
)

cat("summary_auc_engine tests passed\n")

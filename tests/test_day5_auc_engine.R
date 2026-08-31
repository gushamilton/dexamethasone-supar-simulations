source("R/day5_auc_engine.R")

stopifnot(all(eigen(pilot_day5$covariance)$values > 0))
stopifnot(identical(pilot_day5$times, c(0, 1, 3, 5)))
stopifnot(
  abs(
    cov2cor(pilot_day5$covariance)[3, 4] -
      pilot_day5$projected_day3_day5_correlation
  ) < 1e-10
)
stopifnot(
  abs(log_auc_delta_day5(c(pilot_day5$mean, 0))) < 1e-10
)
stopifnot(identical(make_day5_prior()$covariance[5, 5], 1))

calibrated_post <- calibrate_day5_effect(0.20, "immediate")
calibrated_shift <- effect_vector_day5(calibrated_post, "immediate")
achieved_reduction <- 1 -
  standardised_auc_day5(pilot_day5$mean + calibrated_shift) /
  standardised_auc_day5(pilot_day5$mean)
stopifnot(abs(achieved_reduction - 0.20) < 1e-8)

set.seed(22)
large_trial <- simulate_day5_trial(
  n_per_arm = 1500,
  target_auc_reduction = 0.25,
  missing_mechanism = "complete"
)
baseline_means <- aggregate(
  log_supar ~ treatment,
  large_trial[large_trial$time_index == 1, ],
  mean
)
stopifnot(abs(diff(baseline_means$log_supar)) < 0.04)

large_fit <- fit_day5_auc_bayes(large_trial)
stopifnot(large_fit$prob_benefit > 0.99)
stopifnot(large_fit$posterior_mean < 0)
stopifnot(is.finite(large_fit$posterior_sd), large_fit$posterior_sd > 0)
stopifnot(
  is.finite(large_fit$posterior_treatment_mean),
  is.finite(large_fit$posterior_treatment_sd),
  large_fit$posterior_treatment_sd > 0
)
stopifnot(identical(
  large_fit$success,
  large_fit$prob_benefit >= 0.95
))

missing_trial <- simulate_day5_trial(
  n_per_arm = 45,
  target_auc_reduction = 0.20,
  missing_mechanism = "mar_improvement",
  day3_missing_rate = 0.50,
  seed = 23
)
stopifnot(all(missing_trial$observed[missing_trial$time_index != 3]))
stopifnot(any(!missing_trial$observed[missing_trial$time_index == 3]))

auc_fit <- fit_day5_auc_bayes(missing_trial)
day3_fit <- fit_day3_ancova_bayes(missing_trial)
stopifnot(is.finite(auc_fit$posterior_mean))
stopifnot(is.finite(day3_fit$posterior_mean))

cat("day5_auc_engine tests passed\n")

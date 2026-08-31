source("R/nuisance_borrowing_engine.R")

stopifnot(abs(pilot_nuisance$residual_sd - 0.3946) < 0.01)
stopifnot(abs(pilot_nuisance$baseline_slope - 0.976) < 0.02)

set.seed(1)
null_data <- simulate_nuisance_trial(n_per_arm = 500, reduction = 0)
null_fit <- fit_nuisance_bayes(null_data, prior_type = "pilot_conservative")
stopifnot(abs(null_fit$treatment_mean) < 0.08)

set.seed(2)
effect_data <- simulate_nuisance_trial(n_per_arm = 500, reduction = 0.20)
effect_fit <- fit_nuisance_bayes(effect_data, prior_type = "pilot_conservative")
stopifnot(effect_fit$treatment_mean < -0.12)
stopifnot(effect_fit$prob_benefit > 0.95)

small_oc <- simulate_nuisance_oc(
  n_sim = 100,
  seed = 3,
  n_per_arm = 45,
  reduction = 0.20,
  prior_type = "pilot_conservative"
)
stopifnot(small_oc$probability_success > 0.5)
stopifnot(small_oc$coverage_95 > 0.8)

set.seed(4)
joint_data <- simulate_nuisance_trial(n_per_arm = 100, reduction = 0.20)
joint_none <- fit_joint_borrowing_bayes(
  joint_data,
  mean_borrowing = "none"
)
joint_robust <- fit_joint_borrowing_bayes(
  joint_data,
  mean_borrowing = "robust",
  historical_ess = 5,
  robust_weight = 0.25
)
stopifnot(joint_none$historical_posterior_weight == 0)
stopifnot(joint_robust$historical_posterior_weight > 0)
stopifnot(joint_robust$historical_posterior_weight < 1)
stopifnot(joint_robust$prob_benefit > 0.95)

comm_grid <- make_commensurability_grid(
  max_ess = 10,
  gamma_shape = 0.5,
  gamma_rate = 0.1,
  n_grid = 31
)
stopifnot(abs(sum(comm_grid$prior_weight) - 1) < 1e-10)
stopifnot(all(comm_grid$ess > 0), max(comm_grid$ess) <= 10 + 1e-8)

# With a large current sample, compatible controls should retain more
# commensurate information than markedly conflicting controls.
set.seed(5)
compatible_data <- simulate_nuisance_trial(
  n_per_arm = 300,
  reduction = 0,
  current_control_drift = 0
)
set.seed(5)
conflicting_data <- simulate_nuisance_trial(
  n_per_arm = 300,
  reduction = 0,
  current_control_drift = log(1.5)
)
compatible_fit <- fit_commensurate_bayes(compatible_data, max_ess = 10)
conflicting_fit <- fit_commensurate_bayes(conflicting_data, max_ess = 10)
stopifnot(compatible_fit$posterior_mean_ess > conflicting_fit$posterior_mean_ess)
stopifnot(compatible_fit$posterior_mean_ess > 0)

switch_grid <- make_commensurability_grid(
  max_ess = 10,
  prior_family = "scaled_beta",
  beta_shape1 = 0.2,
  beta_shape2 = 0.2,
  n_grid = 31
)
stopifnot(abs(sum(switch_grid$prior_weight) - 1) < 1e-10)
stopifnot(min(switch_grid$ess) > 0, max(switch_grid$ess) < 10)

cat("nuisance_borrowing_engine tests passed\n")

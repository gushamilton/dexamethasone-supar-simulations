library(tidyverse)

source("R/nuisance_borrowing_engine.R")

n_sim <- 20000
posterior_cutoff <- 0.95
treatment_prior_sd <- 0.50

# The conjugate nuisance-prior engine fixes the treatment coefficient prior in
# make_nuisance_prior(). This assertion prevents an unnoticed mismatch with the
# prespecified Normal(0, 0.5^2) prior.
pilot_prior <- make_nuisance_prior("pilot_conservative")
implied_treatment_sd <- sqrt(
  pilot_prior$coefficient_covariance_scale[2, 2]
) * pilot_nuisance$residual_sd
stopifnot(abs(implied_treatment_sd - treatment_prior_sd) < 1e-12)

grid <- crossing(
  drift_pct = c(-20, 0, 20),
  reduction = c(0, 0.20),
  analysis = c("Concurrent controls", "Robust mixture: 10%, ESS <= 5")
)

rows <- split(grid, seq_len(nrow(grid)))
results <- bind_rows(parallel::mclapply(
  seq_along(rows),
  function(index) {
    row <- rows[[index]]
    borrowing <- if (
      row$analysis == "Concurrent controls"
    ) "none" else "robust"
    # The seed depends on the data-generating scenario, not the analysis, so
    # the two analyses use matched simulated trials.
    seed <- 220000 +
      1000 * match(row$drift_pct, c(-20, 0, 20)) +
      round(100 * row$reduction)
    result <- simulate_joint_borrowing_oc(
      n_sim = n_sim,
      seed = seed,
      n_per_arm = 45,
      reduction = row$reduction,
      current_control_drift = log(1 + row$drift_pct / 100),
      nuisance_prior_type = "pilot_conservative",
      mean_borrowing = borrowing,
      historical_ess = 5,
      robust_weight = 0.10,
      posterior_cutoff = posterior_cutoff
    )
    bind_cols(row, as_tibble(result))
  },
  mc.cores = min(4L, parallel::detectCores())
))

write_csv(results, "outputs/stage2_borrowing_sensitivity_20000.csv")

summary <- results |>
  dplyr::select(
    drift_pct,
    reduction,
    analysis,
    probability_success,
    mcse_success
  ) |>
  pivot_wider(
    names_from = analysis,
    values_from = c(probability_success, mcse_success)
  ) |>
  mutate(
    absolute_difference =
      `probability_success_Robust mixture: 10%, ESS <= 5` -
      `probability_success_Concurrent controls`
  )

write_csv(summary, "outputs/stage2_borrowing_sensitivity_summary.csv")
print(summary, n = Inf)

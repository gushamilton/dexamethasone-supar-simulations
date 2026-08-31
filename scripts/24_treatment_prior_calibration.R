library(tidyverse)
library(MASS)

source("R/day5_auc_engine.R")

n_sim <- as.integer(Sys.getenv("N_SIM", unset = "50000"))
n_per_arm <- 40
posterior_cutoff <- 0.95
prior_sds <- c(0.50, 0.75, 1.00, 2.00, 10.00)
reductions <- c(0, 0.20)
control_design <- cbind(diag(4), 0)
treatment_design <- control_design
treatment_design[2:4, 5] <- 1
inverse_covariance <- solve(pilot_day5$covariance)

posterior_setup <- function(prior_sd) {
  prior <- make_day5_prior(prior_sd)
  prior_precision <- solve(prior$covariance)
  information <- n_per_arm * (
    crossprod(control_design, inverse_covariance %*% control_design) +
      crossprod(treatment_design, inverse_covariance %*% treatment_design)
  )
  list(
    prior = prior,
    prior_precision = prior_precision,
    covariance = solve(prior_precision + information)
  )
}

setups <- setNames(lapply(prior_sds, posterior_setup), prior_sds)

posterior_probabilities <- function(control_means, treatment_means, setup) {
  control_score <- control_means %*% inverse_covariance %*% control_design
  treatment_score <- treatment_means %*%
    inverse_covariance %*% treatment_design
  prior_score <- drop(setup$prior_precision %*% setup$prior$mean)
  posterior_means <- (
    n_per_arm * (control_score + treatment_score) +
      matrix(prior_score, nrow(control_means), 5, byrow = TRUE)
  ) %*% setup$covariance
  treatment_sd <- sqrt(setup$covariance[5, 5])
  pnorm(0, mean = posterior_means[, 5], sd = treatment_sd)
}

simulate_reduction <- function(reduction, seed) {
  set.seed(seed)
  post_reduction <- calibrate_day5_effect(reduction, "immediate")
  treatment_shift <- effect_vector_day5(post_reduction, "immediate")
  control_means <- MASS::mvrnorm(
    n_sim,
    pilot_day5$mean,
    pilot_day5$covariance / n_per_arm
  )
  treatment_means <- MASS::mvrnorm(
    n_sim,
    pilot_day5$mean + treatment_shift,
    pilot_day5$covariance / n_per_arm
  )
  map_dfr(prior_sds, function(prior_sd) {
    probabilities <- posterior_probabilities(
      control_means,
      treatment_means,
      setups[[as.character(prior_sd)]]
    )
    probability_success <- mean(probabilities >= posterior_cutoff)
    tibble(
      prior_sd = prior_sd,
      reduction = reduction,
      probability_success = probability_success,
      mcse = sqrt(
        probability_success * (1 - probability_success) / n_sim
      ),
      null_calibrated_cutoff = if_else(
        reduction == 0,
        unname(quantile(probabilities, 0.95)),
        NA_real_
      )
    )
  })
}

results <- bind_rows(
  simulate_reduction(0, 250000),
  simulate_reduction(0.20, 250020)
) |>
  mutate(
    n_sim = n_sim,
    n_per_arm = n_per_arm,
    posterior_cutoff = posterior_cutoff,
    decision_rule = "P(Delta < 0 | data) >= 0.95"
  )

write_csv(results, "outputs/treatment_prior_calibration.csv")
print(results, n = Inf)

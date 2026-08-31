library(MASS)

# Design-stage longitudinal distribution on log suPAR.
#
# Days 0, 1 and 3 are taken directly from the n = 29 pilot cohort. The Stage 2
# schedule adds Day 5. As no pilot Day 5 sample is available, Day 5 is projected
# conservatively to have the same mean and marginal SD as Day 3 and correlation
# 0.80 with Day 3. Covariances with earlier measurements follow from the
# corresponding first-order conditional extension, which guarantees a positive
# definite covariance matrix.
pilot_day5 <- local({
  observed_covariance <- matrix(
    c(
      0.2070, 0.2155, 0.1827,
      0.2155, 0.5011, 0.4864,
      0.1827, 0.4864, 0.5637
    ),
    nrow = 3,
    byrow = TRUE
  )
  day3_day5_correlation <- 0.80
  day5_sd <- sqrt(observed_covariance[3, 3])
  multiplier <- day3_day5_correlation *
    day5_sd / sqrt(observed_covariance[3, 3])
  cov_day5_earlier <- multiplier * observed_covariance[3, ]
  covariance <- rbind(
    cbind(observed_covariance, cov_day5_earlier),
    c(cov_day5_earlier, day5_sd^2)
  )

  list(
    times = c(0, 1, 3, 5),
    mean = c(4.6967, 4.3111, 4.1880, 4.1880),
    covariance = covariance,
    observed_pilot_n = 29,
    projected_day3_day5_correlation = day3_day5_correlation
  )
})

log_sum_exp_day5 <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

standardised_auc_day5 <- function(log_means, times = pilot_day5$times) {
  raw_means <- exp(log_means)
  segment_areas <- diff(times) *
    (raw_means[-length(raw_means)] + raw_means[-1]) / 2
  sum(segment_areas) / (max(times) - min(times))
}

log_auc_delta_day5 <- function(beta) {
  control <- beta[1:4]
  treatment <- control + c(0, beta[5:7])
  log(standardised_auc_day5(treatment)) -
    log(standardised_auc_day5(control))
}

numeric_gradient_day5 <- function(fun, x, step = 1e-5) {
  vapply(seq_along(x), function(index) {
    upper <- lower <- x
    upper[index] <- upper[index] + step
    lower[index] <- lower[index] - step
    (fun(upper) - fun(lower)) / (2 * step)
  }, numeric(1))
}

effect_vector_day5 <- function(
    post_reduction,
    shape = c("immediate", "delayed", "waning")) {
  shape <- match.arg(shape)
  full <- log(1 - post_reduction)
  switch(
    shape,
    immediate = c(0, full, full, full),
    delayed = c(0, 0.5 * full, full, full),
    waning = c(0, full, 0.75 * full, 0.5 * full)
  )
}

calibrate_day5_effect <- function(
    target_auc_reduction,
    shape = c("immediate", "delayed", "waning"),
    control_log_mean = pilot_day5$mean) {
  shape <- match.arg(shape)
  if (target_auc_reduction == 0) return(0)
  objective <- function(post_reduction) {
    treatment_shift <- effect_vector_day5(post_reduction, shape)
    achieved <- 1 -
      standardised_auc_day5(control_log_mean + treatment_shift) /
      standardised_auc_day5(control_log_mean)
    achieved - target_auc_reduction
  }
  uniroot(objective, interval = c(0, 0.95), tol = 1e-10)$root
}

simulate_day5_trial <- function(
    n_per_arm = 45,
    target_auc_reduction = 0.20,
    effect_shape = c("immediate", "delayed", "waning"),
    missing_mechanism = c(
      "complete", "mcar", "mar_improvement", "mar_deterioration"
    ),
    day3_missing_rate = 0,
    missing_strength = 1.25,
    seed = NULL) {
  effect_shape <- match.arg(effect_shape)
  missing_mechanism <- match.arg(missing_mechanism)
  if (!is.null(seed)) set.seed(seed)

  n <- 2 * n_per_arm
  treatment <- rep(c(0, 1), each = n_per_arm)
  values <- MASS::mvrnorm(
    n,
    mu = pilot_day5$mean,
    Sigma = pilot_day5$covariance
  )
  post_reduction <- calibrate_day5_effect(
    target_auc_reduction,
    effect_shape
  )
  treatment_shift <- effect_vector_day5(post_reduction, effect_shape)
  values[treatment == 1, ] <- sweep(
    values[treatment == 1, , drop = FALSE],
    2,
    treatment_shift,
    `+`
  )

  observed <- matrix(TRUE, nrow = n, ncol = 4)
  if (missing_mechanism != "complete" && day3_missing_rate > 0) {
    day1_standardised <- as.numeric(scale(values[, 2]))
    intercept <- qlogis(day3_missing_rate)
    probability_day3_missing <- switch(
      missing_mechanism,
      mcar = rep(day3_missing_rate, n),
      mar_improvement = plogis(
        intercept - missing_strength * day1_standardised
      ),
      mar_deterioration = plogis(
        intercept + missing_strength * day1_standardised
      )
    )
    observed[, 3] <- runif(n) > probability_day3_missing
  }

  long <- do.call(rbind, lapply(seq_len(n), function(index) {
    data.frame(
      id = index,
      treatment = treatment[index],
      time_index = 1:4,
      time = pilot_day5$times,
      log_supar = values[index, ],
      observed = observed[index, ]
    )
  }))
  attr(long, "true_delta") <- log_auc_delta_day5(
    c(pilot_day5$mean, treatment_shift[2:4])
  )
  long
}

make_day5_prior <- function(treatment_prior_sd = 0.50) {
  list(
    mean = c(pilot_day5$mean, 0, 0, 0),
    covariance = diag(c(rep(2.00, 4), rep(treatment_prior_sd, 3))^2)
  )
}

fit_day5_auc_bayes <- function(
    long_data,
    treatment_prior_sd = 0.50,
    posterior_cutoff = 0.95) {
  covariance <- pilot_day5$covariance
  prior <- make_day5_prior(treatment_prior_sd)
  prior_precision <- solve(prior$covariance)
  information <- matrix(0, 7, 7)
  score <- numeric(7)

  for (id_value in unique(long_data$id)) {
    subject <- long_data[long_data$id == id_value, ]
    subject <- subject[subject$observed, ]
    indices <- subject$time_index
    design <- matrix(0, nrow(subject), 7)
    design[cbind(seq_len(nrow(subject)), indices)] <- 1
    if (subject$treatment[1] == 1) {
      post_baseline <- subject$time_index > 1
      design[cbind(
        which(post_baseline),
        3 + subject$time_index[post_baseline]
      )] <- 1
    }
    inverse_covariance <- solve(covariance[indices, indices, drop = FALSE])
    information <- information +
      crossprod(design, inverse_covariance %*% design)
    score <- score +
      crossprod(design, inverse_covariance %*% subject$log_supar)[, 1]
  }

  posterior_covariance <- solve(prior_precision + information)
  posterior_mean <- posterior_covariance %*%
    (prior_precision %*% prior$mean + score)
  delta_mean <- log_auc_delta_day5(drop(posterior_mean))
  gradient <- numeric_gradient_day5(
    log_auc_delta_day5,
    drop(posterior_mean)
  )
  delta_sd <- sqrt(drop(crossprod(
    gradient,
    posterior_covariance %*% gradient
  )))
  probability_benefit <- pnorm(0, delta_mean, delta_sd)

  list(
    posterior_mean = delta_mean,
    posterior_sd = delta_sd,
    prob_benefit = probability_benefit,
    success = probability_benefit >= posterior_cutoff
  )
}

fit_day3_ancova_bayes <- function(
    long_data,
    treatment_prior_sd = 0.50,
    posterior_cutoff = 0.95) {
  observed <- long_data[long_data$observed, ]
  baseline <- observed[observed$time_index == 1, c("id", "log_supar")]
  names(baseline)[2] <- "baseline"
  day3 <- observed[observed$time_index == 3, c(
    "id", "treatment", "log_supar"
  )]
  names(day3)[3] <- "outcome"
  analysis <- merge(day3, baseline, by = "id", all = FALSE)
  analysis$baseline_c <- analysis$baseline - pilot_day5$mean[1]

  design <- model.matrix(~ treatment + baseline_c, data = analysis)
  outcome <- analysis$outcome
  covariance_03 <- pilot_day5$covariance[c(1, 3), c(1, 3)]
  baseline_slope <- covariance_03[1, 2] / covariance_03[1, 1]
  residual_sd <- sqrt(
    covariance_03[2, 2] -
      covariance_03[1, 2]^2 / covariance_03[1, 1]
  )
  prior_mean <- c(pilot_day5$mean[3], 0, baseline_slope)
  prior_covariance <- diag(c(2.00, treatment_prior_sd, 1.00)^2)
  prior_precision <- solve(prior_covariance)
  posterior_covariance <- solve(
    prior_precision + crossprod(design) / residual_sd^2
  )
  posterior_mean <- posterior_covariance %*% (
    prior_precision %*% prior_mean +
      crossprod(design, outcome) / residual_sd^2
  )
  treatment_index <- match("treatment", colnames(design))
  delta_mean <- posterior_mean[treatment_index]
  delta_sd <- sqrt(
    posterior_covariance[treatment_index, treatment_index]
  )
  probability_benefit <- pnorm(0, delta_mean, delta_sd)

  list(
    posterior_mean = delta_mean,
    posterior_sd = delta_sd,
    prob_benefit = probability_benefit,
    success = probability_benefit >= posterior_cutoff
  )
}

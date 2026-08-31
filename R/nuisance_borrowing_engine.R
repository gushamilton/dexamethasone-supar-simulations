# Fast conjugate Bayesian engine for nuisance-prior simulations.
#
# The endpoint is baseline-adjusted log standardised AUC. Historical data may
# inform the residual variance and baseline prognostic slope, but never the
# control mean or treatment effect. A normal-inverse-chi-squared approximation
# permits large operating-characteristic simulations without fixing sigma.

source("R/summary_auc_engine.R")

pilot_nuisance <- list(
  n_historical = 29,
  residual_sd = pilot_summary$residual_sd,
  baseline_slope = pilot_summary$baseline_slope,
  baseline_sd = pilot_summary$sd_log_baseline
)

make_nuisance_prior <- function(
    type = c("weak", "pilot_conservative", "pilot_strong")) {
  type <- match.arg(type)

  settings <- switch(
    type,
    weak = list(
      variance_df = 1,
      variance_scale = 0.75,
      baseline_mean = 0,
      baseline_sd = 1.50
    ),
    pilot_conservative = list(
      variance_df = 4,
      variance_scale = pilot_nuisance$residual_sd,
      baseline_mean = pilot_nuisance$baseline_slope,
      baseline_sd = 0.50
    ),
    pilot_strong = list(
      variance_df = 15,
      variance_scale = pilot_nuisance$residual_sd,
      baseline_mean = pilot_nuisance$baseline_slope,
      baseline_sd = 0.20
    )
  )

  # Parameter order: intercept, treatment, centred baseline, Kenya.
  coefficient_mean <- c(0, 0, settings$baseline_mean, 0)
  coefficient_sd <- c(2.0, 1.00, settings$baseline_sd, 0.75)

  # Conditional prior: beta | sigma^2 ~ N(m0, sigma^2 V0).
  # Keep coefficient-prior strength comparable across nuisance-variance
  # scenarios. In particular, the treatment prior remains equivalent to
  # N(0, 1.00^2) for the log-scale treatment coefficient.
  v0 <- diag((coefficient_sd / pilot_nuisance$residual_sd)^2)

  c(settings, list(
    coefficient_mean = coefficient_mean,
    coefficient_covariance_scale = v0,
    label = type
  ))
}

simulate_nuisance_trial <- function(
    n_per_arm = 45,
    reduction = 0.20,
    residual_sd_multiplier = 1,
    baseline_slope_multiplier = 1,
    current_control_drift = 0,
    allocation_kenya = 0.50,
    kenya_baseline_mean_shift_sd = 0,
    kenya_control_shift = 0,
    kenya_residual_sd_ratio = 1,
    kenya_baseline_slope_multiplier = 1,
    kenya_treatment_log_interaction = 0,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_total <- 2 * n_per_arm
  treatment <- rep(c(0, 1), each = n_per_arm)
  n_ke_arm <- round(n_per_arm * allocation_kenya)
  kenya <- unlist(lapply(c(0, 1), function(z) {
    sample(c(rep(1, n_ke_arm), rep(0, n_per_arm - n_ke_arm)))
  }))

  baseline_c <- rnorm(
    n_total,
    mean = 0,
    sd = pilot_nuisance$baseline_sd
  ) + kenya * kenya_baseline_mean_shift_sd * pilot_nuisance$baseline_sd
  true_slope <- pilot_nuisance$baseline_slope * baseline_slope_multiplier
  true_sigma <- pilot_nuisance$residual_sd * residual_sd_multiplier
  true_delta <- log(1 - reduction)

  log_sauc_centered <-
    current_control_drift +
    true_slope * ifelse(
      kenya == 1,
      kenya_baseline_slope_multiplier,
      1
    ) * baseline_c +
    true_delta * treatment +
    kenya_control_shift * kenya +
    kenya_treatment_log_interaction * kenya * treatment +
    rnorm(
      n_total,
      0,
      true_sigma * ifelse(kenya == 1, kenya_residual_sd_ratio, 1)
    )

  data.frame(
    log_sauc_centered = log_sauc_centered,
    treatment = treatment,
    baseline_c = baseline_c,
    kenya = kenya
  )
}

fit_nig_component <- function(x, y, prior, intercept_mean, intercept_sd) {
  m0 <- prior$coefficient_mean
  m0[1] <- intercept_mean
  v0 <- prior$coefficient_covariance_scale
  v0[1, 1] <- (intercept_sd / prior$variance_scale)^2
  p0 <- solve(v0)
  vn <- solve(p0 + crossprod(x))
  mn <- vn %*% (p0 %*% m0 + crossprod(x, y))

  nu_n <- prior$variance_df + nrow(x)
  sum_squares <-
    prior$variance_df * prior$variance_scale^2 +
    sum(y^2) +
    drop(crossprod(m0, p0 %*% m0)) -
    drop(crossprod(mn, solve(vn, mn)))
  sigma2_scale_n <- sum_squares / nu_n

  treatment_index <- match("treatment", colnames(x))
  treatment_mean <- drop(mn[treatment_index])
  treatment_scale <- sqrt(
    sigma2_scale_n * vn[treatment_index, treatment_index]
  )

  a0 <- prior$variance_df / 2
  b0 <- prior$variance_df * prior$variance_scale^2 / 2
  an <- nu_n / 2
  bn <- sum_squares / 2
  log_marginal <-
    -0.5 * as.numeric(determinant(v0, logarithm = TRUE)$modulus) +
    0.5 * as.numeric(determinant(vn, logarithm = TRUE)$modulus) +
    a0 * log(b0) - an * log(bn) + lgamma(an) - lgamma(a0)

  list(
    treatment_mean = treatment_mean,
    treatment_scale = treatment_scale,
    treatment_df = nu_n,
    prob_benefit = pt(
      (0 - treatment_mean) / treatment_scale,
      df = nu_n
    ),
    log_marginal = log_marginal,
    posterior_sigma_mean = sqrt(sum_squares / max(nu_n - 2, 1))
  )
}

fit_joint_borrowing_bayes <- function(
    data,
    nuisance_prior_type = c("pilot_conservative", "weak", "pilot_strong"),
    mean_borrowing = c("none", "robust"),
    historical_ess = 5,
    robust_weight = 0.25,
    posterior_cutoff = 0.95) {
  nuisance_prior_type <- match.arg(nuisance_prior_type)
  mean_borrowing <- match.arg(mean_borrowing)
  prior <- make_nuisance_prior(nuisance_prior_type)
  x <- model.matrix(~ treatment + baseline_c + kenya, data = data)
  y <- data$log_sauc_centered

  weak_component <- fit_nig_component(
    x, y, prior,
    intercept_mean = 0,
    intercept_sd = 2
  )
  weak_component$label <- "weak"

  if (mean_borrowing == "none" || historical_ess <= 0 || robust_weight <= 0) {
    components <- list(weak_component)
    initial_weights <- 1
  } else {
    historical_component <- fit_nig_component(
      x, y, prior,
      intercept_mean = 0,
      intercept_sd = prior$variance_scale / sqrt(historical_ess)
    )
    historical_component$label <- "historical"
    components <- list(historical_component, weak_component)
    initial_weights <- c(robust_weight, 1 - robust_weight)
  }

  log_weights <- log(initial_weights) + vapply(
    components, `[[`, numeric(1), "log_marginal"
  )
  posterior_weights <- exp(log_weights - log_sum_exp(log_weights))
  treatment_means <- vapply(
    components, `[[`, numeric(1), "treatment_mean"
  )
  probability <- sum(posterior_weights * vapply(
    components, `[[`, numeric(1), "prob_benefit"
  ))
  historical_weight <- sum(posterior_weights[vapply(
    components, function(z) z$label == "historical", logical(1)
  )])

  list(
    treatment_mean = sum(posterior_weights * treatment_means),
    prob_benefit = probability,
    success = probability >= posterior_cutoff,
    historical_posterior_weight = historical_weight,
    nominal_borrowing_index = historical_weight * historical_ess,
    posterior_sigma_mean = sum(posterior_weights * vapply(
      components, `[[`, numeric(1), "posterior_sigma_mean"
    ))
  )
}

simulate_joint_borrowing_oc <- function(
    n_sim = 5000,
    seed = 20260717,
    n_per_arm = 45,
    reduction = 0.20,
    residual_sd_multiplier = 1,
    baseline_slope_multiplier = 1,
    current_control_drift = 0,
    nuisance_prior_type = "pilot_conservative",
    mean_borrowing = c("none", "robust"),
    historical_ess = 5,
    robust_weight = 0.25,
    posterior_cutoff = 0.95,
    ...) {
  mean_borrowing <- match.arg(mean_borrowing)
  set.seed(seed)
  true_delta <- log(1 - reduction)
  results <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    dat <- simulate_nuisance_trial(
      n_per_arm = n_per_arm,
      reduction = reduction,
      residual_sd_multiplier = residual_sd_multiplier,
      baseline_slope_multiplier = baseline_slope_multiplier,
      current_control_drift = current_control_drift,
      ...
    )
    fit <- fit_joint_borrowing_bayes(
      dat,
      nuisance_prior_type = nuisance_prior_type,
      mean_borrowing = mean_borrowing,
      historical_ess = historical_ess,
      robust_weight = robust_weight,
      posterior_cutoff = posterior_cutoff
    )
    results[[i]] <- data.frame(
      success = fit$success,
      estimate = fit$treatment_mean,
      historical_weight = fit$historical_posterior_weight,
      borrowing_index = fit$nominal_borrowing_index,
      sigma_estimate = fit$posterior_sigma_mean
    )
  }
  raw <- do.call(rbind, results)
  probability_success <- mean(raw$success)

  data.frame(
    n_sim = n_sim,
    probability_success = probability_success,
    mcse_success = sqrt(
      probability_success * (1 - probability_success) / n_sim
    ),
    bias_delta = mean(raw$estimate) - true_delta,
    rmse_delta = sqrt(mean((raw$estimate - true_delta)^2)),
    mean_historical_posterior_weight = mean(raw$historical_weight),
    mean_nominal_borrowing_index = mean(raw$borrowing_index),
    mean_sigma_estimate = mean(raw$sigma_estimate)
  )
}

make_commensurability_grid <- function(
    max_ess = 10,
    prior_family = c("gamma", "scaled_beta"),
    gamma_shape = 0.5,
    gamma_rate = 0.1,
    beta_shape1 = 0.2,
    beta_shape2 = 0.2,
    min_ess = 0.01,
    n_grid = 21) {
  prior_family <- match.arg(prior_family)
  stopifnot(max_ess > min_ess, gamma_shape > 0, gamma_rate > 0)
  if (prior_family == "gamma") {
    log_ess <- seq(log(min_ess), log(max_ess), length.out = n_grid)
    ess <- exp(log_ess)
    # Numerical integration is performed on log(ESS); include its Jacobian.
    log_prior_mass <- dgamma(
      ess,
      shape = gamma_shape,
      rate = gamma_rate,
      log = TRUE
    ) + log_ess
  } else {
    stopifnot(beta_shape1 > 0, beta_shape2 > 0)
    lower_u <- min_ess / max_ess
    upper_u <- 1 - lower_u
    logit_u <- seq(qlogis(lower_u), qlogis(upper_u), length.out = n_grid)
    u <- plogis(logit_u)
    ess <- max_ess * u
    # Integrate on logit(u); du/dlogit(u) = u(1-u).
    log_prior_mass <- dbeta(
      u,
      shape1 = beta_shape1,
      shape2 = beta_shape2,
      log = TRUE
    ) + log(u) + log1p(-u)
  }
  prior_weight <- exp(log_prior_mass - log_sum_exp(log_prior_mass))

  data.frame(
    ess = ess,
    prior_weight = prior_weight
  )
}

fit_commensurate_bayes <- function(
    data,
    nuisance_prior_type = c("pilot_conservative", "weak", "pilot_strong"),
    max_ess = 10,
    commensurability_prior = c("gamma", "scaled_beta"),
    gamma_shape = 0.5,
    gamma_rate = 0.1,
    beta_shape1 = 0.2,
    beta_shape2 = 0.2,
    min_ess = 0.01,
    n_grid = 21,
    posterior_cutoff = 0.95) {
  nuisance_prior_type <- match.arg(nuisance_prior_type)
  commensurability_prior <- match.arg(commensurability_prior)
  prior <- make_nuisance_prior(nuisance_prior_type)
  grid <- make_commensurability_grid(
    max_ess = max_ess,
    prior_family = commensurability_prior,
    gamma_shape = gamma_shape,
    gamma_rate = gamma_rate,
    beta_shape1 = beta_shape1,
    beta_shape2 = beta_shape2,
    min_ess = min_ess,
    n_grid = n_grid
  )
  x <- model.matrix(~ treatment + baseline_c + kenya, data = data)
  y <- data$log_sauc_centered

  components <- lapply(grid$ess, function(ess) {
    fit_nig_component(
      x, y, prior,
      intercept_mean = 0,
      intercept_sd = prior$variance_scale / sqrt(ess)
    )
  })
  log_weights <- log(grid$prior_weight) + vapply(
    components, `[[`, numeric(1), "log_marginal"
  )
  posterior_weights <- exp(log_weights - log_sum_exp(log_weights))
  treatment_means <- vapply(
    components, `[[`, numeric(1), "treatment_mean"
  )
  probability <- sum(posterior_weights * vapply(
    components, `[[`, numeric(1), "prob_benefit"
  ))

  list(
    treatment_mean = sum(posterior_weights * treatment_means),
    prob_benefit = probability,
    success = probability >= posterior_cutoff,
    posterior_mean_ess = sum(posterior_weights * grid$ess),
    posterior_prob_ess_below_1 = sum(posterior_weights[grid$ess < 1]),
    posterior_prob_ess_above_half_max = sum(
      posterior_weights[grid$ess > max_ess / 2]
    ),
    posterior_sigma_mean = sum(posterior_weights * vapply(
      components, `[[`, numeric(1), "posterior_sigma_mean"
    ))
  )
}

simulate_commensurate_oc <- function(
    n_sim = 2000,
    seed = 20260717,
    n_per_arm = 45,
    reduction = 0.20,
    residual_sd_multiplier = 1,
    baseline_slope_multiplier = 1,
    current_control_drift = 0,
    nuisance_prior_type = "pilot_conservative",
    max_ess = 10,
    commensurability_prior = c("gamma", "scaled_beta"),
    gamma_shape = 0.5,
    gamma_rate = 0.1,
    beta_shape1 = 0.2,
    beta_shape2 = 0.2,
    min_ess = 0.01,
    n_grid = 21,
    posterior_cutoff = 0.95,
    ...) {
  commensurability_prior <- match.arg(commensurability_prior)
  set.seed(seed)
  true_delta <- log(1 - reduction)
  results <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    dat <- simulate_nuisance_trial(
      n_per_arm = n_per_arm,
      reduction = reduction,
      residual_sd_multiplier = residual_sd_multiplier,
      baseline_slope_multiplier = baseline_slope_multiplier,
      current_control_drift = current_control_drift,
      ...
    )
    fit <- fit_commensurate_bayes(
      dat,
      nuisance_prior_type = nuisance_prior_type,
      max_ess = max_ess,
      commensurability_prior = commensurability_prior,
      gamma_shape = gamma_shape,
      gamma_rate = gamma_rate,
      beta_shape1 = beta_shape1,
      beta_shape2 = beta_shape2,
      min_ess = min_ess,
      n_grid = n_grid,
      posterior_cutoff = posterior_cutoff
    )
    results[[i]] <- data.frame(
      success = fit$success,
      estimate = fit$treatment_mean,
      posterior_mean_ess = fit$posterior_mean_ess,
      prob_ess_below_1 = fit$posterior_prob_ess_below_1,
      prob_ess_above_half_max = fit$posterior_prob_ess_above_half_max,
      sigma_estimate = fit$posterior_sigma_mean
    )
  }
  raw <- do.call(rbind, results)
  probability_success <- mean(raw$success)

  data.frame(
    n_sim = n_sim,
    probability_success = probability_success,
    mcse_success = sqrt(
      probability_success * (1 - probability_success) / n_sim
    ),
    bias_delta = mean(raw$estimate) - true_delta,
    rmse_delta = sqrt(mean((raw$estimate - true_delta)^2)),
    mean_posterior_ess = mean(raw$posterior_mean_ess),
    mean_prob_ess_below_1 = mean(raw$prob_ess_below_1),
    mean_prob_ess_above_half_max = mean(raw$prob_ess_above_half_max),
    mean_sigma_estimate = mean(raw$sigma_estimate)
  )
}

fit_nuisance_bayes <- function(
    data,
    prior_type = c("weak", "pilot_conservative", "pilot_strong"),
    posterior_cutoff = 0.95) {
  prior_type <- match.arg(prior_type)
  prior <- make_nuisance_prior(prior_type)

  x <- model.matrix(~ treatment + baseline_c + kenya, data = data)
  y <- data$log_sauc_centered
  n <- nrow(x)

  m0 <- prior$coefficient_mean
  v0 <- prior$coefficient_covariance_scale
  p0 <- solve(v0)
  vn <- solve(p0 + crossprod(x))
  mn <- vn %*% (p0 %*% m0 + crossprod(x, y))

  nu_n <- prior$variance_df + n
  sum_squares <-
    prior$variance_df * prior$variance_scale^2 +
    sum(y^2) +
    drop(crossprod(m0, p0 %*% m0)) -
    drop(crossprod(mn, solve(vn, mn)))
  sigma2_scale_n <- sum_squares / nu_n

  treatment_index <- match("treatment", colnames(x))
  treatment_mean <- drop(mn[treatment_index])
  treatment_scale <- sqrt(
    sigma2_scale_n * vn[treatment_index, treatment_index]
  )
  prob_benefit <- pt(
    (0 - treatment_mean) / treatment_scale,
    df = nu_n
  )
  interval <- treatment_mean +
    qt(c(0.025, 0.975), df = nu_n) * treatment_scale

  baseline_index <- match("baseline_c", colnames(x))
  baseline_mean <- drop(mn[baseline_index])
  posterior_sigma_mean <- sqrt(
    sum_squares / max(nu_n - 2, 1)
  )

  list(
    treatment_mean = treatment_mean,
    treatment_scale = treatment_scale,
    treatment_df = nu_n,
    prob_benefit = prob_benefit,
    success = prob_benefit >= posterior_cutoff,
    ci_lower = interval[1],
    ci_upper = interval[2],
    baseline_slope_mean = baseline_mean,
    posterior_sigma_mean = posterior_sigma_mean
  )
}

simulate_nuisance_oc <- function(
    n_sim = 5000,
    seed = 20260717,
    n_per_arm = 45,
    reduction = 0.20,
    residual_sd_multiplier = 1,
    baseline_slope_multiplier = 1,
    prior_type = c("weak", "pilot_conservative", "pilot_strong"),
    posterior_cutoff = 0.95,
    ...) {
  prior_type <- match.arg(prior_type)
  set.seed(seed)
  true_delta <- log(1 - reduction)
  true_sigma <- pilot_nuisance$residual_sd * residual_sd_multiplier
  true_slope <- pilot_nuisance$baseline_slope * baseline_slope_multiplier

  results <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    dat <- simulate_nuisance_trial(
      n_per_arm = n_per_arm,
      reduction = reduction,
      residual_sd_multiplier = residual_sd_multiplier,
      baseline_slope_multiplier = baseline_slope_multiplier,
      ...
    )
    fit <- fit_nuisance_bayes(
      dat,
      prior_type = prior_type,
      posterior_cutoff = posterior_cutoff
    )
    results[[i]] <- data.frame(
      success = fit$success,
      estimate = fit$treatment_mean,
      covered = fit$ci_lower <= true_delta && fit$ci_upper >= true_delta,
      sigma_estimate = fit$posterior_sigma_mean,
      slope_estimate = fit$baseline_slope_mean
    )
  }
  raw <- do.call(rbind, results)
  probability_success <- mean(raw$success)

  data.frame(
    n_sim = n_sim,
    probability_success = probability_success,
    mcse_success = sqrt(
      probability_success * (1 - probability_success) / n_sim
    ),
    bias_delta = mean(raw$estimate) - true_delta,
    rmse_delta = sqrt(mean((raw$estimate - true_delta)^2)),
    coverage_95 = mean(raw$covered),
    mean_sigma_estimate = mean(raw$sigma_estimate),
    sigma_bias = mean(raw$sigma_estimate) - true_sigma,
    mean_slope_estimate = mean(raw$slope_estimate),
    slope_bias = mean(raw$slope_estimate) - true_slope
  )
}

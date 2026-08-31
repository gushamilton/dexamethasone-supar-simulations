# Fast Bayesian operating-characteristic engine for log standardised AUC.
#
# This is the first-stage design model. It uses the pilot-calibrated joint
# distribution of baseline log suPAR and log sAUC. It deliberately applies the
# treatment effect only to the post-randomisation outcome, never to baseline.
# The later longitudinal engine will generate irregular Day 0/1/3(/5) samples.

pilot_summary <- list(
  n_historical = 29,
  mean_log_baseline = 4.6967,
  sd_log_baseline = 0.4550,
  mean_log_sauc = 4.385132,
  sd_log_sauc = 0.5944852,
  rho_baseline_sauc = 0.7473633
)

pilot_summary$baseline_slope <- with(
  pilot_summary,
  rho_baseline_sauc * sd_log_sauc / sd_log_baseline
)
pilot_summary$residual_sd <- with(
  pilot_summary,
  sd_log_sauc * sqrt(1 - rho_baseline_sauc^2)
)

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

make_prior_components <- function(
    method = c("none", "fixed", "robust"),
    historical_ess = 10,
    robust_weight = 0.50,
    treatment_prior_sd = 1.00,
    historical_intercept = pilot_summary$mean_log_sauc,
    sigma = pilot_summary$residual_sd) {
  method <- match.arg(method)

  # Parameter order: UK control intercept, treatment, centred baseline, Kenya.
  common_mean <- c(
    historical_intercept,
    0,
    0,
    0
  )

  vague_sd <- c(1.00, treatment_prior_sd, 1.00, 0.75)
  vague <- list(
    weight = 1,
    mean = common_mean,
    covariance = diag(vague_sd^2),
    label = "weak"
  )

  if (method == "none" || historical_ess <= 0) {
    return(list(vague))
  }

  historical_sd <- vague_sd
  historical_sd[1] <- sigma / sqrt(historical_ess)
  historical <- list(
    weight = 1,
    mean = common_mean,
    covariance = diag(historical_sd^2),
    label = paste0("historical_ess_", historical_ess)
  )

  if (method == "fixed") {
    return(list(historical))
  }

  historical$weight <- robust_weight
  vague$weight <- 1 - robust_weight
  list(historical, vague)
}

simulate_summary_trial <- function(
    n_per_arm = 45,
    reduction = 0.25,
    allocation_kenya = 0.50,
    current_control_drift = 0,
    kenya_control_shift = 0,
    kenya_treatment_log_interaction = 0,
    residual_sd_multiplier = 1,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_total <- 2 * n_per_arm
  treatment <- rep(c(0, 1), each = n_per_arm)

  # Keep country approximately balanced within treatment for this design-stage
  # engine, reflecting stratified randomisation.
  n_ke_arm <- round(n_per_arm * allocation_kenya)
  kenya <- unlist(lapply(c(0, 1), function(z) {
    sample(c(rep(1, n_ke_arm), rep(0, n_per_arm - n_ke_arm)))
  }))

  baseline <- rnorm(
    n_total,
    mean = pilot_summary$mean_log_baseline,
    sd = pilot_summary$sd_log_baseline
  )
  baseline_c <- baseline - pilot_summary$mean_log_baseline

  delta <- log(1 - reduction)
  residual_sd <- pilot_summary$residual_sd * residual_sd_multiplier

  log_sauc <-
    pilot_summary$mean_log_sauc +
    current_control_drift +
    pilot_summary$baseline_slope * baseline_c +
    kenya_control_shift * kenya +
    delta * treatment +
    kenya_treatment_log_interaction * kenya * treatment +
    rnorm(n_total, 0, residual_sd)

  data.frame(
    log_sauc = log_sauc,
    treatment = treatment,
    baseline_c = baseline_c,
    kenya = kenya
  )
}

posterior_component <- function(y, X, sigma, prior) {
  v0 <- prior$covariance
  p0 <- solve(v0)
  vn <- solve(p0 + crossprod(X) / sigma^2)
  mn <- vn %*% (p0 %*% prior$mean + crossprod(X, y) / sigma^2)

  # Marginal log density with constants common to all mixture components
  # omitted. This is sufficient for posterior mixture weights.
  quad <-
    sum(y^2) / sigma^2 +
    drop(crossprod(prior$mean, p0 %*% prior$mean)) -
    drop(crossprod(mn, solve(vn, mn)))
  determinant_term <-
    as.numeric(determinant(vn, logarithm = TRUE)$modulus) -
    as.numeric(determinant(v0, logarithm = TRUE)$modulus)

  list(
    mean = drop(mn),
    covariance = vn,
    log_weight_unnormalized = log(prior$weight) + 0.5 * determinant_term - 0.5 * quad,
    label = prior$label
  )
}

fit_summary_bayes <- function(
    data,
    prior_method = c("none", "fixed", "robust"),
    historical_ess = 10,
    robust_weight = 0.50,
    treatment_prior_sd = 1.00,
    historical_intercept = pilot_summary$mean_log_sauc,
    posterior_cutoff = 0.95,
    clinically_relevant_ratio = 0.80,
    analysis_residual_sd = pilot_summary$residual_sd) {
  prior_method <- match.arg(prior_method)
  X <- model.matrix(
    ~ treatment + baseline_c + kenya,
    data = data
  )
  y <- data$log_sauc

  priors <- make_prior_components(
    method = prior_method,
    historical_ess = historical_ess,
    robust_weight = robust_weight,
    treatment_prior_sd = treatment_prior_sd,
    historical_intercept = historical_intercept,
    sigma = analysis_residual_sd
  )
  components <- lapply(
    priors,
    function(prior) posterior_component(
      y = y,
      X = X,
      sigma = analysis_residual_sd,
      prior = prior
    )
  )

  log_weights <- vapply(
    components,
    `[[`,
    numeric(1),
    "log_weight_unnormalized"
  )
  weights <- exp(log_weights - log_sum_exp(log_weights))

  treatment_index <- match("treatment", colnames(X))
  component_means <- vapply(
    components,
    function(z) z$mean[treatment_index],
    numeric(1)
  )
  component_sds <- vapply(
    components,
    function(z) sqrt(z$covariance[treatment_index, treatment_index]),
    numeric(1)
  )

  prob_benefit <- sum(
    weights * pnorm(0, mean = component_means, sd = component_sds)
  )
  clinical_threshold <- log(clinically_relevant_ratio)
  prob_clinically_relevant <- sum(
    weights * pnorm(
      clinical_threshold,
      mean = component_means,
      sd = component_sds
    )
  )

  posterior_mean <- sum(weights * component_means)
  posterior_variance <- sum(
    weights * (component_sds^2 + component_means^2)
  ) - posterior_mean^2

  historical_weight <- sum(
    weights[vapply(components, function(z) z$label != "weak", logical(1))]
  )

  list(
    posterior_mean = posterior_mean,
    posterior_sd = sqrt(posterior_variance),
    prob_benefit = prob_benefit,
    prob_clinically_relevant = prob_clinically_relevant,
    success = prob_benefit >= posterior_cutoff,
    historical_posterior_weight = historical_weight,
    nominal_borrowed_ess = historical_weight * historical_ess
  )
}

simulate_operating_characteristics <- function(
    n_sim = 1000,
    seed = 20260717,
    ...) {
  set.seed(seed)
  dots <- list(...)
  true_delta <- log(1 - dots$reduction)

  results <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    dat <- do.call(
      simulate_summary_trial,
      dots[intersect(names(dots), names(formals(simulate_summary_trial)))]
    )
    fit <- do.call(
      fit_summary_bayes,
      c(
        list(data = dat),
        dots[intersect(names(dots), names(formals(fit_summary_bayes)))]
      )
    )
    results[[i]] <- data.frame(
      posterior_mean = fit$posterior_mean,
      posterior_sd = fit$posterior_sd,
      prob_benefit = fit$prob_benefit,
      prob_clinically_relevant = fit$prob_clinically_relevant,
      success = fit$success,
      historical_posterior_weight = fit$historical_posterior_weight,
      nominal_borrowed_ess = fit$nominal_borrowed_ess
    )
  }

  raw <- do.call(rbind, results)
  data.frame(
    n_sim = n_sim,
    probability_success = mean(raw$success),
    mcse_success = sqrt(mean(raw$success) * (1 - mean(raw$success)) / n_sim),
    mean_posterior_delta = mean(raw$posterior_mean),
    bias_delta = mean(raw$posterior_mean) - true_delta,
    rmse_delta = sqrt(mean((raw$posterior_mean - true_delta)^2)),
    mean_posterior_sd = mean(raw$posterior_sd),
    mean_prob_clinically_relevant = mean(raw$prob_clinically_relevant),
    mean_historical_posterior_weight = mean(raw$historical_posterior_weight),
    mean_nominal_borrowed_ess = mean(raw$nominal_borrowed_ess)
  )
}

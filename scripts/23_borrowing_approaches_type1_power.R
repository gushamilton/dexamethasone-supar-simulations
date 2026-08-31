library(tidyverse)
library(patchwork)

source("R/day5_auc_engine.R")
source("R/nuisance_borrowing_engine.R")

# Paired complete-data simulations using the same Day 0-5 Gaussian trajectory
# and common post-baseline treatment effect as Panel A. The covariance is fixed
# at the pilot/projected design value. Only historical control-trajectory mean
# borrowing varies between approaches. Borrowing uses observed pilot Days 0,
# 1 and 3 only; the projected Day-5 mean retains a weak prior.
n_sim <- as.integer(Sys.getenv("N_SIM", unset = "20000"))
n_per_arm <- 40
historical_n <- pilot_day5$observed_pilot_n
treatment_prior_sd <- 1.00
posterior_cutoff <- 0.95
decision_rule <- "P(Delta < 0 | data) >= 0.95"
drift_grid <- c(-20, -10, 0, 10, 20)
reduction_grid <- c(0, 0.20)

control_design <- cbind(diag(4), 0)
treatment_design <- control_design
treatment_design[2:4, 5] <- 1
design <- rbind(control_design, treatment_design)
observation_covariance <- as.matrix(Matrix::bdiag(
  pilot_day5$covariance / n_per_arm,
  pilot_day5$covariance / n_per_arm
))
observation_precision <- solve(observation_covariance)
prior_mean <- c(pilot_day5$mean, 0)

make_component <- function(ess = 0, label = "weak") {
  control_covariance <- if (ess <= 0) {
    diag(rep(2.00^2, 4))
  } else {
    as.matrix(Matrix::bdiag(
      pilot_day5$covariance[1:3, 1:3, drop = FALSE] / ess,
      2.00^2
    ))
  }
  prior_covariance <- as.matrix(Matrix::bdiag(
    control_covariance,
    treatment_prior_sd^2
  ))
  prior_precision <- solve(prior_covariance)
  posterior_covariance <- solve(
    prior_precision + crossprod(design, observation_precision %*% design)
  )
  marginal_covariance <- observation_covariance +
    design %*% prior_covariance %*% t(design)
  list(
    ess = ess,
    label = label,
    prior_precision = prior_precision,
    posterior_covariance = posterior_covariance,
    marginal_mean = drop(design %*% prior_mean),
    marginal_precision = solve(marginal_covariance),
    marginal_log_determinant = as.numeric(
      determinant(marginal_covariance, logarithm = TRUE)$modulus
    )
  )
}

weak_component <- make_component(label = "weak")
historical_component_5 <- make_component(ess = 5, label = "historical")

switch_grid <- make_commensurability_grid(
  max_ess = 5,
  prior_family = "scaled_beta",
  beta_shape1 = 0.10,
  beta_shape2 = 0.50,
  min_ess = 0.01,
  n_grid = 31
)

approaches <- list(
  "No historical mean borrowing (primary)" = list(
    components = list(weak_component),
    weights = 1
  ),
  "Robust mixture (10%; ESS cap 5)" = list(
    components = list(historical_component_5, weak_component),
    weights = c(0.10, 0.90)
  ),
  "Switch-commensurate (ESS cap 5)" = list(
    components = lapply(
      switch_grid$ess,
      function(value) make_component(ess = value, label = "historical")
    ),
    weights = switch_grid$prior_weight
  ),
  "Power prior a0 = 0.2" = list(
    components = list(make_component(0.2 * historical_n, "historical")),
    weights = 1
  ),
  "Power prior a0 = 0.4" = list(
    components = list(make_component(0.4 * historical_n, "historical")),
    weights = 1
  ),
  "Power prior a0 = 0.6" = list(
    components = list(make_component(0.6 * historical_n, "historical")),
    weights = 1
  )
)

row_log_density <- function(values, component) {
  residual <- sweep(values, 2, component$marginal_mean, `-`)
  quadratic <- rowSums(
    (residual %*% component$marginal_precision) * residual
  )
  -0.5 * (
    ncol(values) * log(2 * pi) +
      component$marginal_log_determinant + quadratic
  )
}

component_posterior <- function(values, component) {
  score <- values %*% observation_precision %*% design
  prior_score <- drop(component$prior_precision %*% prior_mean)
  posterior_means <- (
    score + matrix(prior_score, nrow(values), 5, byrow = TRUE)
  ) %*% component$posterior_covariance
  treatment_sd <- sqrt(component$posterior_covariance[5, 5])
  list(
    probability_benefit = pnorm(
      0,
      mean = posterior_means[, 5],
      sd = treatment_sd
    ),
    log_marginal = row_log_density(values, component)
  )
}

fit_approach <- function(values, approach) {
  component_fits <- lapply(
    approach$components,
    function(component) component_posterior(values, component)
  )
  if (length(component_fits) == 1) {
    probability_benefit <- component_fits[[1]]$probability_benefit
  } else {
    log_weights <- vapply(
      seq_along(component_fits),
      function(index) {
        log(approach$weights[index]) + component_fits[[index]]$log_marginal
      },
      numeric(nrow(values))
    )
    row_maximum <- apply(log_weights, 1, max)
    posterior_weights <- exp(log_weights - row_maximum)
    posterior_weights <- posterior_weights / rowSums(posterior_weights)
    component_probabilities <- vapply(
      component_fits,
      `[[`,
      numeric(nrow(values)),
      "probability_benefit"
    )
    probability_benefit <- rowSums(
      posterior_weights * component_probabilities
    )
  }
  success <- probability_benefit >= posterior_cutoff
  probability_success <- mean(success)
  tibble(
    probability_success = probability_success,
    mcse_success = sqrt(
      probability_success * (1 - probability_success) / n_sim
    )
  )
}

simulate_scenario <- function(drift_pct, reduction, seed) {
  set.seed(seed)
  # Baseline remains exchangeable; drift is applied to the post-baseline
  # control trajectory, where historical-current non-exchangeability matters
  # for the treatment contrast.
  current_control_mean <- pilot_day5$mean +
    c(0, rep(log1p(drift_pct / 100), 3))
  post_reduction <- calibrate_day5_effect(
    reduction,
    shape = "immediate",
    control_log_mean = current_control_mean
  )
  treatment_shift <- effect_vector_day5(post_reduction, "immediate")
  control_means <- MASS::mvrnorm(
    n_sim,
    mu = current_control_mean,
    Sigma = pilot_day5$covariance / n_per_arm
  )
  treatment_means <- MASS::mvrnorm(
    n_sim,
    mu = current_control_mean + treatment_shift,
    Sigma = pilot_day5$covariance / n_per_arm
  )
  values <- cbind(control_means, treatment_means)

  bind_rows(lapply(names(approaches), function(approach_name) {
    fit_approach(values, approaches[[approach_name]]) |>
      mutate(approach = approach_name)
  })) |>
    mutate(
      drift_pct = drift_pct,
      reduction = reduction,
      n_sim = n_sim,
      n_per_arm = n_per_arm,
      effect_shape = "immediate sustained",
      day3_day5_correlation = pilot_day5$projected_day3_day5_correlation,
      treatment_prior_sd = treatment_prior_sd,
      posterior_cutoff = posterior_cutoff,
      decision_rule = decision_rule
    )
}

scenario_grid <- crossing(
  drift_pct = drift_grid,
  reduction = reduction_grid
)

results <- bind_rows(parallel::mclapply(
  seq_len(nrow(scenario_grid)),
  function(index) {
    row <- scenario_grid[index, ]
    simulate_scenario(
      row$drift_pct,
      row$reduction,
      seed = 240000 +
        1000 * match(row$drift_pct, drift_grid) +
        round(100 * row$reduction)
    )
  },
  mc.cores = min(4L, parallel::detectCores())
))

write_csv(results, "outputs/borrowing_approaches_type1_power.csv")

approach_levels <- names(approaches)
approach_colours <- c(
  "No historical mean borrowing (primary)" = "#222222",
  "Robust mixture (10%; ESS cap 5)" = "#2E74B5",
  "Switch-commensurate (ESS cap 5)" = "#009E73",
  "Power prior a0 = 0.2" = "#E69F00",
  "Power prior a0 = 0.4" = "#D55E00",
  "Power prior a0 = 0.6" = "#CC79A7"
)
approach_linetypes <- c(
  "No historical mean borrowing (primary)" = "solid",
  "Robust mixture (10%; ESS cap 5)" = "solid",
  "Switch-commensurate (ESS cap 5)" = "solid",
  "Power prior a0 = 0.2" = "dashed",
  "Power prior a0 = 0.4" = "dashed",
  "Power prior a0 = 0.6" = "dashed"
)

plot_data <- results |>
  mutate(approach = factor(approach, levels = approach_levels))

theme_borrowing <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5, colour = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = .35),
      legend.position = "bottom",
      legend.text = element_text(size = 8.3),
      legend.title = element_blank()
    )
}

panel_type1 <- plot_data |>
  filter(reduction == 0) |>
  ggplot(aes(drift_pct, probability_success, colour = approach, linetype = approach)) +
  geom_hline(yintercept = .05, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = .95) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = approach_colours) +
  scale_linetype_manual(values = approach_linetypes) +
  scale_x_continuous(breaks = drift_grid) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "A. Type-I error",
    subtitle = "No treatment effect",
    x = "Current post-baseline control trajectory vs historical (%)",
    y = "Type-I error"
  ) +
  theme_borrowing() +
  theme(legend.position = "none")

panel_power <- plot_data |>
  filter(reduction == .20) |>
  ggplot(aes(drift_pct, probability_success, colour = approach, linetype = approach)) +
  geom_hline(yintercept = .80, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = .95) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = approach_colours) +
  scale_linetype_manual(values = approach_linetypes) +
  scale_x_continuous(breaks = drift_grid) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "B. Power",
    subtitle = "True 20% reduction in AUC",
    x = "Current post-baseline control trajectory vs historical (%)",
    y = "Bayesian power"
  ) +
  theme_borrowing()

figure <- panel_type1 | panel_power
figure <- figure + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Historical control-trajectory mean borrowing: Type-I error and power",
    subtitle = paste0(
      "Day 0-5 model-based AUC; conservative contributing N = 80; ",
      "decision rule P(treatment reduces AUC | data) >= 0.95; ",
      "common treatment-effect prior N(0, 1^2)"
    ),
    caption = paste0(
      scales::comma(n_sim),
      " paired simulations/scenario. Fixed power priors use pilot n = ",
      historical_n,
      " (nominal ESS 5.8, 11.6 and 17.4 for observed pilot Days 0, 1 and 3; Day 5 is not borrowed). Negative values mean the current post-baseline control trajectory is below the historical cohort; baseline is matched."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.5, colour = "grey35")
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  "figures/borrowing_approaches_type1_power.png",
  figure,
  width = 14,
  height = 6.8,
  dpi = 300
)
ggsave(
  "figures/borrowing_approaches_type1_power.pdf",
  figure,
  width = 14,
  height = 6.8
)

print(results, n = Inf)

library(tidyverse)
library(patchwork)

source("R/day5_auc_engine.R")

n_sim_panel_a <- 10000
n_sim_panel_b <- 5000
posterior_cutoff <- 0.95
treatment_prior_sd <- 0.50
central_auc_reduction <- 0.20

run_rows <- function(grid, fun) {
  rows <- split(grid, seq_len(nrow(grid)))
  bind_rows(parallel::mclapply(
    seq_along(rows),
    function(index) bind_cols(rows[[index]], fun(rows[[index]], index)),
    mc.cores = min(4L, parallel::detectCores())
  ))
}

simulate_complete_day5_auc <- function(
    n_sim,
    seed,
    n_per_arm,
    target_auc_reduction) {
  set.seed(seed)
  covariance <- pilot_day5$covariance
  inverse_covariance <- solve(covariance)
  control_design <- cbind(diag(4), matrix(0, 4, 3))
  treatment_design <- control_design
  treatment_design[2, 5] <- 1
  treatment_design[3, 6] <- 1
  treatment_design[4, 7] <- 1

  prior <- make_day5_prior(treatment_prior_sd)
  prior_precision <- solve(prior$covariance)
  information <- n_per_arm * (
    crossprod(control_design, inverse_covariance %*% control_design) +
      crossprod(treatment_design, inverse_covariance %*% treatment_design)
  )
  posterior_covariance <- solve(prior_precision + information)
  post_reduction <- calibrate_day5_effect(
    target_auc_reduction,
    shape = "immediate"
  )
  treatment_shift <- effect_vector_day5(post_reduction, "immediate")

  success <- logical(n_sim)
  estimate <- numeric(n_sim)
  for (simulation in seq_len(n_sim)) {
    control_mean <- MASS::mvrnorm(
      1,
      pilot_day5$mean,
      covariance / n_per_arm
    )
    treatment_mean <- MASS::mvrnorm(
      1,
      pilot_day5$mean + treatment_shift,
      covariance / n_per_arm
    )
    score <- n_per_arm * (
      crossprod(
        control_design,
        inverse_covariance %*% control_mean
      )[, 1] +
        crossprod(
          treatment_design,
          inverse_covariance %*% treatment_mean
        )[, 1]
    )
    posterior_mean <- posterior_covariance %*% (
      prior_precision %*% prior$mean + score
    )
    estimate[simulation] <- log_auc_delta_day5(posterior_mean)
    gradient <- numeric_gradient_day5(
      log_auc_delta_day5,
      posterior_mean
    )
    posterior_sd <- sqrt(drop(crossprod(
      gradient,
      posterior_covariance %*% gradient
    )))
    success[simulation] <-
      pnorm(0, estimate[simulation], posterior_sd) >= posterior_cutoff
  }

  true_delta <- log(1 - target_auc_reduction)
  probability_success <- mean(success)
  tibble(
    probability_success = probability_success,
    mcse_success = sqrt(
      probability_success * (1 - probability_success) / n_sim
    ),
    bias_delta = mean(estimate) - true_delta,
    n_sim = n_sim
  )
}

grid_a <- crossing(
  total_n = seq(60, 100, by = 10),
  reduction = c(0, .15, .20, .25, .30)
)

result_a <- run_rows(grid_a, function(row, index) {
  simulate_complete_day5_auc(
    n_sim = n_sim_panel_a,
    seed = 200000 + 100 * row$total_n + round(100 * row$reduction),
    n_per_arm = row$total_n / 2,
    target_auc_reduction = row$reduction
  )
})

simulate_both_day5_endpoints <- function(
    n_sim,
    seed,
    missing_mechanism,
    day3_missing_rate) {
  set.seed(seed)
  auc_success <- day3_success <- logical(n_sim)
  achieved_missing <- numeric(n_sim)

  for (simulation in seq_len(n_sim)) {
    trial <- simulate_day5_trial(
      n_per_arm = 45,
      target_auc_reduction = central_auc_reduction,
      effect_shape = "immediate",
      missing_mechanism = missing_mechanism,
      day3_missing_rate = day3_missing_rate
    )
    auc_fit <- fit_day5_auc_bayes(
      trial,
      treatment_prior_sd = treatment_prior_sd,
      posterior_cutoff = posterior_cutoff
    )
    day3_fit <- fit_day3_ancova_bayes(
      trial,
      treatment_prior_sd = treatment_prior_sd,
      posterior_cutoff = posterior_cutoff
    )
    auc_success[simulation] <- auc_fit$success
    day3_success[simulation] <- day3_fit$success
    achieved_missing[simulation] <- mean(
      !trial$observed[trial$time_index == 3]
    )
  }

  probabilities <- c(mean(auc_success), mean(day3_success))
  tibble(
    analysis = c("Model-based AUC", "Day-3 ANCOVA"),
    probability_success = probabilities,
    mcse_success = sqrt(probabilities * (1 - probabilities) / n_sim),
    achieved_day3_missingness = mean(achieved_missing),
    n_sim = n_sim
  )
}

grid_b <- crossing(
  missing_mechanism = c("mcar", "mar_improvement", "mar_deterioration"),
  day3_missing_rate = seq(0, .50, by = .10)
)

missingness_output <- "outputs/stage2_day5_panel_b_missingness.csv"
if (
  file.exists(missingness_output) &&
    Sys.getenv("RERUN_DAY5_MISSINGNESS", unset = "0") != "1"
) {
  result_b <- read_csv(missingness_output, show_col_types = FALSE)
} else {
  result_b <- run_rows(grid_b, function(row, index) {
    simulate_both_day5_endpoints(
      n_sim = n_sim_panel_b,
      seed = 210000 +
        1000 * match(
          row$missing_mechanism,
          c("mcar", "mar_improvement", "mar_deterioration")
        ) +
        round(100 * row$day3_missing_rate),
      missing_mechanism = row$missing_mechanism,
      day3_missing_rate = row$day3_missing_rate
    )
  })
}

write_csv(result_a, "outputs/stage2_day5_panel_a_power.csv")
write_csv(result_b, missingness_output)

missing_labels <- c(
  mcar = "MCAR",
  mar_improvement = "MAR: improving patients",
  mar_deterioration = "MAR: deteriorating patients"
)

panel_a <- result_a |>
  filter(reduction > 0) |>
  mutate(effect = factor(
    scales::percent(reduction, accuracy = 1),
    levels = scales::percent(c(.15, .20, .25, .30), accuracy = 1)
  )) |>
  ggplot(aes(total_n, probability_success, colour = effect)) +
  geom_hline(yintercept = .80, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = .9) +
  geom_point(size = 2.1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(60, 100, 10)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "A. Power under different treatment effects",
    subtitle = "Day 0-5 model-based AUC; concurrent randomised controls",
    x = "Total randomised N",
    y = "Bayesian power",
    colour = "True AUC reduction"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

panel_b <- result_b |>
  mutate(
    mechanism = factor(
      missing_labels[missing_mechanism],
      levels = unname(missing_labels)
    ),
    analysis = factor(
      analysis,
      levels = c("Model-based AUC", "Day-3 ANCOVA")
    )
  ) |>
  ggplot(aes(
    100 * achieved_day3_missingness,
    probability_success,
    colour = mechanism,
    linetype = analysis
  )) +
  geom_hline(yintercept = .80, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = .9) +
  geom_point(size = 1.8) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 50, 10), limits = c(0, 52)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "B. Robustness when Day-3 samples are missing",
    subtitle = "20% AUC reduction; total N = 90",
    x = "Achieved Day-3 missingness (%)",
    y = "Bayesian power",
    colour = "Missingness mechanism",
    linetype = "Analysis"
  ) +
  theme_bw(base_size = 11) +
  guides(
    linetype = guide_legend(order = 1, nrow = 1),
    colour = guide_legend(order = 2, nrow = 2)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 8.5)
  )

figure <- panel_a | panel_b
figure <- figure + plot_annotation(
  title = "Dexamethasone-suPAR trial: primary operating characteristics",
  subtitle = paste0(
    "Decision rule: P(treatment reduces AUC | data) >= 0.95; ",
    "treatment prior N(0, 0.5^2); no historical control-mean borrowing"
  ),
  caption = paste0(
    "Panel A: 10,000 simulations/cell. Panel B: 5,000 paired ",
    "simulations/cell.\nFocused loss of the Day-3 sample with Days 0, 1 ",
    "and 5 retained. Day-5 covariance is a conservative projection from ",
    "the Day 0-3 pilot. Immediate sustained effect. Dotted lines denote ",
    "80% power."
  ),
  theme = theme(
    plot.caption = element_text(size = 8.5, hjust = 0, margin = margin(t = 8))
  )
)

ggsave(
  "figures/stage2_day5_two_panel.png",
  figure,
  width = 14,
  height = 6.8,
  dpi = 300
)
ggsave(
  "figures/stage2_day5_two_panel.pdf",
  figure,
  width = 14,
  height = 6.8
)

print(result_a |>
  filter(total_n == 90) |>
  dplyr::select(total_n, reduction, probability_success, mcse_success))
print(result_b, n = Inf)

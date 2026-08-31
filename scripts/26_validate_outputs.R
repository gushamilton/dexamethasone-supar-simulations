library(readr)

expected_a <- as.integer(Sys.getenv("EXPECTED_A", unset = "20000"))
expected_b <- as.integer(Sys.getenv("EXPECTED_B", unset = "10000"))
expected_cd <- as.integer(Sys.getenv("EXPECTED_CD", unset = "20000"))

panel_a <- read_csv("outputs/stage2_day5_panel_a_power.csv", show_col_types = FALSE)
panel_b <- read_csv("outputs/stage2_day5_panel_b_missingness.csv", show_col_types = FALSE)
borrowing <- read_csv(
  "outputs/borrowing_approaches_type1_power.csv",
  show_col_types = FALSE
)

required <- c(
  "probability_success", "mcse_success", "n_sim",
  "n_per_arm", "effect_shape", "day3_day5_correlation",
  "treatment_prior_sd", "posterior_cutoff", "decision_rule"
)
stopifnot(
  all(required %in% names(panel_a)),
  all(required %in% names(panel_b)),
  all(required %in% names(borrowing)),
  all(panel_a$n_sim == expected_a),
  all(panel_b$n_sim == expected_b),
  all(borrowing$n_sim == expected_cd),
  all(panel_a$treatment_prior_sd == 1),
  all(panel_b$treatment_prior_sd == 1),
  all(borrowing$treatment_prior_sd == 1),
  all(panel_a$posterior_cutoff == .95),
  all(panel_b$posterior_cutoff == .95),
  all(borrowing$posterior_cutoff == .95)
)

expected_rule <- "P(Delta < 0 | data) >= 0.95"
stopifnot(
  all(panel_a$decision_rule == expected_rule),
  all(panel_b$decision_rule == expected_rule),
  all(borrowing$decision_rule == expected_rule),
  all(panel_a$n_per_arm == panel_a$total_n / 2),
  max(panel_a$total_n) == 90,
  all(panel_b$n_per_arm == 40),
  all(borrowing$n_per_arm == 40),
  all(panel_a$effect_shape == "immediate sustained"),
  all(panel_b$effect_shape == "immediate sustained"),
  all(borrowing$effect_shape == "immediate sustained"),
  all(panel_a$day3_day5_correlation == .80),
  all(panel_b$day3_day5_correlation == .80),
  all(borrowing$day3_day5_correlation == .80)
)

for (data in list(panel_a, panel_b, borrowing)) {
  stopifnot(
    all(is.finite(data$probability_success)),
    all(data$probability_success >= 0 & data$probability_success <= 1),
    all(is.finite(data$mcse_success)),
    all(data$mcse_success >= 0)
  )
}

stopifnot(
  file.exists("figures/stage2_final_four_panel.png"),
  file.info("figures/stage2_final_four_panel.png")$size > 10000
)

message("Generated-output validation passed")

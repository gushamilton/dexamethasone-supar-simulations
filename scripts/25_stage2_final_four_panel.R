library(tidyverse)
library(patchwork)

power_n <- read_csv(
  "outputs/stage2_day5_panel_a_power.csv",
  show_col_types = FALSE
)
missingness <- read_csv(
  "outputs/stage2_day5_panel_b_missingness.csv",
  show_col_types = FALSE
)
borrowing <- read_csv(
  "outputs/borrowing_approaches_type1_power.csv",
  show_col_types = FALSE
)

stopifnot(
  all(power_n$treatment_prior_sd == 1),
  all(missingness$treatment_prior_sd == 1),
  all(borrowing$treatment_prior_sd == 1),
  all(power_n$posterior_cutoff == .95),
  all(missingness$posterior_cutoff == .95),
  all(borrowing$posterior_cutoff == .95)
)

theme_stage2 <- function() {
  theme_bw(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.3, colour = "grey30"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = .35),
      legend.text = element_text(size = 7.9),
      legend.title = element_text(size = 8.2),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8.7)
    )
}

effect_colours <- c(
  "15%" = "#F8766D",
  "20%" = "#7CAE00",
  "25%" = "#00BFC4",
  "30%" = "#C77CFF"
)

panel_a <- power_n |>
  filter(reduction > 0) |>
  mutate(effect = factor(
    scales::percent(reduction, accuracy = 1),
    levels = names(effect_colours)
  )) |>
  ggplot(aes(total_n, probability_success, colour = effect)) +
  geom_hline(yintercept = .80, linetype = "dotted", colour = "grey45") +
  geom_vline(xintercept = 80, linetype = "dashed", colour = "grey55") +
  annotate(
    "text", x = 79.2, y = .12,
    label = "Conservative contributing N = 80",
    angle = 90, hjust = 0, size = 2.7, colour = "grey35"
  ) +
  geom_line(linewidth = .9) +
  geom_point(size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = effect_colours) +
  scale_x_continuous(breaks = seq(60, 90, 10)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "A. Power by treatment effect and contributing N",
    subtitle = "Recruitment fixed at N=90; allow five lost to follow-up per arm",
    x = "Participants contributing suPAR data (total N)",
    y = "Bayesian power",
    colour = "True AUC reduction"
  ) +
  theme_stage2() +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")

missing_labels <- c(
  mcar = "MCAR",
  mar_improvement = "MAR: improving patients",
  mar_deterioration = "MAR: deteriorating patients"
)
missing_colours <- c(
  "MCAR" = "#F8766D",
  "MAR: improving patients" = "#00BA38",
  "MAR: deteriorating patients" = "#619CFF"
)

panel_b <- missingness |>
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
  geom_hline(yintercept = .80, linetype = "dotted", colour = "grey45") +
  geom_line(linewidth = .9) +
  geom_point(size = 1.7) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = missing_colours) +
  scale_x_continuous(breaks = seq(0, 50, 10), limits = c(0, 52)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "B. Model-based AUC versus Day-3 ANCOVA",
    subtitle = "20% AUC reduction; conservative contributing N = 80; focused Day-3 loss",
    x = "Achieved Day-3 missingness (%)",
    y = "Bayesian power",
    colour = "Missingness",
    linetype = "Analysis"
  ) +
  theme_stage2() +
  guides(
    linetype = guide_legend(order = 1, nrow = 1),
    colour = guide_legend(order = 2, nrow = 2)
  ) +
  theme(legend.position = "bottom", legend.box = "vertical")

approach_labels <- c(
  "No historical mean borrowing (primary)" = "Primary: no historical mean borrowing",
  "Robust mixture (10%; ESS cap 5)" = "Robust mixture: 10%, ESS <=5",
  "Switch-commensurate (ESS cap 5)" = "Switch-commensurate: ESS <=5",
  "Power prior a0 = 0.2" = "Fixed power prior: a0=0.2",
  "Power prior a0 = 0.4" = "Fixed power prior: a0=0.4",
  "Power prior a0 = 0.6" = "Fixed power prior: a0=0.6"
)
approach_colours <- c(
  "Primary: no historical mean borrowing" = "#222222",
  "Robust mixture: 10%, ESS <=5" = "#2E74B5",
  "Switch-commensurate: ESS <=5" = "#009E73",
  "Fixed power prior: a0=0.2" = "#E69F00",
  "Fixed power prior: a0=0.4" = "#D55E00",
  "Fixed power prior: a0=0.6" = "#CC79A7"
)
approach_linetypes <- c(
  "Primary: no historical mean borrowing" = "solid",
  "Robust mixture: 10%, ESS <=5" = "solid",
  "Switch-commensurate: ESS <=5" = "solid",
  "Fixed power prior: a0=0.2" = "dashed",
  "Fixed power prior: a0=0.4" = "dashed",
  "Fixed power prior: a0=0.6" = "dashed"
)

borrowing_plot <- borrowing |>
  mutate(approach = factor(
    approach_labels[approach],
    levels = unname(approach_labels)
  ))

panel_c <- borrowing_plot |>
  filter(reduction == .20) |>
  ggplot(aes(drift_pct, probability_success, colour = approach, linetype = approach)) +
  geom_hline(yintercept = .80, linetype = "dotted", colour = "grey45") +
  geom_line(linewidth = .9) +
  geom_point(size = 1.9) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = approach_colours) +
  scale_linetype_manual(values = approach_linetypes) +
  scale_x_continuous(breaks = seq(-20, 20, 10)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "C. Power with historical mean-borrowing sensitivities",
    subtitle = "True 20% AUC reduction; pilot Days 0, 1 and 3 only; baseline matched",
    x = "Current post-baseline control trajectory vs historical (%)",
    y = "Bayesian power",
    colour = "Analysis",
    linetype = "Analysis"
  ) +
  theme_stage2()

panel_d <- borrowing_plot |>
  filter(reduction == 0) |>
  ggplot(aes(drift_pct, probability_success, colour = approach, linetype = approach)) +
  geom_hline(yintercept = .05, linetype = "dotted", colour = "grey45") +
  geom_line(linewidth = .9) +
  geom_point(size = 1.9) +
  coord_cartesian(ylim = c(0, .20)) +
  scale_colour_manual(values = approach_colours) +
  scale_linetype_manual(values = approach_linetypes) +
  scale_x_continuous(breaks = seq(-20, 20, 10)) +
  scale_y_continuous(
    breaks = seq(0, .20, .05),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "D. Type-I error under historical-current differences",
    subtitle = "No treatment effect; pilot Days 0, 1 and 3 only; baseline matched",
    x = "Current post-baseline control trajectory vs historical (%)",
    y = "Type-I error",
    colour = "Analysis",
    linetype = "Analysis"
  ) +
  theme_stage2()

top <- panel_a | panel_b
bottom <- (panel_c | panel_d) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

figure <- (top / bottom) +
  plot_layout(heights = c(1, 1.12)) +
  plot_annotation(
    title = "Dexamethasone-suPAR trial: operating characteristics supporting the primary analysis",
    subtitle = paste0(
      "Primary decision: P(Delta < 0 | data) >= 0.95, Delta = log(AUC dexamethasone/AUC control); ",
      "common post-baseline treatment-effect prior N(0, 1^2); recruit N=90, with conservative contributing N=80"
    ),
    caption = stringr::str_wrap(paste0(
      "Panel A: ", scales::comma(unique(power_n$n_sim)),
      " simulations/cell. Panel B: ", scales::comma(unique(missingness$n_sim)),
      " paired simulations/cell. Panels C/D: ",
      scales::comma(unique(borrowing$n_sim)),
      " paired simulations/scenario. Every simulated positive decision applies the stated 0.95 posterior-probability rule. Historical mean borrowing uses observed pilot Days 0, 1 and 3 only; projected Day 5 is not borrowed. Known-covariance Gaussian design approximation using pilot Day 0-3 data and a prespecified Day-5 projection. ",
      "Immediate sustained effect. Dotted lines denote 80% power or 5% Type-I error."
    ), width = 185),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.2, colour = "grey35", hjust = 0)
    )
  )

ggsave(
  "figures/stage2_final_four_panel.png",
  figure,
  width = 14,
  height = 10.4,
  dpi = 300
)
ggsave(
  "figures/stage2_final_four_panel.pdf",
  figure,
  width = 14,
  height = 10.4
)

print(figure)

# Bayesian suPAR AUC simulation study

This repository contains the design-stage simulations for a randomised trial
of dexamethasone in pleural infection.

## Primary design

- Estimand: log ratio of model-based suPAR AUCs from Day 0 to Day 5.
- Scheduled samples: Days 0, 1, 3 and 5; actual elapsed times will be used in
  the trial analysis.
- Decision rule: `P(treatment reduces AUC | data) >= 0.95`.
- Treatment-effect prior: `Normal(0, 1^2)` on one common post-baseline
  log-scale treatment coefficient. This is centred on no treatment effect.
- The primary analysis uses concurrent randomised controls and does not borrow
  the historical control mean.
- Planned recruitment is fixed at 90 participants (45 per arm), allowing for
  five participants lost to follow-up per arm. Every randomised participant
  with at least one suPAR measurement contributes to the longitudinal
  likelihood; operating characteristics nevertheless use a conservative 40
  contributing participants per arm.
- Historical data inform conservative priors for nuisance parameters. Robust
  mean borrowing capped at effective sample size 5 and fixed power priors are
  sensitivity analyses only. These mean-borrowing sensitivities use observed
  pilot Days 0, 1 and 3; projected Day-5 information is not borrowed.

The pilot cohort contains Days 0, 1 and 3. For the Day 0-5 design simulation,
Day 5 is transparently projected to have the Day-3 mean and marginal SD and
correlation 0.80 with Day 3. This is a design assumption, not an observed pilot
estimate.

## Reproduce the main figure

Install R packages `MASS`, `Matrix`, `tidyverse` and `patchwork`, then run from the
repository root:

```sh
Rscript tests/test_day5_auc_engine.R
RERUN_DAY5_MISSINGNESS=1 Rscript scripts/20_stage2_day5_two_panel.R
Rscript scripts/23_borrowing_approaches_type1_power.R
Rscript scripts/24_treatment_prior_calibration.R
Rscript scripts/25_stage2_final_four_panel.R
```

The script creates:

- `figures/stage2_final_four_panel.png` and `.pdf`;
- `outputs/stage2_day5_panel_a_power.csv`;
- `outputs/stage2_day5_panel_b_missingness.csv`;
- `outputs/borrowing_approaches_type1_power.csv`;
- `outputs/treatment_prior_calibration.csv` (supplementary prior calibration).

Panel A uses 20,000 simulations per cell. Panel B uses 10,000 paired simulations
per cell. Panels C and D use 20,000 paired simulations per scenario. Seeds are
fixed. At 80% power the Monte Carlo standard error is approximately 0.3
percentage points with 20,000 simulations and 0.4 points with 10,000.

In every simulated replicate, a positive primary decision is made only when
the fitted posterior probability satisfies `P(Delta < 0 | data) >= 0.95`.
Thus the reported power and Type-I error are operating characteristics of the
prespecified Bayesian decision rule, not of a frequentist significance test.

The simulation is a design-stage Gaussian approximation: it fixes the
within-participant covariance at the pilot/projected value. The final analysis
will estimate covariance parameters using conservative priors. The assumptions
register and outline SAP are in [`docs/simulation_plan.md`](docs/simulation_plan.md) and
[`docs/stage2_primary_analysis_and_sap.md`](docs/stage2_primary_analysis_and_sap.md).

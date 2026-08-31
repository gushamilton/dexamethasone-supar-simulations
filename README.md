# Bayesian suPAR AUC simulation study

This repository contains the design-stage simulations for a randomised trial
of dexamethasone in pleural infection.

## Primary design

- Estimand: log ratio of model-based suPAR AUCs from Day 0 to Day 5.
- Scheduled samples: Days 0, 1, 3 and 5; actual elapsed times will be used in
  the trial analysis.
- Decision rule: `P(treatment reduces AUC | data) >= 0.95`.
- Treatment-effect prior: `Normal(0, 0.5^2)` on each post-baseline log-scale
  treatment coefficient.
- The primary analysis uses concurrent randomised controls and does not borrow
  the historical control mean.
- Historical data inform conservative priors for nuisance parameters. Robust
  mean borrowing capped at effective sample size 5 and fixed power priors are
  sensitivity analyses only.

The pilot cohort contains Days 0, 1 and 3. For the Day 0-5 design simulation,
Day 5 is transparently projected to have the Day-3 mean and marginal SD and
correlation 0.80 with Day 3. This is a design assumption, not an observed pilot
estimate.

## Reproduce the main figure

Install R packages `MASS`, `tidyverse` and `patchwork`, then run from the
repository root:

```sh
Rscript tests/test_day5_auc_engine.R
RERUN_DAY5_MISSINGNESS=1 Rscript scripts/20_stage2_day5_two_panel.R
```

The script creates:

- `figures/stage2_day5_two_panel.png` and `.pdf`;
- `outputs/stage2_day5_panel_a_power.csv`;
- `outputs/stage2_day5_panel_b_missingness.csv`.

Panel A uses 10,000 simulations per cell. Panel B uses 5,000 paired simulations
per cell. Seeds are fixed in the script. At 80% power, the Monte Carlo standard
error is approximately 0.6 percentage points for Panel B.

Earlier scripts in `scripts/` reproduce the historical-mean borrowing,
nuisance-prior, prior-conflict, missingness and UK-Kenya sensitivity analyses.
The rationale and outline SAP are in
[`docs/stage2_primary_analysis_and_sap.md`](docs/stage2_primary_analysis_and_sap.md).

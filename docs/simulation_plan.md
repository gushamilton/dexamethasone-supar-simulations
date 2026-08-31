# Simulation plan and assumptions register

## 1. Primary estimand and decision rule

The primary mechanistic estimand is the ratio of model-based pleural-fluid
suPAR AUCs from Day 0 to Day 5:

\[
R_{AUC}=AUC_{dexamethasone}/AUC_{control}, \qquad
\Delta=\log(R_{AUC}).
\]

Values below 1 favour dexamethasone. The primary success rule is

\[
P(\Delta<0\mid data)\ge 0.95.
\]

The posterior probability of at least a 20% reduction,
`P(Delta <= log(0.80) | data)`, will also be reported. Probability cut-offs of
0.90 and 0.975 are sensitivity analyses.

## 2. Historical information

The historical cohort contains 29 complete Day 0, Day 1 and Day 3
trajectories. On the natural-log scale:

| Quantity | Estimate |
|---|---:|
| Mean Day 0 | 4.6967 |
| Mean Day 1 | 4.3111 |
| Mean Day 3 | 4.1880 |
| SD Day 0 | 0.4550 |
| SD Day 1 | 0.7079 |
| SD Day 3 | 0.7508 |
| Correlation Day 0--1 | 0.6691 |
| Correlation Day 0--3 | 0.5347 |
| Correlation Day 1--3 | 0.9151 |
| SD log individual sAUC 0--3 | 0.5945 |
| Correlation baseline--log sAUC | 0.7474 |

The Stage 2 schedule includes Day 5, but the pilot does not. In the main
design simulation, Day 5 is projected to have the Day-3 mean and marginal SD
and correlation 0.80 with Day 3. This is an explicit design assumption rather
than an observed pilot estimate.

## 3. Prior specification

The historical cohort contains no randomised dexamethasone comparison and does
not provide prior evidence for treatment efficacy.

- One common post-baseline log-scale treatment coefficient has prior
  `Normal(0, 1^2)`. This symmetric prior is centred on no treatment effect.
- The treatment difference is constrained to zero at baseline.
- Concurrent randomised controls estimate the control trajectory. The primary
  analysis does not borrow the historical control mean.
- Standard deviations have pilot-centred log-Normal priors with log-scale SD
  0.5.
- Correlation and baseline-association parameters have pilot-centred Normal
  priors on Fisher-z or otherwise unconstrained scales with SD 0.5.
- Centred control-trajectory coefficients have weak `Normal(0, 2^2)` priors.

Historical control-mean borrowing is a sensitivity analysis only. The selected
robust mixture has 90% weak and 10% pilot-informed components, with the
informative component capped at effective sample size 5. Fixed power-prior
weights `a0 = 0.2, 0.4, 0.6` provide further sensitivities. Borrowing applies
only to observed pilot Days 0, 1 and 3; the projected Day-5 mean retains the
weak primary prior.

## 4. Longitudinal analysis model

The intended primary analysis is a Bayesian log-normal longitudinal
mixed-effects model using actual elapsed sampling time. The mean trajectory is
represented flexibly at Days 0, 1, 3 and 5, with one common post-baseline
treatment effect constrained to zero at baseline, site adjustment and a
participant random intercept. A parsimonious continuous-time correlation
structure accounts for within-participant dependence. Time-varying treatment
effects will be examined in sensitivity analysis.

Posterior marginal mean trajectories are integrated over Day 0--5 to obtain
the AUC ratio. An individual complete-case trapezoidal AUC is not the intended
primary analysis.

## 5. Missingness and intercurrent events

AUC does not itself solve missingness. The joint longitudinal likelihood uses
all observed samples and integrates over unobserved outcomes under missing at
random conditional on included predictors and observed outcomes. A separate
imputation of the primary outcome is not required.

Reasons for absent samples will distinguish logistical omission, blocked or
removed drains, discharge, rescue fibrinolysis, surgery, deterioration,
withdrawal and death. Delta-adjusted pattern-mixture analyses will assess
clinically plausible departures from missing at random.

Panel B examines loss of the Day-3 sample under MCAR, MAR after improvement and
MAR after deterioration, from 0% to 50%. Days 0, 1 and 5 are retained in this
focused comparison with Day-3 ANCOVA.

## 6. Main operating-characteristic scenarios

- planned recruitment fixed at 90 participants (45 per arm), allowing for five
  lost to follow-up per arm;
- all randomised participants with at least one suPAR measurement contribute to
  the longitudinal likelihood;
- total participants contributing suPAR data in Panel A: 60, 70, 80 and 90,
  with 80 used as the conservative design point;
- true Day 0--5 AUC reduction: 0%, 15%, 20%, 25% and 30%;
- immediate sustained treatment effect, never applied at baseline;
- primary decision rule `P(Delta < 0 | data) >= 0.95`; every simulated power
  and Type-I error estimate is the probability that this rule is met;
- 20,000 simulations per Panel A cell;
- 10,000 paired simulations per Panel B cell.

The historical-mean sensitivity uses 20,000 paired simulations per scenario,
post-baseline control-trajectory differences from -20% to +20%, and true
reductions of 0% and 20%. Baseline is matched in this focused conflict scenario.

## 7. Planned displays

The application figure contains four panels:

1. power versus total N contributing suPAR data for 15%, 20%, 25% and 30% AUC
   reductions, with recruitment fixed at 90 and 80 as the conservative design
   point;
2. model-based AUC versus Day-3 ANCOVA under MCAR and two MAR mechanisms, with
   Day-3 missingness up to 50%;
3. power under prespecified historical mean-borrowing sensitivities;
4. Type-I error under historical-current post-baseline trajectory differences.

## 8. Reproducibility

The design simulation fixes the covariance at the pilot/projected value and is
therefore an approximation rather than a simulation fit of the full final
hierarchical model. The final scripts use fixed seeds and save cell-level results with Monte Carlo
standard errors. Tests check positive definiteness, exact effect calibration,
no treatment effect at baseline, recovery under large samples and successful
analysis with missing Day-3 data.

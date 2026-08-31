# Simulation plan and assumptions register

## 1. Trial objective and estimand

The primary mechanistic question is whether dexamethasone reduces cumulative
pleural-fluid suPAR over time. The candidate primary estimand is the ratio of
model-estimated, time-standardised AUC between treatment and control over Days
0--3:

\[
R_{AUC} = \frac{AUC_{Tx}}{AUC_{Ctrl}}.
\]

Values below 1 favour dexamethasone. An AUC ratio of 0.80 represents a 20%
reduction in cumulative suPAR burden. A Days 0--5 estimand is a sensitivity
analysis until Day-5 availability and covariance assumptions are supported.

The treatment effect is

\[
\Delta = \log(R_{AUC}).
\]

The primary Bayesian success rule to calibrate is

\[
P(\Delta < 0 \mid data) \ge c,
\]

where `c = 0.95` for the primary rule. Thresholds of 0.90 and 0.975 are
sensitivity analyses. The
clinically relevant probability `P(Delta <= log(0.80) | data)` will be reported
separately and is not a competing primary rule.

## 2. Historical information

The historical cohort contains 29 complete Day 0, Day 1 and Day 3 trajectories.
On the natural-log scale:

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

These estimates are uncertain and are not treated as known truths. Variance and
correlation sensitivity scenarios are required.

## 3. Prior specification

The historical cohort contains no randomized dexamethasone comparison. It must
not inform the treatment-effect prior.

The initial treatment prior is:

\[
\Delta \sim N(0, 0.50^2),
\]

with prior-SD sensitivity values 0.25 and 0.50.

Historical borrowing is restricted to selected control-trajectory or nuisance
parameters. The preferred primary implementation is a robust mixture containing
25% historical and 75% weakly informative components:

- historical component capped at effective sample size (ESS) 10;
- weakly informative component carrying 75% initial prior weight;
- posterior mixture weight and posterior ESS reported;
- borrowing diminished when concurrent controls conflict with history.

Fixed power-prior weights `a0 = 0, 0.2, 0.4, 0.6` are sensitivity analyses.
For historical N=29 their nominal ESS values are 0, 5.8, 11.6 and 17.4.

If the historical participants are UK-derived, historical control means are
borrowed only into the UK control component. Kenyan control means are learned
from concurrent Kenyan participants. Variance/correlation borrowing across
countries is evaluated by sensitivity analysis rather than assumed.

## 4. Longitudinal model

The intended primary model is a Bayesian log-normal longitudinal mixed model
using actual elapsed sampling time. The initial mean model is piecewise linear
with knots at Days 1 and 3, treatment-by-time effects constrained to zero at
baseline, site/country intercept adjustment and a participant random intercept.

A random slope or continuous-time residual process will be added only if model
diagnostics and simulations show it is estimable with N=90.

Posterior marginal mean trajectories are integrated to obtain the AUC ratio.
An individual complete-case trapezoidal AUC is not the intended primary
analysis.

## 5. Missingness and intercurrent events

AUC does not itself solve missingness. The primary longitudinal model relies on
MAR conditional on observed outcomes and included covariates. Reasons for every
absent sample should be recorded, including:

- logistical omission;
- blocked drain;
- drain removal after improvement;
- discharge;
- rescue fibrinolysis or surgery;
- deterioration, withdrawal or death.

Drain removal, surgery and death may be intercurrent or structural events rather
than ordinary missing values. The clinical meaning of suPAR after these events
must be agreed before the final estimand is fixed.

Simulation scenarios:

- complete follow-up;
- MCAR 10%, 20%, 30%, 40%, 50%;
- MAR after improvement;
- MAR after deterioration/rescue;
- differential missingness by treatment;
- informative sampling intensity (odds ratio 1.5 or 2.0 per latent-SD);
- MNAR delta shifts of -0.50, -0.25, +0.25, +0.50 residual SD;
- structural observation termination at 10%, 20%, 30%.

The primary MNAR sensitivity analysis is a delta-adjusted pattern-mixture model.

## 6. Data-generating scenarios

### Sample size and treatment effect

- total N: 60, 70, 80, 90, 100, 120;
- true AUC reduction: 0%, 15%, 20%, 25%, 30%, 35%, 40%;
- immediate sustained, delayed and waning treatment shapes;
- treatment is never applied at baseline.

### Control distribution

- pilot covariance with variance multipliers 0.75, 1.00, 1.25 and 1.50;
- weak, moderate and pilot-like correlation structures;
- historical-to-current control drift from -30% to +30%;
- Day-5 mean: continued decline, plateau or rebound;
- Day-5 SD: 0.65, 0.75, 0.85.

### UK and Kenya

- country allocation 50:50 and 70:30;
- country control shifts 0, 0.3 and 0.6 residual SD;
- Kenya:UK residual-SD ratios 1.0, 1.25 and 1.5;
- homogeneous treatment ratios 0.80/0.80 and 0.70/0.70;
- modest HTE 0.70/0.85;
- strong HTE 0.70/1.00;
- qualitative reverse-direction stress test.

The overall effect is primary. Treatment-by-country HTE is exploratory. With
approximately 22--23 participants per treatment-country cell, the study cannot
reliably rule in or rule out clinically important HTE.

## 7. Operating characteristics

For every scenario report:

- Bayesian assurance: probability of satisfying the primary rule;
- false-positive probability under no treatment effect;
- probability of satisfying the 20% clinically relevant rule;
- bias and RMSE of the posterior effect;
- 95% credible-interval coverage and width;
- model convergence/fitting failure;
- prior/posterior mixture weight and borrowed ESS;
- probability that borrowing changes the go/no-go decision;
- Monte Carlo standard error.

## 8. Planned displays

The main four-panel figure contains:

1. power versus sample size, by true AUC reduction;
2. model-based AUC versus Day-3 ANCOVA under three missingness mechanisms;
3. power with no borrowing versus the selected robust mixture;
4. type-I error and power versus historical-current control conflict.

UK/Kenya HTE will be a supplementary display showing country-specific posterior
intervals and what N=90 can and cannot learn.

## 9. Decisions still requiring clinical agreement

1. Is Days 0--3 or Days 0--5 the clinically meaningful primary window?
2. What is the estimand after drain removal following improvement?
3. What is the estimand after rescue fibrinolysis, surgery or death?
4. Is the historical cohort wholly UK-derived and otherwise comparable to the
   proposed trial controls?
5. Which variables defining aetiology and care pathway will be measured
   consistently in the UK and Kenya?

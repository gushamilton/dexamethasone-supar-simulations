# Stage-2 statistical rationale and outline SAP

## Grant-ready text

The primary analysis will estimate the effect of dexamethasone on pleural-fluid
suPAR over time using a Bayesian longitudinal mixed-effects model. The model
will use each participant's available measurements at their actual sampling
times and will estimate the treatment contrast in model-based area under the
curve (AUC). The treatment effect will be reported as
`Delta = log(AUC_dexamethasone / AUC_control)`, where each AUC is the area under
the relevant group-specific geometric-mean suPAR trajectory. The primary
criterion for evidence that dexamethasone reduces suPAR AUC will be
`P(Delta < 0 | trial data) >= 0.95`; we will also report the posterior
probability of a clinically relevant reduction of at least 20%.

We considered using the previous suPAR cohort to inform the concurrent control
trajectory. In 20,000 simulations per scenario, conservative ESS-capped
borrowing produced only a modest power gain when the historical and current
controls were exchangeable. Greater apparent gains required stronger assumptions and were
accompanied by Type-I error inflation under prior--data conflict,
additional complexity, and uncertain exchangeability, particularly because the
new trial includes both UK and Kenyan participants. We will therefore not
borrow the historical control mean in the primary analysis. Concurrent
randomized controls will determine the control trajectory and treatment
contrast. The previous cohort will instead inform conservative priors for
nuisance parameters, including variability, within-participant correlation and
the association between baseline suPAR and subsequent AUC. These priors will be
sufficiently dispersed to allow the trial data to dominate when the new cohort
differs from the pilot cohort.

The operating-characteristic simulations are a design-stage Gaussian
approximation. They fix the covariance at the pilot value, with a prespecified
projection for Day 5; the final analysis will estimate covariance parameters
using the priors specified below.

Historical control-mean borrowing will be restricted to prespecified
sensitivity analyses. These will include a robust mixture or
switch-commensurate prior capped at an effective sample size of five, and fixed
power priors over a limited range of weights. Results will be presented with
and without borrowing. Mean borrowing will use observed pilot Days 0, 1 and 3
only; projected Day-5 information will not be borrowed. This makes the primary conclusion depend only on the
randomized comparison while still showing whether cautious borrowing would
materially alter the inference.

Potential UK--Kenya differences will be represented by prespecified fixed
country/site terms. Country-specific residual variability will be examined in
sensitivity analysis, and country-specific treatment contrasts will remain
exploratory rather than carrying separate confirmatory decision rules.

The longitudinal likelihood uses all observed suPAR measurements and naturally
accommodates unequal sampling schedules and incomplete outcome vectors under a
missing-at-random assumption conditional on observed treatment, country/site,
baseline suPAR and previous suPAR measurements. A separate imputation of the
primary outcome is therefore not required. Informative missingness will be
examined using pattern-mixture sensitivity analyses in which unobserved values
are shifted from their missing-at-random predictions under prespecified
clinically plausible scenarios. Sampling that stops because of drain removal,
discharge, surgery or death will be identified separately rather than treated
as ordinary intermittent missingness.

## Outline statistical analysis plan

### 1. Population and estimand

- Primary population: all randomized participants analysed according to their
  allocated treatment. Participants with any suPAR measurement remain in the
  longitudinal likelihood; a participant with baseline data only can inform
  nuisance parameters but provides no direct post-randomization treatment
  contrast. The number with no post-randomization sample will be reported by
  arm, with a prespecified sensitivity analysis for this group.
- Recruitment target: 90 participants (45 per arm), allowing for five lost to
  follow-up per arm. The design operating characteristics conservatively use
  40 contributing participants per arm, although every randomised participant
  with at least one suPAR measurement enters the longitudinal likelihood.
- Outcome window: Day 0 to Day 5, using actual elapsed time from randomization
  and sample collection.
- Summary measure: area under each group-specific model-based geometric-mean
  pleural-fluid suPAR trajectory over the outcome window.
- Treatment contrast: `Delta = log(AUC_dexamethasone / AUC_control)`, marginalized
  over the randomized trial's country/site distribution. Negative values favour
  dexamethasone.
- Primary decision: evidence of benefit if `P(Delta < 0 | data) >= 0.95`.
- Additional summaries: posterior median and 95% credible interval for Delta;
  AUC ratio and percentage reduction; `P(Delta <= log(0.80) | data)`. The
  supportive clinically relevant criterion is
  `P(Delta <= log(0.80) | data) >= 0.90`; it is not a second primary rule.

### 2. Primary longitudinal model

- Analyse log suPAR because the pilot distributions are right-skewed and the
  treatment effect is naturally interpreted as a ratio.
- Represent time flexibly using prespecified piecewise-linear terms or effects
  at the principal collection times, evaluated using actual elapsed sampling
  time.
- Include flexible time effects and one common post-baseline treatment effect,
  with no treatment difference at the pre-randomization baseline by
  construction. Time-varying treatment effects will be assessed in sensitivity
  analysis.
- Include baseline in the longitudinal response and constrain the randomized
  groups to have a common pre-randomization mean. Adjust for the randomization
  stratification variables; do not also enter the same baseline measurement as
  a separate covariate in the primary constrained longitudinal model.
- Include country or site as a fixed effect. With only two countries and four
  hospitals, a random country effect is not reliably estimable. If site is
  included as a fixed effect, a separate country main effect is redundant.
- Include a participant-specific random intercept to account for repeated
  measurements. Do not include a random treatment slope. A random time slope
  is also unnecessary unless supported by model diagnostics and prespecified
  as a sensitivity analysis.
- Use a parsimonious within-participant covariance model, such as continuous-time
  AR(1), with an alternative unstructured covariance sensitivity analysis if
  computationally stable.
- Derive posterior AUCs by integrating each posterior mean trajectory over the
  prespecified time window, then derive Delta from the treatment-specific AUCs.

### 3. Priors

- Treatment effect: the common post-baseline treatment coefficient has the
  fixed zero-centred prior `Normal(0, 1^2)` on the log scale.
- Control trajectory/mean: weakly informative priors in the primary analysis;
  no direct historical mean borrowing.
- Residual and random-effect standard deviations: pilot-centred log-Normal
  priors with log-scale SD 0.5.
- Correlation: pilot-centred Normal priors on the Fisher-z scale with SD 0.5.
- Baseline association: pilot-centred Normal prior on an unconstrained scale
  with SD 0.5.
- Prior predictive checks will be completed before unblinding and documented in
  the final SAP.

### 4. UK--Kenya heterogeneity

- Allow different mean suPAR levels between countries/sites through fixed
  intercepts in the primary model.
- Estimate one common randomized treatment effect for the primary estimand.
- Assess country-by-treatment interaction as exploratory because approximately
  45 participants per country cannot support a precise heterogeneous treatment
  effect estimate.
- Examine country-specific residual variance and country-by-time terms in
  sensitivity analyses if exploratory plots indicate substantial differences.
- Report country-specific posterior treatment contrasts descriptively, without
  separate confirmatory decision rules.

### 5. Missing outcomes and intercurrent events

- The primary Bayesian likelihood integrates over missing suPAR outcomes under
  MAR; no preliminary single or multiple imputation of suPAR is performed.
- Describe missingness by arm, time, country/site and reason.
- Include observed predictors of missingness in the outcome model where
  prespecified and measured before the missing value.
- Perform delta-adjusted pattern-mixture analyses for values missing after
  improvement and after deterioration. The direction and magnitude of each
  delta should be clinically justified and applied separately by reason for
  missingness.
- Prespecify how drain removal, surgery, discharge and death define the primary
  estimand. These events should not be silently treated as interchangeable MAR
  dropout. A hypothetical continuation analysis may be the primary mechanistic
  estimand, supported by observed-data and composite/worst-case sensitivities
  where clinically meaningful.
- If baseline covariates are missing, use a separate joint-model or multiple-
  imputation procedure; this is distinct from incomplete longitudinal outcomes.

### 6. Sensitivity analyses

1. Weak nuisance priors instead of pilot-centred nuisance priors.
2. Robust mixture or switch-commensurate historical control-mean prior, maximum
   effective sample size five.
3. Fixed power priors over the prespecified range.
4. Decision thresholds of 0.90 and 0.975.
5. Clinically relevant decision target:
   `P(AUC reduction >= 20% | data) >= 0.90`.
6. MNAR pattern-mixture scenarios.
7. Alternative covariance structure and country-specific residual variance.
8. Conventional Day-3 baseline-adjusted analysis for comparison with the
   model-based AUC analysis.
9. Time-varying treatment-by-time effects in place of the common post-baseline
   effect.

### 7. DMC and interim/go-no-go decisions

The independent DMC will review recruitment, protocol adherence, completeness
and timing of biomarker sampling, safety events and prespecified clinical
futility or feasibility criteria. The DMC will not select a prior after seeing
unblinded outcomes. Prior specifications, borrowing caps and Bayesian decision
thresholds will be fixed in the SAP before database lock. The primary decision
will be based on the no-historical-mean-borrowing analysis; borrowing analyses
will be supportive sensitivity analyses.

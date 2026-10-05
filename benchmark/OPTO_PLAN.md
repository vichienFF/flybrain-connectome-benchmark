# Pre-registration: incomplete optogenetic silencing in the silencing panel

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-10-06, before any simulation with partial silencing. Public timestamp = the GitHub commit adding this file.
Code: `benchmark/optosil.py`.

## Question
The fly silencing data in panel B (Shiu et al. 2024, GtACR1) come from optogenetic silencing, which is incomplete in real
neurons, whereas the model silences cells completely. Does incomplete silencing improve the model's silencing predictions?

## Rule (a mimic of incomplete silencing)
A silenced neuron may still spike, but each spike is transmitted with probability (1 - eff). eff = 1 reproduces the
model's complete silencing for all other neurons (verified with an unused seed, 999). A separate random stream is used, so
the sensory Poisson input is identical across conditions.

## Design
Model D0. The 17 panel-B tests (8 sugar, 9 water; labels, exclusions and the >= 20% rule as in benchmark v2) at the
input rates chosen for D0 in benchmark v4 (sugar 50 Hz, water 180 Hz), test seeds 70-79 (paired). eff = 0.9, 0.7, 0.5 are
simulated; eff = 1.0 and the baselines are the existing benchmark v4 D0 data (same seeds, rates and code path, identical
by construction).

## Criteria (fixed now)
- **Primary efficacy: eff = 0.7** (fixed in advance; the other levels are reported as a dose-response, not used for
  the verdict).
- **IMPROVED** if the number of passed tests at eff = 0.7 is >= 2 higher than at eff = 1.0 (11/17).
  **WORSE** if >= 2 lower. Otherwise **NO CLEAR CHANGE**.
- Pre-stated expectation: incomplete silencing weakens all knockout effects, so "not required" tests (e.g., Roundup under
  sugar, -89% at eff 1.0) may move toward passing, while "required" tests may drop below the 20% threshold; the net
  effect is uncertain.

## Limitations
The model's water baseline at these rates (21.4 Hz) was slightly outside the v4 matching range. Efficacy is uniform
across neurons; real efficacy varies with expression and light. 10 seeds. Execution: Kaggle (private, internet off).

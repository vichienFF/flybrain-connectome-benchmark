# Pre-registration: Benchmark v3 - silencing panel at a matched operating point

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-09-29, before any simulation of this design. Public timestamp = the GitHub commit that adds this file.
Code: `benchmark/bm3.py` (reuses neuron groups, labels and name lookup of `bm2.py`).
Preprint this follows: https://doi.org/10.5281/zenodo.23025560

## Question
In benchmark v2 (panel B), D0 and D0U (Usnea outputs excitatory) each passed 11/17 silencing tests, but at unmatched
operating points (water baseline MN9: D0 8.3 Hz, D0U 56.4 Hz). Does the Usnea re-signing still cost other tests when
both models are compared at the same baseline MN9 activity?

## Design
Models: D0, D0U (unchanged). Tests and labels: exactly the 17 panel-B tests of benchmark v2 (8 sugar, 9 water;
GtACR1 silencing results of Shiu et al. 2024 Suppl. Tables 2 and 5; same exclusions).
1. **Rate selection (seeds 40-44):** for each model and stimulus, simulate the stimulus alone at every rate of a fixed
   grid (sugar GRNs: 30, 35, 40, 45, 50, 55, 60, 70 Hz; water GRNs: 40, 60, 80, 100, 120, 160, 200, 250, 300 Hz) and choose
   the rate whose mean MN9 (average of left and right MN9, over seeds) is closest to **15 Hz** (ties: lower rate).
2. **Test (fresh seeds 50-59, 10 paired seeds):** at the chosen rate, baseline and each bilateral silencing.
   Call "required" if mean MN9 decreases by >= 20% (same rule as v2). Paired Wilcoxon p-values reported (descriptive).

## Criteria (fixed now)
- **Validity:** each model's test-seed baseline for each stimulus is within 15 Hz +/- 30% (10.5-19.5 Hz). If not, report
  that the operating points could not be matched with this grid, and treat the verdict as provisional.
- **Verdict on the Usnea hypothesis:**
  - SUPPORTED: D0U passes both Usnea tests AND loses none of the 15 non-Usnea tests that D0 passes.
  - WEAKENED: D0U loses >= 1 non-Usnea test that D0 passes (the re-signing has a cost even at matched activity).
  - NOT SUPPORTED: D0U loses no test but fails at least one Usnea test.
- Totals reported as pass counts out of 17 and out of 15 (excluding the two Usnea tests).

## Pre-stated expectation
From v2, the Rattle-under-sugar loss in D0U is expected to shrink at matched activity (operating-point effect seen in the
preprint), whereas the G2N-1-under-water result is uncertain. No expectation is stated for the final verdict.

## Limitations
Matching is on MN9 rate only; other parts of the circuit may still sit at different operating points. 10 seeds.
Execution: Kaggle (private, internet off), same code as in this commit.

## Disclosure (added before registration)
One code smoke test was run with an unused seed (999; water 100 Hz, Usnea silenced) only to check that the script runs;
its output was not used for any choice. No simulation with seeds 40-44 or 50-59 was run before registration.

Registered publicly: GitHub commit bca37c0, 2026-09-29 01:08:43 UTC (08:08 Bangkok); files verified identical.

## Result (scored 2026-09-29 09:10; Kaggle complete, 275/275 simulations per model; results/bm3_scores.json)
- Chosen rates: D0 sugar 55 Hz, water 160 Hz; D0U sugar 50 Hz, water 80 Hz.
- Test-seed baselines (target 15 Hz, valid range 10.5-19.5): D0 sugar 21.2, water 12.9; D0U sugar 22.8, water 8.8.
  **Validity criterion not met** (3 of 4 baselines outside the range; the grid was too coarse where MN9 rises steeply:
  D0 sugar 50 -> 55 Hz gave 9.0 -> 20.7 Hz; D0U water 80 -> 100 Hz gave 9.0 -> 30.5 Hz). Verdict therefore **provisional**.
- Pass counts: D0 11/17 (11/15 without Usnea tests); D0U 12/17 (10/15).
- **Verdict (registered rule): WEAKENED (provisional).** D0U passed both Usnea tests (sugar -28%, water -96%) and gained
  sugar:Clavicle, but lost two water tests that D0 passes: Bract (-10% vs -44%) and Rattle (-4% vs -45%); both are
  required in flies. The Rattle-under-sugar loss seen in v2 disappeared at closer operating points (as expected).
- Caveat: D0U's water baseline (8.8 Hz) was lower than D0's (12.9 Hz); the two lost tests are in the water condition.

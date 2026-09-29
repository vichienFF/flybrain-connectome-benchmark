# Pre-registration: Benchmark v4 - repeat of v3 with a finer rate grid (final model test of the Usnea hypothesis)

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-09-29, before any simulation of this design. Public timestamp = the GitHub commit that adds this file.
Code: `benchmark/bm4.py` (identical to `bm3.py` except the grid, the seeds and output names).

## Why
Benchmark v3 (commit bca37c0) gave the verdict WEAKENED, but its validity criterion failed: the rate grid was too coarse
where MN9 rises steeply, so baselines were not matched (D0 sugar 21.2, water 12.9; D0U sugar 22.8, water 8.8 Hz; target
15 Hz). Both tests lost by D0U (water: Bract, Rattle) were in the condition where D0U had the lower baseline.
The v3 grid results informed only the range of the new grid (seeds 40-44 are not reused). **This is the last model-based
test of the Usnea hypothesis in this project; the next step is the wet-lab test.**

## Design (same as v3 except)
- Grid: sugar GRNs 44, 46, 48, 50, 52, 54, 56, 58 Hz; water GRNs 80, 84, 88, 92, 96, 100, 120, 140, 160, 180, 200 Hz.
- Selection seeds 60-67 (8 seeds, fewer noise-driven choices); test seeds 70-79 (fresh, 10 paired seeds).
- Same 17 tests, labels, target (15 Hz), selection rule (closest; ties -> lower rate), call rule (>= 20% decrease).

## Criteria (unchanged from v3)
- Validity: each model's test-seed baseline for each stimulus within 10.5-19.5 Hz; otherwise the verdict is provisional.
- SUPPORTED: D0U passes both Usnea tests and loses none of the 15 non-Usnea tests that D0 passes.
- WEAKENED: D0U loses >= 1 non-Usnea test that D0 passes.
- NOT SUPPORTED: D0U loses no test but fails at least one Usnea test.
- The v4 verdict replaces the provisional v3 verdict if v4 is valid; if v4 is also not valid, both are reported and the
  model-based evidence is called inconclusive. No further grid refinement will be run.

## Disclosure
The v3 selection curves were seen before choosing this grid (they motivated the finer spacing). No simulation with seeds
60-67 or 70-79 was run before registration.

Registered publicly: GitHub commit 99a0ae4, 2026-09-29 02:24:49 UTC (09:24 Bangkok); files verified identical.

## Result (scored 2026-09-29 10:35; Kaggle complete, 342/342 simulations per model; results/bm4_scores.json)
- Chosen rates: D0 sugar 50, water 180 Hz; D0U sugar 46, water 92 Hz.
- Test-seed baselines (valid range 10.5-19.5 Hz): D0 sugar 11.5, water **21.4**; D0U sugar 14.2, water 17.1.
  **Validity criterion not met** (D0 water 21.4 Hz, just above the range). Per the registered rule, the model-based
  evidence on the Usnea hypothesis is called **INCONCLUSIVE** and no further refinement is run.
- Pass counts: D0 11/17 (11/15 without Usnea tests); D0U 11/17 (9/15). Rule outcome would be WEAKENED.
- Consistent across v2, v3 and v4: D0U passes both Usnea tests and gains sugar:Clavicle, but loses water:Bract
  (-7% vs -40%) and water:Rattle (-15% vs -40%), which are required in flies. sugar:Rattle, lost here at -19.4%
  (threshold 20%), was passed in v3; it is borderline.
- Conclusion for the project: model evidence neither confirms nor refutes an excitatory Usnea; a consistent cost in
  the water pathway (Bract, Rattle) must be reported alongside the Usnea gains. The decisive test is immunostaining.

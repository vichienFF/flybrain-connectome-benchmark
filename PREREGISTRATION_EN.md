# Pre-registration summary (English) — original time-stamped files are in Thai (`benchmark/BENCHMARK_PLAN.md`, `benchmark/D_PLAN.md`, `benchmark/NT_FLIP_PLAN.md`, `malecns/REPLICATE_PLAN.md`, `REVISION_PLAN.md`)
Times are local (UTC+7), self-recorded; amendments are appended with reasons. Not externally registered.

| Plan | Written before | Criteria (summary) | Outcome |
|---|---|---|---|
| Benchmark (16 tests) | any benchmark simulation (2026-09-25) | Table 1 of the preprint; F5 direction fixed from literature before viewing F5 | D0 11/16 |
| Model variants D1–D3 (round 1) | running any variant | better = more passes and no loss of D0 passes; synaptic scale calibrated on low-salt + taste-peg GRNs (cold dropped before any variant ran: runaway in D0) | D2 12/16; D1 5; D3 4 |
| D2 robustness + D4 (round 2) | running | D2 confirmed if new seeds ≥12 and ≥2/3 neighbours ≥12; D4 better if > D2 without losses | D2 confirmed; D4 worse |
| Near-threshold knockouts (round 3) | running (hypotheses from exploratory runs, overlapping seeds) | Rattle ≤−20% p<0.05 (D0); Usnea: D0U ≤−20% p<0.05 and D0 >−10% | both met |
| Water prediction (round 4) | running | D0U Usnea KO during water ≤−30% p<0.05; D0 not reduced | met |
| NT re-signing screen (44 types) | running | hit = score > D0 without losing guard tests; positive control Usnea ≥7/8 (threshold mis-specified; Usnea 6/8) | only Usnea improves; specificity only |
| Male CNS replication | running; amended after first run (runaway) → calibration on sugar response; extended 5→20 seeds after interim p=0.043 | DNge031 KO ≥+20% p<0.05; robustness: at every non-runaway, non-silent k | calibrated k: +59% (seeds 0–19), +55% fresh seeds 20–29; robustness across k failed |
| Revision analyses (after review) | running | DNge031 bilateral ≥+20% p<0.05 both stimulus sides; matched operating point on fresh seeds 20–39; male k-sensitivity and Roundup KO | DNge031 met (ipsilateral expectation not met); operating point: Rattle and Usnea criteria met; male robustness failed; male Roundup −94% |

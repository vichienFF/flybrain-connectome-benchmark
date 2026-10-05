# Pre-registration: a neuropeptide-like co-transmission rule for Usnea (D0P)

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-10-05, before any simulation with a non-zero peptide gain. Public timestamp = the GitHub commit adding this file.
Code: `benchmark/peptide.py` (rule), `benchmark/pep_pipeline.py` (calibration, scoring), `benchmark.py` (variant D0P), `bm4.py` (panel B).

## Motivation
Re-signing Usnea from GABAergic to excitatory (D0U) repaired the water test but consistently broke two water-pathway
silencing tests (Bract, Rattle) in benchmarks v2-v4. G. Sterne (personal communication, 2026-10-01) noted that knockdown
of Amon (a neuropeptide-processing enzyme) in Usnea reproduces the Usnea silencing phenotype, suggesting a neuropeptide
co-transmitter that the model lacks.

## Rule (a mimic, not a chemical simulation)
Usnea (CB0008, 2 neurons) keeps its GABA synapses unchanged. Each Usnea spike adds 1 to a shared peptide level c, which
decays with tau = 500 ms (fixed a priori as a typical slow-modulation time scale). Each direct postsynaptic partner j of
Usnea (336 neurons) receives a slow depolarising drive dv_j = gain * c * w_j * dt, with w_j = synapses from Usnea to j
divided by the maximum. Gain 0 reproduces D0 bit-for-bit (verified with an unused seed, 999).

## Steps and criteria
1. **Calibration (seeds 80-84).** Gains tested in increasing order: 0.0005, 0.001, 0.002, 0.004, 0.008, 0.016 mV/ms per
   unit. The smallest gain giving mean water->MN9 (F2 condition, 150 Hz) >= 20 Hz is chosen. If none does, the registered
   outcome is "the rule cannot rescue F2 within this range" and the study stops. F2 is therefore a calibration target and
   no longer an independent test.
2. **16-test benchmark (seeds 0-4).** Criterion "no cost": D0P loses none of the 11 tests passed by D0.
3. **Silencing panel B at a matched operating point** (identical to benchmark v4: same grid, selection seeds 60-67, test
   seeds 70-79, target MN9 15 Hz, valid range 10.5-19.5 Hz; D0 comparison = the v4 D0 data, same seeds).
   - **SUPPORTED:** D0P passes both Usnea tests and loses no non-Usnea test that D0 passes.
   - **WEAKENED:** D0P loses >= 1 non-Usnea test that D0 passes.
   - **NOT SUPPORTED:** no loss, but at least one Usnea test fails.
   - **Key prediction:** water:Bract and water:Rattle pass (they failed in D0U in every round).
   - If D0P baselines are outside the valid range, or D0's are (v4 D0 water baseline was 21.4 Hz), the verdict is
     provisional.

## Interpretation limits
Targets, time constant and the form of the drive are assumptions (receptor locations of any Usnea peptide are unknown).
A positive result means a slow excitatory co-transmitter is sufficient in principle, not that Usnea uses a specific
peptide. Execution: Kaggle (private, internet off), same code as this commit.

# Pre-registration: Benchmark v2 (expanded test set from published fly experiments)

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-09-28, before any simulation of this panel. Public timestamp = the GitHub commit that adds this file.
Code: `benchmark/bm2.py`. Neuron IDs: `literature/shiu2024_supp/sez_neurons_783.json` (see Mapping).

## Source of ground truth
Shiu et al. 2024, Nature 634:210 (CC BY 4.0), Supplementary Tables 2, 3, 5, 10 (file 41586_2024_7763_MOESM2_ESM.xlsx,
SHA1 4555ce42e3ba1bc3fceb24b58884c19cbeb3f612). Only the columns reporting REAL fly experiments are used as labels
(never the columns with Shiu's simulated predictions).

## Mapping of neuron IDs (done before this registration; no simulation involved)
Shiu's list (`sez_neurons.pickle`, 106 SEZ cell types, 372 neurons) uses FlyWire v630 IDs; 64 are absent in v783.
They were mapped by connectivity fingerprint (cosine over synapse counts with partners present in both versions):
mean cosine 0.996, min 0.884, no duplicates (`map_630_to_783.py`, report csv). One type is ambiguous (pSG1: best 0.884
vs runner-up 0.881 not of the same type) -> **pSG1 excluded**.

## Panel A - sufficiency screen (Suppl. Table 3; Shiu Fig. 2A-B)
Real label: fraction of flies showing MN9 activity/proboscis extension on optogenetic activation of the split-GAL4
line; label positive if > 0 (as Shiu). 106 types, 14 positive.
Exclusions (fixed now): pSG1 (mapping); **Usnea** (cross-lab conflict: positive in Shiu 2024, no PER in Jacobs et al.
2024 Cell Rep Fig. S3A). -> **104 types, 13 positive, 91 negative.**
Simulation: activate all neurons of the type at 50 Hz for 1 s (Poisson input as in the model), seeds 0-4.
Model call: positive if the seed-mean rate of either MN9 (left or right) > 0 Hz (Shiu's rule, primary);
secondary rule: >= 5 Hz.
Metric: **balanced accuracy** (mean of sensitivity and specificity) - raw accuracy is misleading with 91 negatives.

## Panel B - necessity (Suppl. Tables 2 and 5, GtACR1 silencing, real flies)
B-sugar (stimulus: sugar GRNs 50 Hz): required = Clavicle, FMIn, G2N-1, Rattle, Usnea; not required = Bract, Phantom,
Roundup. Excluded: Fdg (conflicting reports within the table), Zorro (no neuron IDs in Shiu's list). -> 8 tests.
B-water (stimulus: water GRNs 160 Hz): required = Bract, Clavicle, Rattle, Roundup, Usnea; not required = G2N-1,
Phantom, Tophat, Tulip. Excluded: Zorro (no IDs), MN9 (not tested). -> 9 tests.
Model call: "required" if silencing (both sides) lowers mean MN9 (both MN9 averaged, seeds 0-4 paired) by >= 20%
(Shiu's rule). A test passes if the model call equals the real label. If the baseline MN9 of a stimulus is < 5 Hz,
all tests of that stimulus are scored "not evaluable" (not pass).

## Panel C - antennal grooming (Suppl. Table 10 rows citing Hampel et al. 2015)
Activate all Johnston's organ neurons in the model groups (JO-CE + JO-F + JO-D_m) at 150 Hz, seeds 0-4.
Pass if seed-mean rate >= 5 Hz for each of: aBN1 (SAD093), aDN1 (DNg62), DN2 (DNge078). -> 3 tests.

## Models compared (unchanged, as in the preprint)
D0 (Shiu model), D0U (Usnea outputs excitatory), D2 (spike-frequency adaptation, k from calib_D2.json).

## Pre-stated expectations and criteria
1. **Replication:** D0 panel-A balanced accuracy >= 0.75 (Shiu's own counts imply about 0.85). If lower, report
   as a replication failure of the model/re-implementation on this panel.
2. **Usnea hypothesis cost check:** D0U panel-A balanced accuracy is not lower than D0 by more than 0.05.
3. **Usnea hypothesis gain:** D0U passes more panel-B tests than D0. Note: B-water Usnea and B-sugar Usnea are the
   same kind of evidence already used in the preprint -> report B with and without the two Usnea tests.
4. D2: no directional prediction (exploratory); report all numbers.
Totals are reported per panel (A as balanced accuracy; B and C as pass counts), never as one merged score.

## Caveat on circularity
Panel A labels come from a screen of an existing split-GAL4 collection (not chosen by the model) -> independent of
the model. But Shiu et al. used the same model (= D0) and already reported this comparison, so D0 on panel A is a
replication, not new evidence. The new information is the comparison between model variants.

## Disclosure (added before registration)
One smoke-test simulation was run to check the code after the criteria above were written: panel C, D0, seed 0 ->
aBN1 26, aDN1 5, DN2 5 spikes/s. Criteria were not changed. No other panel-A/B/C simulation was run before registration.

# Pre-registration: directed simplices in the Drosophila feeding circuit and their link to model activity

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-10-19, before any simplex count on real connectome data (only tool tests on random graphs and the size of
the circuit, 406 neurons / 6,712 edges, have been computed). Public timestamp = the GitHub commit adding this file.
Code: `topology/simplex.py` (pure Python; verified against brute force on 200 random graphs).
Idea source: Reimann et al. 2017, Front. Comput. Neurosci., doi:10.3389/fncom.2017.00048. Literature check:
`literature/TOPOLOGY_LIT_CHECK.md` (no prior simplex/Betti analysis of the adult FlyWire connectome found; closest:
Lin et al. 2024 Nature, 2-3 node motifs on FlyWire; Hess group, simplex counts in the larval Drosophila connectome).

## Definitions (fixed now)
- Graph: FlyWire v783 neurons of the feeding circuit used in our preprint (406 neurons, `benchmark/results/nt_audit.csv`),
  edge pre -> post if >= 5 synapses (`flywire/W_783_min5.npz`). Self-loops removed. Edge sign ignored (structure only).
- Directed k-simplex: k+1 neurons ordered so that every earlier neuron connects to every later one (Reimann 2017).
- Null model: 1,000 graphs from directed double-edge swaps (10 swaps per edge), which keep every neuron's in-degree and
  out-degree and the number of edges. Seeds 0-999.
- Model activity: the published LIF whole-brain model (D0, FlyWire v783), sugar input 50 Hz (rate used in benchmark v4),
  10 seeds (100-109), 10 s each, spikes of the 406 neurons recorded in 20 ms bins. Execution: Kaggle (private, internet off).

## Questions and criteria
**Q1 (structure). Is the feeding circuit enriched in high-dimensional directed simplices?**
For each dimension d >= 2, ratio R_d = real count / mean null count, and empirical p = fraction of nulls >= real.
- ENRICHED if R_d > 1 with p < 0.01 for every d from 3 to the highest dimension present in the real graph.
- NOT ENRICHED if p >= 0.01 at d = 3. Otherwise PARTIAL (report the dimensions).
Expectation (from other species): ENRICHED.

**Q3 (structure -> activity; primary new question). Do connected neurons that share higher-dimensional simplices fire
more synchronously in the model?**
For every edge u -> v with both neurons active (>= 1 Hz mean), dim(u,v) = highest dimension of a simplex containing
both; sync(u,v) = Pearson correlation of binned spike counts, averaged over seeds.
- Primary statistic: partial Spearman correlation between dim(u,v) and sync(u,v), controlling for synapse count of the
  edge and the two firing rates.
- SUPPORTED if the partial correlation > 0 and its 95% interval (bootstrap over the 10 seeds, 2,000 resamples) excludes 0.
  CONTRADICTED if the interval lies entirely below 0. Otherwise INCONCLUSIVE.
- Control: the same statistic with dim(u,v) taken from one null graph (seed 0) should give an interval including 0.

**Q2 (exploratory, no verdict).** Report, for the 17 silencing tests of benchmark v2 (panel B), each tested cell type's
highest simplex dimension and number of simplices containing it, next to the real-fly label and the model's result.
Used only to generate hypotheses (structure-only scores cannot distinguish sugar from water tests).

## Limitations
One circuit (selected from model activity in our earlier screen), one synapse threshold, one model, sugar input only.
Secondary analyses (excitatory-only edges; threshold >= 10 synapses; water input) will be labelled exploratory.
Betti numbers (cavities) are not computed in this round (no homology tool on this machine); planned for a later round.
Cross-species comparison is a later round (needs external data downloads).

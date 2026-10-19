# Literature check: directed simplices / Betti numbers on FlyWire (2026-10-16)
Idea (user, from Reimann et al. 2017, doi:10.3389/fncom.2017.00048): count directed cliques (simplices) and cavities (Betti
numbers) in FlyWire sub-circuits vs a degree-preserving null; uses: (1) prioritise circuits for virtual experiments,
(2) cross-species comparison, (3) test whether neurons in high-dimensional simplices respond more synchronously.

## Coverage (5-rule protocol)
- Europe PMC (incl. bioRxiv/medRxiv preprints): 3 queries (FlyWire/Drosophila/hemibrain/FAFB/male CNS x simplices/flag complex/Betti/clique).
- OpenAlex: all 404 works citing Reimann 2017, filtered by "Drosophila" (15), "FlyWire" (0), "connectome" (106 listed, titles screened).
- Web search: 5 queries (general, FlyWire-specific, directed q-analysis, Time-Optimal q-analysis, bioRxiv/arXiv/journal domains).
- Full texts read (keyword scan): Lin et al. 2024 Nature (FlyWire network statistics, PMC10402125); Unger, Krebs & Mueller
  arXiv 2205.06520; Windisch/Unger arXiv 2501.04596 and 2508.15583.
- Not done: arXiv API (rate-limited, 429); Google Scholar full search; authors not contacted.

## Findings
| Work | What it did | Overlap with our idea |
|---|---|---|
| Lin et al. 2024, Nature (FlyWire) | 2- and 3-node motifs, rich club, reciprocity, clustering; ER and other null models | **Partial**: only 2-3 node structures; no simplices of dimension >= 3, no Betti numbers |
| Unger et al. 2022/23 (arXiv 2205.06520) | simplex closing probabilities; data: C. elegans, BBP, mouse V1 | cites hemibrain only; no fly analysis |
| Windisch/Unger 2025 (arXiv 2501.04596) | fast directed q-analysis; benchmarks C. elegans, BBP | states max simplex size "25 (Drosophila)" in passing, source not specified (cites Phelps 2021 adult VNC); no fly analysis reported |
| Hess-Bellwald talk 2025 (NTU) | simplex overexpression grows with dimension: C. elegans, **Drosophila larva**, MICrONS, BBP | larva only (Winding 2023), not adult FlyWire |
| arXiv 2601.14482 (2026) | directed clique counts incl. Drosophila larva | larva only |
| Reimann 2017 + BBP follow-ups | simplices/cavities and correlated activity (structure-function) in cortical model | method source for use (3) |

## Verdict (INFERENCE, coverage above)
- No published study found that counts directed simplices / Betti numbers of the **adult FlyWire** connectome against a
  degree-preserving null, or links simplex membership to knockout effects or synchrony in a whole-brain model.
- Closest prior art: Lin 2024 (low-order motifs on FlyWire) and larva simplex counts (Hess group). Any claim of novelty must
  say "to our knowledge" and cite these. Re-check arXiv (API) and Google Scholar before pre-registration.
- Side finding (relevant to the feeding work): bioRxiv 10.64898/2025.12.14.694122 (Dec 2025) - FlyWire v783 spiking model of
  feeding suppression during escape. Must be read before any further feeding-circuit claims.

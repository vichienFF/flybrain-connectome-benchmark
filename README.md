# flybrain-connectome-benchmark

Code, pre-registrations and results for the preprint:
**"Benchmarking a whole-brain connectome model of *Drosophila* against experimental data: diagnosing knockout-prediction failures and a candidate excitatory role for the water-taste neuron Usnea"** — Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670). Preprint (Zenodo): https://doi.org/10.5281/zenodo.23025560.
## Contents
| Folder | What |
|---|---|
| `model/` | `flybrain.py` — NumPy/SciPy re-implementation of the Shiu et al. (2024) LIF whole-brain model (validated against the original Brian2 code) |
| `benchmark/` | Pre-registered benchmark (`BENCHMARK_PLAN.md`, `D_PLAN.md`, `NT_FLIP_PLAN.md` — criteria and time-stamped amendments), model variants (`flybrain_d.py`), calibration, scoring, NT audit, NT-flip screen; `results/` = raw simulation outputs |
| `knockout_screen/` | Whole-brain single-neuron knockout screen of the feeding circuit (round 1 all active neurons; round 2 10 seeds + FDR); `results/` = ranked tables |
| `malecns/` | Male CNS v1.0 model build, calibration and replication (`REPLICATE_PLAN.md`) |
| `make_figures.py`, `figures/` | Figures, generated only from result files |

## Data (not redistributed — download from the original sources)
- FlyWire v783 connectivity as used by Shiu et al. (2024): https://github.com/philshiu/Drosophila_brain_model (MIT) → place `Connectivity_783.parquet`, `Completeness_783.csv` in `shiu_model/`
- FlyWire neuron annotations v3.1.0: https://github.com/flyconnectome/flywire_annotations → `flywire/neuron_annotations.tsv`
- Male CNS v1.0 (Janelia, CC-BY): https://male-cns.janelia.org/download/ → `malecns/` (checksums in `malecns/SHA256SUMS.txt`)
Scripts currently use absolute paths (`D:\Research_FlyBrain\...`) in a few places; edit the path constants at the top of each script.

## Requirements
Python ≥ 3.10, numpy, scipy, pandas, pyarrow, matplotlib. One 1-s whole-brain simulation ≈ 50 s single-core, ~2.5 GB RAM.

## Credits and licenses
Code: MIT (this repository). Connectome data: FlyWire Consortium (CC-BY 4.0), Janelia FlyEM Male CNS (CC-BY). Model design: Shiu et al. (2024). Analyses and code were developed with the assistance of an AI system (Claude, Anthropic) under the author's direction.

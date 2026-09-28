# Pre-registration: NT-flip screen B (olfaction O1 and navigation N1)

Author: Vichien Fugsukjit (Independent Researcher, Bangkok; ORCID 0009-0005-2468-6670)
Written: 2026-09-28, before any simulation of this screen. Public timestamp = the GitHub commit that adds this file.
Code: `benchmark/nt_flip_b.py` (same commit). Model: D0 = Shiu et al. (2024) LIF model on FlyWire v783, unchanged.

## Question
In screen A (feeding circuit, 44 uncertain cell types), reversing the transmitter sign of one cell type (Usnea, CB0008) and no other repaired a failing benchmark test (water -> MN9). Screen B asks the same question for the two remaining D0 failures outside feeding:
- **O1 (olfaction):** DM1 ORN stimulation should drive DM1 PNs >= 20 Hz and >= 3x the mean of other active uniglomerular PNs. D0 fails (activity spreads to other PNs).
- **N1 (navigation):** after a 300 ms pulse to 6 neighbouring EPG neurons, their mean rate at 400-1000 ms should stay >= 5 Hz. D0 fails (0 Hz).

## Candidate cell types (fixed rule, applied once, before the screen)
1. Run D0 once (seed 0) with the O1 and N1 stimuli; collect every neuron with >= 1 spike.
2. Keep cell types with >= 1 active neuron that are flagged exactly as in screen A (`nt_audit.py`): FlyWire top_nt_conf < 0.7, OR FlyWire/Male CNS sign conflict.
3. Exclude the stimulated types (ORN_DM1, EPG) and cell types with no annotation.
4. If more than 80 types remain, keep the 80 with the most active neurons.
5. Controls: "none" (D0) and Usnea CB0008 as a negative control (feeding-circuit type; expected to change neither O1 nor N1).
The resulting list is saved to `results/flip_types_B.json` before the screen starts.

## Screen panel
Each flip (all output synapses of the type change sign) x 6 conditions x seeds 0-2, using the exact conditions and readouts of `benchmark.py`:
C1_none, F1_sugar, G1_jonCE, E1_lc4 (protection tests D0 passes), O1_ornDM1, N1_epg6_pulse (target tests). Scored with the `score.py` definitions (mean over 3 seeds).

## Criteria
- **Expected D0 on the panel:** C1, F1, G1, E1 pass; O1, N1 fail. If not, the screen is invalid and is reported as such.
- **Candidate:** a flip that makes O1 or N1 pass AND keeps C1, F1, G1, E1 passing.
- **Specificity:** if more than 10% of screened types rescue the same test, that rescue is called non-specific (the test is sensitive to global excitation/inhibition balance) and no single type is named for it.
- **Confirmation (required before any claim):** at most the 3 strongest specific candidates run the full 17-test benchmark on fresh seeds 20-24. A candidate is "supported" only if it scores >= 12/16 (D0 = 11) and loses no test that D0 passes. Then a literature check on that cell type is done before any public statement.
- **Prediction stated in advance:** N1 is probably not rescuable by one sign flip (a persistent bump needs slow recurrent dynamics that this model lacks). O1 may be rescued by flipping inhibitory local neurons, but such rescues may be non-specific.
- A null result (no specific candidate) will be reported as a result.

## Limitations
3 seeds per condition = screening, not inference. Dopamine/serotonin/octopamine are treated as excitatory in the model, so flipping them tests a model assumption, not biology. Execution: Kaggle (private, internet off) or local, same code.

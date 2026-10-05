"""
Neuropeptide-rule test (PEPTIDE_PLAN.md). Steps, run in order (Kaggle or local):
  python pep_pipeline.py calib   -> calib_D0P.json  (smallest gain in GRID with F2 water->MN9 >= 20 Hz, seeds 80-84)
  then: BM_VARIANT=D0P python benchmark.py 4 ; python score.py D0P          (16-test benchmark, seeds 0-4)
  then: BM_VARIANT=D0P python bm4.py run                                    (panel B, matched operating point, seeds 60-79)
  python pep_pipeline.py score DIR   -> results/pep_scores.json (D0 panel-B data = v4 file bm4_D0.jsonl in DIR)
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
GRID = [0.0005, 0.001, 0.002, 0.004, 0.008, 0.016]   # mV per ms per unit peptide level
SEEDS = range(80, 85)
TAU = 500.0


def calib():
    os.environ["BM_VARIANT"] = "D0"
    import benchmark as bm
    from peptide import BrainP
    bm.setup()
    W, ids = bm.B["b"].W, bm.B["b"].ids
    usnea = bm.B["g"]["usnea"]
    log = []
    for gain in GRID:
        bm.B["b"] = BrainP(W, ids, src=usnea, gain=gain, tau_ms=TAU)
        vals = [float(np.mean(bm.run_one(("F2_water", s))["mn9"])) for s in SEEDS]
        log.append(dict(gain=gain, f2_mean=round(float(np.mean(vals)), 2), per_seed=vals))
        print(log[-1], flush=True)
        if np.mean(vals) >= 20:
            break
    ok = log[-1]["f2_mean"] >= 20
    out = dict(gain=log[-1]["gain"] if ok else None, rescued=bool(ok), tau_ms=TAU, seeds=list(SEEDS), log=log)
    json.dump(out, open(os.path.join(os.environ.get("PEP_OUT", HERE), "calib_D0P.json"), "w"), indent=1)
    if not ok:
        raise SystemExit("no gain in grid rescued F2 -> stop (registered outcome)")


def score(d):
    import bm4
    R = {}
    for v in ("D0", "D0P"):
        X = [json.loads(l) for l in open(os.path.join(d, f"bm4_{v}.jsonl"), encoding="utf-8")]
        ch = bm4.select(X); r = dict(chosen_hz=ch, tests={})
        for s, labels in bm4.LABELS.items():
            get = lambda ko: np.array([np.mean(x["mn9"]) for sd in bm4.TEST_SEEDS for x in X
                                       if x["stim"] == s and x["hz"] == ch[s] and x["ko"] == ko and x["seed"] == sd])
            base = get("none"); r[f"baseline_{s}"] = round(float(base.mean()), 1)
            r[f"matched_{s}"] = bool(abs(base.mean() - bm4.TARGET) <= 0.3 * bm4.TARGET)
            for n, req in labels.items():
                k = get(n); chg = k.mean() / base.mean() - 1
                r["tests"][f"{s}:{n}"] = dict(real_required=req, change_pct=round(100 * chg, 1), pass_=bool((chg <= -0.20) == req))
        t = r["tests"]; r["pass_all"] = sum(x["pass_"] for x in t.values())
        r["pass_non_usnea"] = sum(x["pass_"] for k, x in t.items() if not k.endswith("usnea")); R[v] = r
    a, b = R["D0"]["tests"], R["D0P"]["tests"]
    lost = [k for k in a if not k.endswith("usnea") and a[k]["pass_"] and not b[k]["pass_"]]
    usnea_ok = all(b[k]["pass_"] for k in b if k.endswith("usnea"))
    R["verdict"] = dict(D0P_lost_non_usnea=lost, D0P_passes_both_usnea=usnea_ok,
                        key_water_bract_rattle_pass=[b["water:bract"]["pass_"], b["water:rattle"]["pass_"]],
                        D0P_valid=all(R["D0P"][f"matched_{s}"] for s in bm4.LABELS), D0_valid=all(R["D0"][f"matched_{s}"] for s in bm4.LABELS),
                        outcome="SUPPORTED" if (usnea_ok and not lost) else "WEAKENED" if lost else "NOT SUPPORTED")
    json.dump(R, open(os.path.join(HERE, "results", "pep_scores.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    for v in ("D0", "D0P"):
        x = R[v]; print(v, x["chosen_hz"], "baselines", x["baseline_sugar"], x["baseline_water"], "pass", x["pass_all"], "/17; non-Usnea", x["pass_non_usnea"], "/15")
    print("verdict:", R["verdict"])


if __name__ == "__main__":
    calib() if sys.argv[1] == "calib" else score(sys.argv[2])

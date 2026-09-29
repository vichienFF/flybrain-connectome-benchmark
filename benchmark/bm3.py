"""
Benchmark v3 (2026-09-29): silencing panel B at a matched operating point - criteria in BENCHMARK_V3_PLAN.md
Step 1 (selection, seeds 40-44): for each stimulus, pick the input rate from a fixed grid whose mean MN9 (both sides
averaged) is closest to 15 Hz, separately for each model. Step 2 (test, fresh seeds 50-59): 17 silencing tests.
  BM_VARIANT=D0|D0U BM3_WORKERS=4 python bm3.py run   -> bm3_<variant>.jsonl
  python bm3.py score                                  -> results/bm3_scores.json
"""
import os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bm2  # noqa: E402  (same neuron groups, labels and name lookup as benchmark v2)

TARGET = 15.0
GRID = {"sugar": [30, 35, 40, 45, 50, 55, 60, 70], "water": [40, 60, 80, 100, 120, 160, 200, 250, 300]}
SEL_SEEDS = range(40, 45)
TEST_SEEDS = range(50, 60)
LABELS = {"sugar": bm2.B_SUGAR, "water": bm2.B_WATER}


def sim(task):
    stim, hz, ko, seed = task
    bm2._setup(); S = bm2.S; b = S["b"]
    ix = S["grp"][stim]; r = np.full(len(ix), float(hz))
    sil = S["grp"][f"sez:{bm2.sez_key(ko)}"] if ko else None
    b.rng = np.random.default_rng(seed)
    c, _ = b.run(1000.0, rate_fn=lambda s: (ix, r), silence=sil)
    return dict(variant=bm2.VARIANT, stim=stim, hz=hz, ko=ko or "none", seed=seed, mn9=c[S["mn9"]].tolist())


def run():
    out = os.path.join(os.environ.get("BM3_OUT", HERE), f"bm3_{bm2.VARIANT}.jsonl")
    workers = int(os.environ.get("BM3_WORKERS", "4"))
    R = [json.loads(l) for l in open(out, encoding="utf-8")] if os.path.exists(out) else []
    done = {(x["stim"], x["hz"], x["ko"], x["seed"]) for x in R}
    f = open(out, "a", encoding="utf-8")

    def batch(tasks):
        todo = [t for t in tasks if (t[0], t[1], t[2] or "none", t[3]) not in done]
        with ProcessPoolExecutor(workers) as ex:
            for x in ex.map(sim, todo):
                f.write(json.dumps(x) + chr(10)); f.flush(); R.append(x)
        print(f"batch done: {len(todo)} sims", flush=True)

    batch([(s, hz, None, sd) for s in GRID for hz in GRID[s] for sd in SEL_SEEDS])
    chosen = select(R)
    print("chosen rates:", chosen, flush=True)
    batch([(s, chosen[s], ko, sd) for s in LABELS for ko in [None] + list(LABELS[s]) for sd in TEST_SEEDS])
    f.close()


def select(R):
    ch = {}
    for s in GRID:
        m = {hz: np.mean([np.mean(x["mn9"]) for x in R if x["stim"] == s and x["hz"] == hz and x["ko"] == "none" and x["seed"] in SEL_SEEDS])
             for hz in GRID[s]}
        ch[s] = min(GRID[s], key=lambda hz: (abs(m[hz] - TARGET), hz))
    return ch


def score():
    from scipy.stats import wilcoxon
    d = os.environ.get("BM3_OUT", HERE); res = {}
    for fn in sorted(os.listdir(d)):
        if not (fn.startswith("bm3_") and fn.endswith(".jsonl")):
            continue
        R = [json.loads(l) for l in open(os.path.join(d, fn), encoding="utf-8")]
        v = R[0]["variant"]; ch = select(R); r = dict(chosen_hz=ch, tests={})
        for s, labels in LABELS.items():
            get = lambda ko: np.array([np.mean(x["mn9"]) for sd in TEST_SEEDS for x in R
                                       if x["stim"] == s and x["hz"] == ch[s] and x["ko"] == ko and x["seed"] == sd])
            base = get("none"); r[f"baseline_{s}"] = round(float(base.mean()), 1)
            r[f"matched_{s}"] = bool(abs(base.mean() - TARGET) <= 0.3 * TARGET)
            for n, req in labels.items():
                ko = get(n); ch_pct = ko.mean() / base.mean() - 1 if base.mean() > 0 else float("nan")
                p = float(wilcoxon(ko, base).pvalue) if np.any(ko != base) else 1.0
                r["tests"][f"{s}:{n}"] = dict(real_required=req, change_pct=round(100 * ch_pct, 1), p=round(p, 4),
                                              pass_=bool((ch_pct <= -0.20) == req))
        t = r["tests"]
        r["pass_all"] = sum(x["pass_"] for x in t.values()); r["pass_non_usnea"] = sum(x["pass_"] for k, x in t.items() if not k.endswith("usnea"))
        res[v] = r
    if {"D0", "D0U"} <= set(res):
        a, b = res["D0"]["tests"], res["D0U"]["tests"]
        lost = [k for k in a if not k.endswith("usnea") and a[k]["pass_"] and not b[k]["pass_"]]
        gained = [k for k in a if not k.endswith("usnea") and not a[k]["pass_"] and b[k]["pass_"]]
        usnea_ok = all(b[k]["pass_"] for k in b if k.endswith("usnea"))
        res["verdict"] = dict(D0U_lost_non_usnea=lost, D0U_gained_non_usnea=gained, D0U_passes_both_usnea=usnea_ok,
                              valid=all(res[m][f"matched_{s}"] for m in ("D0", "D0U") for s in LABELS),
                              outcome="SUPPORTED" if (usnea_ok and not lost) else "WEAKENED" if lost else "NOT SUPPORTED")
    json.dump(res, open(os.path.join(HERE, "results", "bm3_scores.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    for v in ("D0", "D0U"):
        if v in res:
            x = res[v]; print(v, "rates", x["chosen_hz"], "baselines", x["baseline_sugar"], x["baseline_water"], "pass", x["pass_all"], "/17 non-Usnea", x["pass_non_usnea"], "/15")
    print("verdict:", res.get("verdict"))


if __name__ == "__main__":
    {"run": run, "score": score}[sys.argv[1]]()

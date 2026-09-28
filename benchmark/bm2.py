"""
Benchmark v2 (2026-09-28) - criteria in BENCHMARK_V2_PLAN.md (registered before running)
  BM_VARIANT=D0|D0U|D2  BM2_WORKERS=3 python bm2.py run     -> bm2_<variant>.jsonl (checkpoint per simulation)
  python bm2.py score                                        -> results/bm2_scores.json + table
Local runs pause during 07:05-07:45.
"""
import os, sys, json, time
from datetime import datetime
import numpy as np
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import benchmark as bm  # noqa: E402  (reads BM_VARIANT)

VARIANT = bm.VARIANT
SEEDS = range(5)
T_MS = 1000.0
SEZ = json.load(open(os.path.join(ROOT, "literature", "shiu2024_supp", "sez_neurons_783.json")))
LOW = {k.lower().replace("-", "_"): k for k in SEZ}
A_EXCLUDE = {"psg1", "usnea"}
A_LABELS_FN = os.path.join(ROOT, "literature", "shiu2024_supp", "csv", "Sup_Table_3_Predicted_MN9_vs._o.csv")
B_SUGAR = {"clavicle": True, "FMIn": True, "G2N_1": True, "rattle": True, "usnea": True,
           "bract": False, "phantom": False, "roundup": False}
B_WATER = {"bract": True, "clavicle": True, "rattle": True, "roundup": True, "usnea": True,
           "G2N_1": False, "phantom": False, "tophat": False, "tulip": False}
C_READ = {"aBN1": 720575940630907434, "aDN1": 720575940616185531, "DN2": 720575940629806974}


def a_labels():
    import csv
    R = list(csv.reader(open(A_LABELS_FN, encoding="utf-8")))[1:]
    return {r[0]: float(r[1]) > 0 for r in R if r[0].lower() not in A_EXCLUDE}


def sez_key(name):  # Table 3 names differ in case/hyphen from the ID list (bug fix 2026-09-28, see plan amendment)
    return LOW[name.lower().replace("-", "_")]


def conditions():
    C = {}
    for t in a_labels():
        C[f"A|{t}"] = ([(f"sez:{sez_key(t)}", 50.0)], [])
    C["B|sugar50|none"] = ([("sugar", 50.0)], [])
    for n in B_SUGAR:
        C[f"B|sugar50|{n}"] = ([("sugar", 50.0)], [f"sez:{sez_key(n)}"])
    C["B|water160|none"] = ([("water", 160.0)], [])
    for n in B_WATER:
        C[f"B|water160|{n}"] = ([("water", 160.0)], [f"sez:{sez_key(n)}"])
    C["C|jon150"] = ([("jon_all", 150.0)], [])
    return C


S = {}


def _setup():
    if S:
        return
    bm.setup(); b = bm.B["b"]; g = bm.B["g"]
    idx = {int(r): i for i, r in enumerate(b.ids.tolist())} if hasattr(b, "ids") else None
    if idx is None:
        from flybrain import load_shiu_783
        _, ids, _ = load_shiu_783(); idx = {int(r): i for i, r in enumerate(ids.tolist())}
    grp = dict(sugar=g["sugar"], water=g["water"])
    L = json.load(open(os.path.join(HERE, "ids_shiu_figures.json")))
    grp["jon_all"] = np.array([idx[x] for x in L["neu_JON_CE"] + L["neu_JON_F"] + L["neu_JON_D_m"] if x in idx])
    for k, v in SEZ.items():
        grp[f"sez:{k}"] = np.array([idx[int(x)] for x in v if int(x) in idx])
    S.update(b=b, grp=grp, mn9=bm.B["rd"]["mn9b"], cread={k: idx[v] for k, v in C_READ.items()}, C=conditions())


def run_one(task):
    cond, seed = task
    _setup(); b = S["b"]; stims, ko = S["C"][cond]
    ix = np.concatenate([S["grp"][k] for k, _ in stims]); r = np.concatenate([np.full(len(S["grp"][k]), hz) for k, hz in stims])
    sil = np.concatenate([S["grp"][k] for k in ko]) if ko else None
    b.rng = np.random.default_rng(seed); t = time.time()
    c, _ = b.run(T_MS, rate_fn=lambda s: (ix, r), silence=sil)
    return dict(variant=VARIANT, cond=cond, seed=seed, mn9=c[S["mn9"]].tolist(),
                grooming={k: int(c[i]) for k, i in S["cread"].items()}, n_active=int((c > 0).sum()), sec=round(time.time() - t, 1))


def _paused():
    if os.environ.get("KAGGLE_KERNEL_RUN_TYPE"):
        return False
    hm = datetime.now().hour * 60 + datetime.now().minute
    return 7 * 60 + 5 <= hm < 7 * 60 + 45


def run():
    out = os.path.join(os.environ.get("BM2_OUT", HERE), f"bm2_{VARIANT}.jsonl")
    workers = int(os.environ.get("BM2_WORKERS", "3")); deadline = float(os.environ.get("BM2_DEADLINE_S", str(10.3 * 3600)))
    done = set()
    if os.path.exists(out):
        for l in open(out, encoding="utf-8"):
            try:
                x = json.loads(l); done.add((x["cond"], x["seed"]))
            except Exception:
                pass
    todo = [(c, s) for c in conditions() for s in SEEDS if (c, s) not in done]
    print(f"{VARIANT}: todo {len(todo)} workers {workers}", flush=True)
    t0 = time.time(); n = 0
    with open(out, "a", encoding="utf-8") as f, ProcessPoolExecutor(workers) as ex:
        it, pend = iter(todo), set()
        while True:
            while len(pend) < workers * 2 and time.time() - t0 < deadline and not _paused():
                t = next(it, None)
                if t is None:
                    break
                pend.add(ex.submit(run_one, t))
            if not pend:
                if _paused():
                    time.sleep(60); continue
                break
            fin, pend = wait(pend, return_when=FIRST_COMPLETED)
            for fu in fin:
                f.write(json.dumps(fu.result()) + chr(10)); f.flush(); n += 1
                if n % 50 == 0:
                    print(f"{n}/{len(todo)} | {(time.time() - t0) / 60:.0f} min", flush=True)
    print(f"finished {n}/{len(todo)} in {(time.time() - t0) / 60:.1f} min", flush=True)


def score():
    d = os.environ.get("BM2_OUT", HERE); res = {}
    for fn in sorted(os.listdir(d)):
        if not (fn.startswith("bm2_") and fn.endswith(".jsonl")):
            continue
        R = [json.loads(l) for l in open(os.path.join(d, fn), encoding="utf-8")]
        v = R[0]["variant"]; by = {}
        for x in R:
            by.setdefault(x["cond"], {})[x["seed"]] = x
        mean_lr = lambda c: np.mean([by[c][s]["mn9"] for s in SEEDS], axis=0)  # [left, right] seed-mean
        lab = a_labels(); pred = {}; pred5 = {}
        for t in lab:
            m = mean_lr(f"A|{t}"); pred[t] = bool(m.max() > 0); pred5[t] = bool(m.max() >= 5)
        def ba(p):
            pos = [t for t in lab if lab[t]]; neg = [t for t in lab if not lab[t]]
            sens = np.mean([p[t] for t in pos]); spec = np.mean([not p[t] for t in neg])
            return dict(sensitivity=round(float(sens), 3), specificity=round(float(spec), 3), balanced_accuracy=round(float((sens + spec) / 2), 3),
                        false_neg=[t for t in pos if not p[t]], false_pos=[t for t in neg if p[t]])
        B = {}
        for stim, labels in (("sugar50", B_SUGAR), ("water160", B_WATER)):
            base = np.array([np.mean(by[f"B|{stim}|none"][s]["mn9"]) for s in SEEDS])
            for n, req in labels.items():
                ko = np.array([np.mean(by[f"B|{stim}|{n}"][s]["mn9"]) for s in SEEDS])
                if base.mean() < 5:
                    B[f"{stim}:{n}"] = dict(pass_=None, note="baseline < 5 Hz"); continue
                ch = ko.mean() / base.mean() - 1; call = ch <= -0.20
                B[f"{stim}:{n}"] = dict(real_required=req, change_pct=round(100 * ch, 1), pass_=bool(call == req))
        g = {k: float(np.mean([by["C|jon150"][s]["grooming"][k] for s in SEEDS])) for k in C_READ}
        Cp = {k: dict(hz=round(x, 1), pass_=x >= 5) for k, x in g.items()}
        bp = [x["pass_"] for x in B.values()]; bpn = [x["pass_"] for k, x in B.items() if not k.endswith("usnea")]
        res[v] = dict(A_primary=ba(pred), A_secondary_5Hz=ba(pred5), B=B, B_pass=f"{sum(p is True for p in bp)}/{len(bp)}",
                      B_pass_without_usnea=f"{sum(p is True for p in bpn)}/{len(bpn)}", C=Cp, C_pass=f"{sum(x['pass_'] for x in Cp.values())}/3")
    json.dump(res, open(os.path.join(HERE, "results", "bm2_scores.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    for v, r in res.items():
        print(f"== {v}: A BA {r['A_primary']['balanced_accuracy']} (sens {r['A_primary']['sensitivity']}, spec {r['A_primary']['specificity']}) "
              f"| A@5Hz BA {r['A_secondary_5Hz']['balanced_accuracy']} | B {r['B_pass']} (no Usnea {r['B_pass_without_usnea']}) | C {r['C_pass']}")


if __name__ == "__main__":
    {"run": run, "score": score}[sys.argv[1]]()

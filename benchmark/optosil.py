"""
optosil.py (2026-10-06) - incomplete optogenetic silencing (OPTO_PLAN.md).
Real GtACR1 silencing is incomplete; the model's default silencing holds the cell below threshold (complete).
Rule: a silenced neuron may still spike, but each of its spikes is transmitted with probability (1 - eff); eff = 1
reproduces complete silencing for all other neurons (verified). A separate RNG is used so the Poisson input stream is unchanged.
  BM_VARIANT=D0 OPTO_WORKERS=4 python optosil.py run   -> opto_D0.jsonl
  python optosil.py score DIR                           -> results/opto_scores.json (eff 1.0 = benchmark v4 D0 data, bm4_D0.jsonl in DIR)
"""
import os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bm2  # noqa: E402
from flybrain import Brain  # noqa: E402

EFFS = [0.9, 0.7, 0.5]          # 1.0 taken from benchmark v4 (identical by construction)
PRIMARY = 0.7
RATES = {"sugar": 50, "water": 180}   # D0 rates chosen in benchmark v4 (matched-operating-point design)
SEEDS = range(70, 80)
LABELS = {"sugar": bm2.B_SUGAR, "water": bm2.B_WATER}


class BrainO(Brain):
    def run(self, t_ms=1000.0, rate_fn=None, silence=None, eff=1.0, seed=0, **kw):
        p, n = self.p, self.W.shape[0]
        dt = p["dt"]; steps = int(t_ms / dt)
        dly = int(round(p["t_dly"] / dt)); ref = int(round(p["t_ref"] / dt))
        v = np.full(n, p["v0"]); g = np.zeros(n)
        ref_left = np.zeros(n, dtype=np.int32); counts = np.zeros(n, dtype=np.int32)
        buf = [np.empty(0, dtype=np.int64)] * dly
        a_s = np.exp(-dt / p["tau_syn"]); e_m = np.exp(-dt / p["tau_m"])
        k_vg = p["tau_syn"] / (p["tau_syn"] - p["tau_m"]) * (a_s - e_m)
        poi_w = p["w_syn"] * p["f_poi"]
        is_sil = np.zeros(n, bool)
        if silence is not None:
            is_sil[silence] = True
        drop_rng = np.random.default_rng(10_000_000 + seed)
        no_ref = np.zeros(n, bool)
        for s in range(steps):
            srcs = [rate_fn(s)] if rate_fn is not None else []
            active = ref_left == 0
            g_a = g[active]
            v[active] = p["v0"] + (v[active] - p["v0"]) * e_m + g_a * k_vg
            g[active] = g_a * a_s
            ref_left[~active] -= 1
            for pidx, r in srcs:
                no_ref[pidx] = True
                hit = self.rng.random(len(pidx)) < np.asarray(r) * dt / 1000.0
                if hit.any():
                    v[pidx[hit]] += poi_w
            spk = np.flatnonzero((v > p["v_th"]) & ((ref_left == 0) | no_ref))
            v[spk] = p["v_reset"]; g[spk] = 0.0
            ref_left[spk] = np.where(no_ref[spk], 0, ref)
            counts[spk] += 1
            tx = spk
            if is_sil[spk].any():
                sil = is_sil[spk]
                keep = ~sil | (drop_rng.random(len(spk)) >= eff)
                tx = spk[keep]
            arr = buf[s % dly]
            if arr.size:
                g += p["w_syn"] * np.asarray(self.WT[arr].sum(axis=0)).ravel()
            buf[s % dly] = tx
        return counts, None


S = {}


def _setup():
    if not S:
        bm2._setup(); b0 = bm2.S["b"]
        S["b"] = BrainO(b0.W, b0.ids)


def sim(task):
    stim, ko, eff, seed = task
    _setup(); b = S["b"]; G = bm2.S["grp"]
    ix = G[stim]; r = np.full(len(ix), float(RATES[stim]))
    b.rng = np.random.default_rng(seed)
    c, _ = b.run(1000.0, rate_fn=lambda s: (ix, r), silence=G[f"sez:{bm2.sez_key(ko)}"], eff=eff, seed=seed)
    return dict(stim=stim, ko=ko, eff=eff, seed=seed, mn9=c[bm2.S["mn9"]].tolist())


def run():
    out = os.path.join(os.environ.get("OPTO_OUT", HERE), "opto_D0.jsonl")
    done = set()
    if os.path.exists(out):
        done = {(x["stim"], x["ko"], x["eff"], x["seed"]) for x in map(json.loads, open(out, encoding="utf-8"))}
    todo = [(s, k, e, sd) for e in EFFS for s in LABELS for k in LABELS[s] for sd in SEEDS if (s, k, e, sd) not in done]
    print("todo", len(todo), flush=True)
    with open(out, "a", encoding="utf-8") as f, ProcessPoolExecutor(int(os.environ.get("OPTO_WORKERS", "4"))) as ex:
        for i, x in enumerate(ex.map(sim, todo), 1):
            f.write(json.dumps(x) + chr(10)); f.flush()
            if i % 50 == 0:
                print(i, "/", len(todo), flush=True)


def score(d):
    V4 = [json.loads(l) for l in open(os.path.join(d, "bm4_D0.jsonl"), encoding="utf-8")]
    O = [json.loads(l) for l in open(os.path.join(d, "opto_D0.jsonl"), encoding="utf-8")]
    mean = lambda X: float(np.mean([np.mean(x["mn9"]) for x in X]))
    res = {}
    for eff in [1.0] + EFFS:
        r = {}
        for s, labels in LABELS.items():
            base = mean([x for x in V4 if x["stim"] == s and x["hz"] == RATES[s] and x["ko"] == "none" and x["seed"] in SEEDS])
            for k, req in labels.items():
                if eff == 1.0:
                    X = [x for x in V4 if x["stim"] == s and x["hz"] == RATES[s] and x["ko"] == k and x["seed"] in SEEDS]
                else:
                    X = [x for x in O if x["stim"] == s and x["ko"] == k and x["eff"] == eff]
                assert len(X) == len(SEEDS), (eff, s, k, len(X))
                ch = mean(X) / base - 1
                r[f"{s}:{k}"] = dict(real_required=req, change_pct=round(100 * ch, 1), pass_=bool((ch <= -0.20) == req))
        res[str(eff)] = dict(tests=r, pass_all=sum(x["pass_"] for x in r.values()))
    p1, pp = res["1.0"]["pass_all"], res[str(PRIMARY)]["pass_all"]
    res["verdict"] = dict(complete=p1, primary_eff=PRIMARY, primary=pp,
                          outcome="IMPROVED" if pp >= p1 + 2 else "WORSE" if pp <= p1 - 2 else "NO CLEAR CHANGE")
    json.dump(res, open(os.path.join(HERE, "results", "opto_scores.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print({e: res[e]["pass_all"] for e in res if e != "verdict"}, res["verdict"])


if __name__ == "__main__":
    run() if sys.argv[1] == "run" else score(sys.argv[2])

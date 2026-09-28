"""
NT-flip screen B (2026-09-28) — criteria in NT_FLIP_B_PLAN.md (registered before running)
Reverse the transmitter sign of one uncertain cell type at a time; score the olfaction (O1) and navigation (N1)
tests plus 4 protection tests, using the exact conditions/readouts of benchmark.py (D0 model).
Steps:
  python nt_flip_b.py scope                     -> results/flip_types_B.json (fixed candidate rule, D0 seed 0)
  NB_SHARD=i/N NB_WORKERS=3 python nt_flip_b.py screen   -> nt_flip_b_s<i>of<N>.jsonl (checkpoint per simulation)
  python nt_flip_b.py score                     -> results/ntflip_b_scores.csv
Local runs pause during 07:05-07:45 (bot schedule).
"""
import os, sys, json, time
from datetime import datetime
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
os.environ.setdefault("BM_VARIANT", "D0")
import benchmark as bm  # noqa: E402
from flybrain import load_shiu_783, Brain  # noqa: E402

PANEL = ["C1_none", "F1_sugar", "G1_jonCE", "E1_lc4", "O1_ornDM1", "N1_epg6_pulse"]
SEEDS = range(3)
SIGN = {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "dopamine": 1, "serotonin": 1, "octopamine": 1}  # as nt_audit.py
TYPES_FN = os.path.join(HERE, "results", "flip_types_B.json")
S = {}


def _base():
    if "W" not in S:
        bm.setup()  # D0 brain + groups + readouts (benchmark.B)
        W, ids, ann = load_shiu_783()
        a = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids)
        S.update(W=W.tocsc(), ids=ids, ct=a.cell_type.values, d0=bm.B["b"], flip="none")
    return S


def _set_flip(flip):
    s = _base()
    if s["flip"] != flip:
        if flip == "none":
            bm.B["b"] = s["d0"]
        else:
            Wc = s["W"].copy()
            for j in np.flatnonzero(s["ct"] == flip):
                Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]] *= -1
            bm.B["b"] = Brain(Wc.tocsr(), s["ids"])
        s["flip"] = flip


def run_one(task):
    flip, cond, seed = task
    _set_flip(flip)
    r = bm.run_one((cond, seed)); r["flip"] = flip
    return r


# ---------- step 1: candidate types (rule fixed in the plan) ----------
def _active(cond, seed=0):
    b, g = bm.B["b"], bm.B["g"]; stims, _ = bm.CONDS[cond]; dt = b.p["dt"]
    def rate_fn(s):
        t = s * dt; idx, r = [], []
        for grp, hz, t0, t1 in stims:
            if t0 <= t < t1:
                idx.append(g[grp]); r.append(np.full(len(g[grp]), float(hz)))
        return (np.concatenate(idx), np.concatenate(r)) if idx else (np.empty(0, np.int64), np.empty(0))
    b.rng = np.random.default_rng(seed)
    c, _ = b.run(bm.T_MS, rate_fn=rate_fn)
    return np.flatnonzero(c)


def scope():
    s = _base()
    act = np.union1d(_active("O1_ornDM1"), _active("N1_epg6_pulse"))
    ids = s["ids"]
    df = pd.DataFrame({"root_id": ids[act].astype(str), "cell_type": s["ct"][act]}).dropna(subset=["cell_type"])
    df = df[~df.cell_type.isin(["ORN_DM1", "EPG"])]
    fw = pd.read_csv(os.path.join(ROOT, "flywire", "neuron_annotations.tsv"), sep="\t", low_memory=False,
                     dtype={"root_id": str})[["root_id", "top_nt", "top_nt_conf"]]
    df = df.merge(fw, on="root_id", how="left")
    mc = pd.read_feather(os.path.join(ROOT, "malecns", "body-annotations-male-cns-v1.0-minconf-0.5.feather"))[["bodyId", "flywireType"]]
    nt = pd.read_feather(os.path.join(ROOT, "malecns", "body-neurotransmitters-male-cns-v1.0.feather"))[["body", "consensus_nt"]]
    m = mc.merge(nt, left_on="bodyId", right_on="body").dropna(subset=["flywireType"])
    mt = m.groupby("flywireType").consensus_nt.agg(lambda x: x.mode().iat[0] if len(x.dropna()) else None).rename("male_nt")
    df = df.merge(mt, left_on="cell_type", right_index=True, how="left")
    s_fw = df.top_nt.map(SIGN); s_m = df.male_nt.map(SIGN)
    df["flag"] = (df.top_nt_conf < 0.7) | (s_fw.notna() & s_m.notna() & (s_fw != s_m))
    per = df.groupby("cell_type").agg(n_active=("root_id", "size"), flagged=("flag", "any"),
                                      top_nt=("top_nt", lambda x: x.mode().iat[0] if len(x.dropna()) else None),
                                      conf=("top_nt_conf", "mean"), male_nt=("male_nt", "first"))
    cand = per[per.flagged].sort_values("n_active", ascending=False)
    n_all = len(cand); cand = cand.head(80)
    types = cand.index.tolist()
    os.makedirs(os.path.dirname(TYPES_FN), exist_ok=True)
    json.dump(dict(types=types, negative_control="CB0008", n_flagged_before_cap=n_all, n_active_neurons=int(len(act)),
                   made=datetime.now().isoformat(timespec="minutes")), open(TYPES_FN, "w"), indent=1)
    cand.to_csv(os.path.join(HERE, "results", "flip_types_B_detail.csv"), encoding="utf-8-sig")
    print(f"active neurons {len(act)} | flagged types {n_all} | screened {len(types)}", flush=True)


# ---------- step 2: screen ----------
def _paused():
    if os.environ.get("KAGGLE_KERNEL_RUN_TYPE"):
        return False
    hm = datetime.now().hour * 60 + datetime.now().minute
    return 7 * 60 + 5 <= hm < 7 * 60 + 45


def screen():
    si, sn = map(int, os.environ.get("NB_SHARD", "0/1").split("/"))
    workers = int(os.environ.get("NB_WORKERS", "3"))
    deadline = float(os.environ.get("NB_DEADLINE_S", str(10.3 * 3600)))
    out = os.path.join(os.environ.get("NB_OUT", HERE), f"nt_flip_b_s{si}of{sn}.jsonl")
    T = json.load(open(TYPES_FN))
    flips = ["none", T["negative_control"]] + T["types"]
    mine = flips[si::sn]
    done = set()
    if os.path.exists(out):
        for l in open(out, encoding="utf-8"):
            try:
                r = json.loads(l); done.add((r["flip"], r["cond"], r["seed"]))
            except Exception:
                pass
    todo = [(f, c, s) for f in mine for c in PANEL for s in SEEDS if (f, c, s) not in done]
    print(f"shard {si}/{sn}: flips {len(mine)} todo {len(todo)} workers {workers}", flush=True)
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
                if n % 36 == 0:
                    print(f"{n}/{len(todo)} | {(time.time() - t0) / 60:.0f} min", flush=True)
    print(f"finished {n}/{len(todo)} in {(time.time() - t0) / 60:.1f} min", flush=True)


# ---------- step 3: score (definitions copied from score.py) ----------
def _panel_scores(R):
    by = lambda c: [r for r in R if r["cond"] == c]
    mean = lambda c, key: float(np.mean([np.mean(r[key]) for r in by(c)]))
    f1 = mean("F1_sugar", "mn9"); g1 = mean("G1_jonCE", "abn1")
    e1 = float(np.mean([max(r["gf"]) for r in by("E1_lc4")]))
    o = by("O1_ornDM1"); dm1 = float(np.mean([max(r["dm1pn"]) for r in o]))
    oth = float(np.mean([np.mean([x for x in r["upn"] if x > 0] or [0]) for r in o]))
    late = float(np.mean([np.sum(np.array(r["epg_late"])[np.array(r["epg6"])]) / 6 / 0.6 for r in by("N1_epg6_pulse")]))
    c1 = max(r["n_active"] for r in by("C1_none"))
    return dict(C1=c1 == 0, F1=f1 >= 20, G1=g1 >= 20, E1=e1 >= 20, O1=dm1 >= 20 and dm1 >= 3 * max(oth, 1e-9), N1=late >= 5,
                f1=round(f1, 1), g1=round(g1, 1), e1=round(e1, 1), dm1=round(dm1, 1), upn_other=round(oth, 1), n1_late=round(late, 2))


def score():
    R = []
    for fn in os.listdir(os.environ.get("NB_OUT", HERE)):
        if fn.startswith("nt_flip_b_s") and fn.endswith(".jsonl"):
            R += [json.loads(l) for l in open(os.path.join(os.environ.get("NB_OUT", HERE), fn), encoding="utf-8")]
    rows = []
    for flip in sorted({r["flip"] for r in R}):
        rr = [r for r in R if r["flip"] == flip]
        if len(rr) < len(PANEL) * len(SEEDS):
            continue
        rows.append(dict(flip=flip, **_panel_scores(rr)))
    df = pd.DataFrame(rows)
    df["protect_ok"] = df[["C1", "F1", "G1", "E1"]].all(axis=1)
    df["cand_O1"] = df.O1 & df.protect_ok; df["cand_N1"] = df.N1 & df.protect_ok
    df.to_csv(os.path.join(HERE, "results", "ntflip_b_scores.csv"), index=False, encoding="utf-8-sig")
    scr = df[~df.flip.isin(["none", "CB0008"])]
    print(df[df.flip.isin(["none", "CB0008"])].to_string(index=False))
    for t in ("O1", "N1"):
        k = int(scr[f"cand_{t}"].sum())
        print(f"{t}: rescued by {k}/{len(scr)} types ({100 * k / max(len(scr), 1):.0f}%) -> "
              f"{'NON-SPECIFIC (>10%)' if k > 0.1 * len(scr) else 'specific' if k else 'none'}:",
              scr.loc[scr[f'cand_{t}'], 'flip'].tolist())


if __name__ == "__main__":
    {"scope": scope, "screen": screen, "score": score}[sys.argv[1]]()

"""ทดสอบซ้ำ DNge031 + Usnea ในสมองตัวผู้ (ดู REPLICATE_PLAN.md)"""
import os, sys, json
import numpy as np, pandas as pd, scipy.sparse as sp
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from flybrain import Brain, P
K = float(os.environ["MALE_K"]) if os.environ.get("MALE_K") else json.load(open(os.path.join(HERE, "calib_male.json")))["k"] if os.environ.get("MALE_CALIB", "1") == "1" and os.path.exists(os.path.join(HERE, "calib_male.json")) else 1.0
PK = dict(P); PK["w_syn"] = P["w_syn"] * K; PK["f_poi"] = P["f_poi"] / K  # ปรับเทียบ (ดู REPLICATE_PLAN.md) — อินพุตรับรู้คงเดิม
S = {}

def setup():
    W = sp.load_npz(os.path.join(HERE, "model", "W_male.npz")).tocsr(); ids = np.load(os.path.join(HERE, "model", "ids.npy"))
    a = pd.read_parquet(os.path.join(HERE, "model", "ann.parquet"))
    full = pd.read_feather(os.path.join(HERE, "body-annotations-male-cns-v1.0-minconf-0.5.feather")).set_index("bodyId")
    side = full.rootSide.reindex(a.bodyId).fillna(full.somaSide.reindex(a.bodyId)).astype(str).values
    ft = a.flywireType.fillna("").astype(str).values; ty = a.type.fillna("").astype(str).values
    g = dict(sugar=np.flatnonzero((ft == "LB3") & (ty == "LB3c") & (side == "R")),
             water=np.flatnonzero((ft == "LB3") & (ty == "LB3a") & (side == "R")),
             dnge031=np.flatnonzero(ty == "DNge031"), usnea=np.flatnonzero(ft == "CB0008"), roundup=np.flatnonzero(ft == "CB0553"))
    if not len(g["sugar"]):
        g["sugar"] = np.flatnonzero((ft == "LB3") & (ty == "LB3c"))[::2]; g["water"] = np.flatnonzero((ft == "LB3") & (ty == "LB3a"))[::2]
    Wu = W.tocsc(copy=True)
    for j in g["usnea"]:
        Wu.data[Wu.indptr[j]:Wu.indptr[j + 1]] = np.abs(Wu.data[Wu.indptr[j]:Wu.indptr[j + 1]])
    S.update(W=W, Wu=Wu.tocsr(), ids=ids, g=g, mn9=np.flatnonzero(ty == "MN9"))

CONDS = {"M_sugar": ("orig", "sugar", []), "M_sugar_koDNge031": ("orig", "sugar", ["dnge031"]), "M_sugar_koRoundup": ("orig", "sugar", ["roundup"]),
         "M_water_orig": ("orig", "water", []), "M_water_orig_koUsnea": ("orig", "water", ["usnea"]),
         "M_water_uexc": ("uexc", "water", []), "M_water_uexc_koUsnea": ("uexc", "water", ["usnea"])}

def run_one(task):
    if not S: setup()
    cond, seed = task; model, stim, ko = CONDS[cond]
    b = S.get("b_" + model)
    if b is None:
        b = Brain(S["W"] if model == "orig" else S["Wu"], S["ids"], p=PK); S["b_" + model] = b
    st = S["g"][stim]; b.rng = np.random.default_rng(seed)
    c, _ = b.run(1000.0, rate_fn=lambda s: (st, 150.0), silence=np.concatenate([S["g"][k] for k in ko]) if ko else None)
    return dict(cond=cond, seed=seed, mn9=c[S["mn9"]].tolist(), n_active=int((c > 0).sum()), n_stim=int(len(st)))

if __name__ == "__main__":
    out = os.path.join(HERE, (f"replicate_k{K:.3f}.jsonl" if os.environ.get("MALE_K") else "replicate_calib.jsonl") if K != 1.0 else "replicate.jsonl")
    print("k =", K, flush=True)
    only = os.environ.get("REP_CONDS"); sa, sb = map(int, os.environ.get("REP_SEEDS", "0,5").split(","))
    tasks = [(c, s) for c in CONDS for s in range(sa, sb) if not only or c in only.split(",")]
    with open(out, "a" if os.environ.get("REP_SEEDS") else "w") as f, ProcessPoolExecutor(int(os.environ.get('REP_WORKERS', '2'))) as ex:
        for r in ex.map(run_one, tasks):
            f.write(json.dumps(r) + chr(10)); f.flush(); print(r["cond"], r["seed"], r["mn9"], r["n_active"], flush=True)

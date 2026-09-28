"""
Map Shiu 2024 SEZ neuron root IDs (FlyWire v630, sez_neurons.pickle) that are absent in v783, by connectivity fingerprint.
For each missing v630 neuron: vector of signed synapse counts to/from partners whose root ID exists in BOTH versions;
compare (cosine) with every v783 neuron's vector over the same partner set; take the best match.
Output: shiu2024_supp/sez_neurons_783.json (+ mapping report csv). Matches with cosine < 0.8 or a non-unique best are flagged.
"""
import os, sys, json, pickle
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
d = pickle.load(open(os.path.join(ROOT, "shiu_model", "sez_neurons.pickle"), "rb"))
c6 = pd.read_parquet(os.path.join(ROOT, "shiu_model", "2023_03_23_connectivity_630_final.parquet"))
c7 = pd.read_parquet(os.path.join(ROOT, "shiu_model", "Connectivity_783.parquet"))
print(c6.columns.tolist(), c7.columns.tolist())
pre, post, w = "Presynaptic_ID", "Postsynaptic_ID", "Connectivity"
ids6 = set(c6[pre]) | set(c6[post]); ids7 = set(c7[pre]) | set(c7[post])
shared = np.array(sorted(ids6 & ids7)); sidx = {r: i for i, r in enumerate(shared)}
print(f"v630 {len(ids6)} | v783 {len(ids7)} | shared IDs {len(shared)}")


def vectors(c, rows):
    """rows: list of root ids -> sparse matrix [len(rows), 2*len(shared)] (outputs then inputs, counts)"""
    ridx = {r: i for i, r in enumerate(rows)}
    o = c[c[pre].isin(ridx) & c[post].isin(sidx)]
    i_ = c[c[post].isin(ridx) & c[pre].isin(sidx)]
    r_ = np.r_[o[pre].map(ridx), i_[post].map(ridx)]
    c_ = np.r_[o[post].map(sidx), i_[pre].map(sidx) + len(shared)]
    v = np.r_[o[w].abs(), i_[w].abs()].astype(float)
    return sp.csr_matrix((v, (r_, c_)), shape=(len(rows), 2 * len(shared)))


missing = sorted({x for v in d.values() for x in v if x not in ids7})
print("missing v630 ids:", len(missing))
A = vectors(c6, missing)
cand = np.array(sorted(ids7))
B = vectors(c7, cand.tolist())
nA = np.sqrt(A.multiply(A).sum(1)).A1; nB = np.sqrt(B.multiply(B).sum(1)).A1
S = (A @ B.T).toarray() / np.outer(nA + 1e-12, nB + 1e-12)
best = S.argmax(1); top = S.max(1)
S2 = S.copy(); S2[np.arange(len(best)), best] = -1; second = S2.max(1)
m = {old: int(cand[b]) for old, b in zip(missing, best)}
rep = pd.DataFrame(dict(v630=missing, v783=[m[x] for x in missing], cos=top.round(3), second=second.round(3),
                        type=[next(k for k, v in d.items() if x in v) for x in missing]))
rep["already_in_type"] = [r.v783 in d[r.type] for r in rep.itertuples()]
rep["flag"] = (rep.cos < 0.8) | (rep.second > rep.cos - 0.05) | rep.already_in_type
new = {k: [x if x in ids7 else m[x] for x in v] for k, v in d.items()}
dup = sum(len(v) - len(set(v)) for v in new.values())
os.makedirs(os.path.join(HERE, "shiu2024_supp"), exist_ok=True)
json.dump({k: [str(x) for x in v] for k, v in new.items()}, open(os.path.join(HERE, "shiu2024_supp", "sez_neurons_783.json"), "w"), indent=1)
rep.to_csv(os.path.join(HERE, "shiu2024_supp", "map_630_to_783_report.csv"), index=False)
print(rep.cos.describe().round(3).to_string())
print(f"flagged {int(rep.flag.sum())}/{len(rep)} | duplicate ids after mapping {dup}")
print(rep[rep.flag].to_string(index=False))

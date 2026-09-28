"""
ขั้น D — ปรับเทียบตัวคูณความแรง synapse (k) ต่อตัวแปร ด้วยสิ่งเร้าที่ไม่อยู่ในชุดทดสอบ (ดู D_PLAN.md)
k คูณ w_syn และหาร f_poi (อินพุตรับรู้คงเดิม) · เป้า: เซลล์ที่ยิง (เฉลี่ยเรขาคณิต) = D0
ใช้: python calibrate_d.py D1|D2|D3  → calib_<variant>.json
"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from flybrain import load_shiu_783, P  # noqa: E402
from flybrain_d import BrainD  # noqa: E402

VARIANTS = {"D0": {}, "D1": dict(std=dict(U=0.5, tau_rec=300.0)), "D2": dict(sfa=dict(dth=2.0, tau_th=100.0)),
            "D3": dict(std=dict(U=0.5, tau_rec=300.0), sfa=dict(dth=2.0, tau_th=100.0)),
            # รอบยืนยัน (D_PLAN.md ส่วน "รอบ 2"): ค่าใกล้เคียงของ D2 + D4 = D2 + STD เฉพาะ ORN
            "D2a": dict(sfa=dict(dth=1.5, tau_th=100.0)), "D2b": dict(sfa=dict(dth=2.5, tau_th=100.0)),
            "D2c": dict(sfa=dict(dth=2.0, tau_th=150.0)),
            "D4": dict(std=dict(U=0.5, tau_rec=300.0, mask="ORN"), sfa=dict(dth=2.0, tau_th=100.0)),
            # รอบ 5 (D_PLAN.md): D2 + Usnea (CB0008) เป็นกระตุ้น
            "D2U": dict(sfa=dict(dth=2.0, tau_th=100.0))}
USNEA_FLIP = {"D2U"}
ANN = {}


def _ct(ids):
    if "ct" not in ANN:
        import pandas as pd
        a = pd.read_csv(os.path.join(ROOT, "flywire", "neuron_annotations.tsv"), sep="	", low_memory=False)
        ANN["ct"] = a.drop_duplicates("root_id").set_index("root_id").cell_type.reindex(ids)
    return ANN["ct"].fillna("").astype(str)


def make(W, ids, variant, k):
    p = dict(P); p["w_syn"] = P["w_syn"] * k; p["f_poi"] = P["f_poi"] / k
    if variant in USNEA_FLIP:  # กลับเครื่องหมาย output ของ Usnea เป็นบวก (เหมือน D0U ใน benchmark.py)
        Wc = W.tocsc(copy=True)
        for j in np.flatnonzero(_ct(ids).values == "CB0008"):
            Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]] = np.abs(Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]])
        W = Wc.tocsr()
    kw = {key: dict(val) for key, val in VARIANTS[variant].items()}
    if kw.get("std", {}).get("mask") == "ORN":  # เซลล์ดมกลิ่น (cell_type ขึ้นต้น ORN_) เท่านั้น
        kw["std"]["mask"] = _ct(ids).str.startswith("ORN_").to_numpy(bool)
    return BrainD(W, ids, p=p, **kw)


def calib_stims(ann):
    a = ann.drop_duplicates("root_id").set_index("root_id")
    return [np.flatnonzero(a.cell_sub_class.values == s) for s in ("low-salt", "taste peg")]  # cold ตัดออก (D0 ลุกลาม) — ดู D_PLAN.md


def score(W, ids, stims, variant, k):
    b = make(W, ids, variant, k); vals = []
    for st in stims:
        for seed in (0, 1):
            b.rng = np.random.default_rng(seed)
            c, _ = b.run(500.0, rate_fn=lambda s, st=st: (st, 150.0)); vals.append(max(1, int((c > 0).sum())))
    return float(np.exp(np.mean(np.log(vals)))), vals


def main():
    v = sys.argv[1]
    W, ids, ann = load_shiu_783()
    ann = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids).reset_index()
    stims = calib_stims(ann)
    tp = os.path.join(HERE, "calib_D0.json")
    if os.path.exists(tp):
        target = json.load(open(tp))["gm"]
    else:
        target, vals = score(W, ids, stims, "D0", 1.0); json.dump(dict(gm=target, vals=vals, k=1.0), open(tp, "w"))
    print(v, "target", round(target, 1), flush=True)
    lo, hi, log = 1.0, 8.0, []
    for it in range(7):
        k = float(np.sqrt(lo * hi)); gm, vals = score(W, ids, stims, v, k); log.append((k, gm, vals))
        print(v, it, f"k={k:.3f} gm={gm:.1f}", vals, flush=True)
        if gm < target: lo = k
        else: hi = k
    best = min(log, key=lambda x: abs(np.log(x[1] / target)))
    json.dump(dict(variant=v, k=best[0], gm=best[1], target=target, log=log), open(os.path.join(HERE, f"calib_{v}.json"), "w"))
    print(v, "BEST k", round(best[0], 3), "gm", round(best[1], 1), "target", round(target, 1), flush=True)


if __name__ == "__main__":
    main()

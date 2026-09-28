"""ปรับเทียบตัวคูณ w_syn (k) ของสมองตัวผู้: เซลล์ที่ยิงเมื่อกระตุ้นหวาน (LB3c ขวา 150 Hz, 500 ms, seed 0-1) ≈ 384 (ตัวเมีย)"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replicate as r
from flybrain import Brain, P
TARGET = 384
if __name__ == "__main__":
    r.setup(); st = r.S["g"]["sugar"]; lo, hi, log = 0.2, 1.0, []
    for it in range(7):
        k = float(np.sqrt(lo * hi)); p = dict(P); p["w_syn"] = P["w_syn"] * k; p["f_poi"] = P["f_poi"] / k
        b = Brain(r.S["W"], r.S["ids"], p=p); n = []
        for seed in (0, 1):
            b.rng = np.random.default_rng(seed); c, _ = b.run(500.0, rate_fn=lambda s: (st, 150.0)); n.append(int((c > 0).sum()))
        m = float(np.mean(n)); log.append((k, m, c[r.S["mn9"]].tolist())); print(it, f"k={k:.3f} active={m:.0f} mn9={c[r.S['mn9']].tolist()}", flush=True)
        if m > TARGET: hi = k
        else: lo = k
    best = min(log, key=lambda x: abs(np.log(x[1] / TARGET)))
    json.dump(dict(k=best[0], active=best[1], target=TARGET, log=log), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "calib_male.json"), "w"))
    print("BEST", best, flush=True)

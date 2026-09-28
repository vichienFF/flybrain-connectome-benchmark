"""
วิเคราะห์รอบ 2 (ยืนยันผู้ต้องสงสัย): KO vs baseline แบบจับคู่ seed ต่อ seed
- ผลต่อ MN9 (ซ้าย+ขวา) เป็น % ของ baseline seed เดียวกัน → ค่าเฉลี่ย, 95% CI, p (t-test จับคู่), q (Benjamini–Hochberg)
- p2_sig = q < 0.05 และ |ผลเฉลี่ย| ≥ 5% (เล็กกว่านี้ถือว่าไม่มีนัยทางชีววิทยา)
- เทียบ 2 สภาวะ: เซลล์ไหนสำคัญ "เฉพาะตอนมีขม" (Welch t-test ระหว่างสภาวะ, BH)
ใช้: python analyze_p2.py [โฟลเดอร์ผล=results/p2] → results/p2_rank_<cond>.csv, results/p2_compare.csv
"""
import os, sys, json, glob
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "p2")


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def per_seed(cond):
    rows = [json.loads(l) for f in glob.glob(os.path.join(SRC, f"ks_{cond}_p2_*.jsonl")) for l in open(f, encoding="utf-8")]
    if not rows:
        return None
    df = pd.DataFrame([dict(ko=r["ko"], seed=r["seed"], m=sum(r["mn9"])) for r in rows]).drop_duplicates(["ko", "seed"])
    base = df[df.ko < 0].set_index("seed").m
    ko = df[df.ko >= 0].copy()
    ko = ko[ko.seed.isin(base.index)]
    ko["d"] = 100 * (ko.m - ko.seed.map(base)) / ko.seed.map(base)
    return ko, base


def main():
    per = {}
    for cond in ("sugar", "sugar_bitter"):
        r = per_seed(cond)
        if r is None:
            print(cond, "no p2 results"); continue
        ko, base = r
        per[cond] = ko
        out = []
        for k, g in ko.groupby("ko"):
            d = g.d.values; n = len(d); m = d.mean(); se = d.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
            p = stats.ttest_1samp(d, 0).pvalue if n > 1 and d.std() > 0 else (1.0 if m == 0 else 0.0)
            t = stats.t.ppf(0.975, n - 1) if n > 1 else np.nan
            out.append(dict(ko=k, p2_n=n, p2_d_pct=round(m, 1), p2_ci_lo=round(m - t * se, 1), p2_ci_hi=round(m + t * se, 1), p2_p=p))
        o = pd.DataFrame(out)
        o["p2_q"] = bh(o.p2_p.fillna(1))
        o["p2_sig"] = (o.p2_q < 0.05) & (o.p2_d_pct.abs() >= 5)
        p1 = os.path.join(RES, f"p1_rank_{cond}.csv")
        if os.path.exists(p1):
            o = o.merge(pd.read_csv(p1, dtype={"root_id": str})[["ko", "root_id", "cell_type", "top_nt", "d_pct"]], on="ko", how="left")
        o = o.sort_values("p2_d_pct", key=lambda x: -x.abs())
        o.to_csv(os.path.join(RES, f"p2_rank_{cond}.csv"), index=False, encoding="utf-8-sig")
        print(f"\n[{cond}] baseline MN9 seeds {base.min()}–{base.max()} (n={len(base)}) · {len(o)} cells · confirmed {int(o.p2_sig.sum())}")
        cols = [c for c in ("root_id", "cell_type", "top_nt", "d_pct", "p2_d_pct", "p2_ci_lo", "p2_ci_hi", "p2_q", "p2_sig") if c in o]
        print(o[cols].head(30).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    if len(per) == 2:
        a, b = per["sugar"], per["sugar_bitter"]
        rows = []
        for k in sorted(set(a.ko) & set(b.ko)):
            da, db = a[a.ko == k].d.values, b[b.ko == k].d.values
            p = stats.ttest_ind(db, da, equal_var=False).pvalue if len(da) > 1 and len(db) > 1 else np.nan
            rows.append(dict(ko=k, sugar=round(da.mean(), 1), sugar_bitter=round(db.mean(), 1), diff=round(db.mean() - da.mean(), 1), p=p))
        c = pd.DataFrame(rows); c["q"] = bh(c.p.fillna(1)); c["bitter_specific"] = (c.q < 0.05) & (c["diff"].abs() >= 10)
        c = c.sort_values("diff", key=lambda x: -x.abs())
        c.to_csv(os.path.join(RES, "p2_compare.csv"), index=False, encoding="utf-8-sig")
        print(f"\n[เทียบสภาวะ] เซลล์ที่ผลต่างกันระหว่าง หวาน vs หวาน+ขม อย่างมีนัย: {int(c.bitter_specific.sum())}")
        print(c.head(12).to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()

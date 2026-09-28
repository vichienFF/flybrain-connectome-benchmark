"""
วิเคราะห์ผลรอบ 1 (knockout screen) → ตารางจัดอันดับ + รายชื่อผู้ต้องสงสัยสำหรับรอบ 2 + สร้าง kernel รอบ 2
ใช้: python analyze_p1.py [TOP=40]   (อ่านจาก results/p1/, เขียน results/p1_rank_<cond>.csv, results/cand_<cond>.json)
"""
import os, sys, json, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 40


def main():
    from flybrain import load_shiu_783
    _, ids, ann = load_shiu_783()
    ann = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids)
    cands = {}
    for cond in ("sugar", "sugar_bitter"):
        rows = [json.loads(l) for f in glob.glob(os.path.join(RES, "p1", f"ks_{cond}_p1_*.jsonl")) for l in open(f, encoding="utf-8")]
        if not rows:
            print(cond, "no results"); continue
        base = json.load(open(os.path.join(RES, "p1", f"ks_{cond}_baseline.json"), encoding="utf-8"))
        b0 = base[0]; bsum = sum(b0["mn9"])
        # ความแปรปรวนระหว่าง seed ของ baseline (ใช้เป็นเกณฑ์หยาบว่า "ต่างจริง" ก่อนรอบ 2)
        bs = [sum(r["mn9"]) for r in base]
        brate = dict(zip(b0["idx"], b0["cnt"]))
        # เซลล์ที่รับอินพุตตรง (ปิดแล้ว MN9 ลดแน่นอน) → ติดธงไว้ ไม่นับเป็นผู้ต้องสงสัยรอบ 2
        from knockout_screen import setup, B
        if not B:
            setup()
        stim = set(B["sugar"].tolist()) | (set(B["bitter"].tolist()) if cond == "sugar_bitter" else set())
        df = pd.DataFrame([dict(ko=r["ko"], mn9_R=r["mn9"][0], mn9_L=r["mn9"][1], n_active=r["n_active"]) for r in rows])
        df["d_mn9"] = df.mn9_R + df.mn9_L - bsum
        df["d_pct"] = (100 * df.d_mn9 / bsum).round(1)
        a = ann.iloc[df.ko.values]
        df["root_id"] = ids[df.ko.values].astype(str)
        for c in ("cell_type", "cell_class", "super_class", "side", "top_nt"):
            df[c] = a[c].values
        df["base_spikes"] = df.ko.map(brate).fillna(0).astype(int)
        # ชื่อในงานวิจัยเดิม (Shiu et al. 2022 eLife / 2024 Nature: shiu_model/sez_neurons.pickle) — ใช้แยก "รู้แล้ว" vs "อาจใหม่"
        import pickle
        sez = pickle.load(open(os.path.join(ROOT, "shiu_model", "sez_neurons.pickle"), "rb"))
        df["lit_name"] = df.root_id.map({str(r): k for k, v in sez.items() for r in v})
        extra = {"CB0051": "sternum", "CB0248": "billiards/specter", "DNge080": "rounddown"}  # จาก Male CNS synonyms (Shiu 2022)
        df["lit_name"] = df.lit_name.fillna(df.cell_type.map(extra))
        df["is_input"] = df.ko.isin(stim)
        df["is_target"] = df.ko.isin(B["mn9"].tolist())  # MN9 เอง (ตัวที่วัด) — ปิดแล้วลดแน่นอน ไม่ใช่การค้นพบ
        df = df.reindex(df.d_mn9.abs().sort_values(ascending=False).index)
        df.to_csv(os.path.join(RES, f"p1_rank_{cond}.csv"), index=False, encoding="utf-8-sig")
        print(f"\n[{cond}] baseline MN9 R+L seed0={bsum} (3 seeds: {bs}) · KO {len(df)} cells")
        print(df.head(25)[["root_id", "cell_type", "super_class", "top_nt", "is_input", "is_target", "base_spikes", "d_mn9", "d_pct"]].to_string(index=False))
        cands[cond] = df[~df.is_input & ~df.is_target].head(TOP).ko.tolist()
    # รอบ 2: ผู้ต้องสงสัยรวมของทั้ง 2 สภาวะ รันในทั้งสองสภาวะ (เทียบได้ว่าเซลล์ไหนสำคัญเฉพาะตอนมีขม)
    union = sorted(set().union(*cands.values()))
    json.dump([int(x) for x in union], open(os.path.join(RES, "cand_p2.json"), "w"))
    print("pass-2 candidates (union):", len(union))
    print("\nwrote ranks + candidates to", RES)


if __name__ == "__main__":
    main()

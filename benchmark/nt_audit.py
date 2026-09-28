"""
ข้อ 2: ตรวจ "สารสื่อประสาทที่อาจทายผิด" ในวงจรการกิน (ต่อยอดจากกรณี Usnea)
ขอบเขต: ทุกเซลล์ที่ยิงใน baseline ของ knockout screen (หวาน / หวาน+ขม) · FlyWire v783 vs Male CNS v1.0 (จับคู่ด้วย flywireType)
ธง: (ก) FlyWire top_nt_conf < 0.7  (ข) เครื่องหมาย (กระตุ้น/ยับยั้ง) ไม่ตรงกันระหว่าง 2 ชุดข้อมูล
จัดอันดับตามความสำคัญ: |ผล knockout| สูงสุดจากรอบ 1 (หวาน หรือ หวาน+ขม)
→ results/nt_audit.csv
"""
import os, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SIGN = {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "dopamine": 1, "serotonin": 1, "octopamine": 1}  # ตามสมมติฐานโมเดล


def main():
    ks = os.path.join(ROOT, "knockout_screen", "results")
    r = pd.concat([pd.read_csv(os.path.join(ks, f"p1_rank_{c}.csv"), dtype={"root_id": str}).assign(cond=c)
                   for c in ("sugar", "sugar_bitter")])
    r = r[~r.is_input & ~r.is_target]
    eff = r.groupby("root_id").agg(cell_type=("cell_type", "first"), side=("side", "first"),
                                    max_abs_ko=("d_pct", lambda x: x.abs().max()), lit_name=("lit_name", "first")).reset_index()
    fw = pd.read_csv(os.path.join(ROOT, "flywire", "neuron_annotations.tsv"), sep="\t", low_memory=False,
                     dtype={"root_id": str})[["root_id", "top_nt", "top_nt_conf", "known_nt"]]
    eff = eff.merge(fw, on="root_id", how="left")
    mc = pd.read_feather(os.path.join(ROOT, "malecns", "body-annotations-male-cns-v1.0-minconf-0.5.feather"))[["bodyId", "flywireType"]]
    nt = pd.read_feather(os.path.join(ROOT, "malecns", "body-neurotransmitters-male-cns-v1.0.feather"))[["body", "consensus_nt", "celltype_predicted_nt_confidence"]]
    m = mc.merge(nt, left_on="bodyId", right_on="body").dropna(subset=["flywireType"])
    mt = m.groupby("flywireType").agg(male_nt=("consensus_nt", lambda x: x.mode().iat[0] if len(x.dropna()) else None),
                                      male_conf=("celltype_predicted_nt_confidence", "mean")).reset_index()
    eff = eff.merge(mt, left_on="cell_type", right_on="flywireType", how="left").drop(columns="flywireType")
    s_fw = eff.top_nt.map(SIGN); s_m = eff.male_nt.map(SIGN)
    eff["low_conf"] = eff.top_nt_conf < 0.7
    eff["sign_conflict"] = s_fw.notna() & s_m.notna() & (s_fw != s_m)
    eff["flag"] = eff.low_conf | eff.sign_conflict
    eff = eff.sort_values(["flag", "max_abs_ko"], ascending=[False, False])
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    eff.to_csv(os.path.join(HERE, "results", "nt_audit.csv"), index=False, encoding="utf-8-sig")
    print(f"เซลล์ในวงจรการกินที่ตรวจ: {len(eff)} · NT มั่นใจต่ำ (<0.7): {int(eff.low_conf.sum())} · เครื่องหมายขัดกัน FlyWire↔Male: {int(eff.sign_conflict.sum())}")
    cols = ["cell_type", "side", "lit_name", "top_nt", "top_nt_conf", "male_nt", "male_conf", "sign_conflict", "max_abs_ko"]
    print(eff[eff.flag].head(25)[cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()

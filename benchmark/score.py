"""
ให้คะแนนขั้น A ตามเกณฑ์ใน BENCHMARK_PLAN.md (ตั้งไว้ก่อนรัน) → scores.json + ตาราง
F5 / G2: ทิศทางที่คาดกำหนดใน EXPECT_DIR (เติมจากผลค้นวรรณกรรมก่อนดูผลจำลอง)
"""
import os, json
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
EXPECT_DIR = {"F5": "lower", "G2": "lower"}  # F5: Shiu 2024 Fig3d,e (แมลงจริง: Ir94e กด PER) — กำหนด 2026-09-25 ~16:05 ก่อนดูผลจำลอง


def main(variant="D0"):
    fn = "results_raw.jsonl" if variant == "D0" else f"results_raw_{variant}.jsonl"
    R = [json.loads(l) for l in open(os.path.join(HERE, fn), encoding="utf-8")]
    df = pd.DataFrame(R)
    hz = lambda cond, key, i=0: np.array([r[key][i] for r in R if r["cond"] == cond], float)  # 1 s → Hz
    m = lambda cond, key="mn9": hz(cond, key).mean()
    out = []
    def add(code, system, value, crit, ok, note=""):
        out.append(dict(test=code, system=system, value=value, criterion=crit, passed=bool(ok) if ok is not None else None, note=note))
    add("C1", "ควบคุม", int(df[df.cond == "C1_none"].n_active.max()), "0 เซลล์ยิง", df[df.cond == "C1_none"].n_active.max() == 0)
    f1 = m("F1_sugar"); add("F1", "การกิน", round(f1, 1), "MN9 ≥ 20 Hz", f1 >= 20)
    f2 = m("F2_water"); add("F2", "การกิน", round(f2, 1), "MN9 ≥ 20 Hz", f2 >= 20)
    f3 = m("F3_bitter"); add("F3", "การกิน", round(f3, 1), "MN9 < 5 Hz", f3 < 5)
    f4 = m("F4_sugar_bitter"); add("F4", "การกิน", f"{f4:.1f} ({100 * (f4 / f1 - 1):+.0f}%)", "ลด ≥ 30% จาก F1", f4 <= 0.7 * f1)
    f5 = m("F5_sugar_ir94e"); d5 = 100 * (f5 / f1 - 1)
    ok5 = None if EXPECT_DIR["F5"] is None else (d5 < 0 if EXPECT_DIR["F5"] == "lower" else d5 > 0)
    add("F5", "การกิน", f"{f5:.1f} ({d5:+.0f}%)", f"ทิศทาง: {EXPECT_DIR['F5'] or 'รอวรรณกรรม'}", ok5)
    xs = [25, 50, 100, 150, 200]; ys = [m(f"F6_sugar_{x}") if x != 150 else f1 for x in xs]
    rho = stats.spearmanr(xs, ys).statistic
    add("F6", "การกิน", f"ρ={rho:.2f} · {[round(y) for y in ys]}", "Spearman ρ ≥ 0.9", rho >= 0.9)
    for code, c in (("F7", "F7_sugar_koRattle"), ("F8", "F8_sugar_koUsnea")):
        v = m(c); add(code, "การกิน", f"{v:.1f} ({100 * (v / f1 - 1):+.0f}%)", "ลด ≥ 20% จาก F1", v <= 0.8 * f1)
    f9 = m("F9_fdg"); add("F9", "การกิน", round(f9, 1), "MN9 ≥ 20 Hz", f9 >= 20)
    if any(r["cond"] == "F10_sugar_koRoundup" for r in R):
        v = m("F10_sugar_koRoundup")
        add("F10*", "การกิน (post-hoc)", f"{v:.1f} ({100 * (v / f1 - 1):+.0f}%)", "ลด < 20% (แมลงจริงไม่ลด)", v > 0.8 * f1)
    g1 = m("G1_jonCE", "abn1"); g2 = m("G2_jonF", "abn1")
    add("G1", "ทำความสะอาด", round(g1, 1), "aBN1 ≥ 20 Hz", g1 >= 20)
    add("G2", "ทำความสะอาด", round(g2, 1), "aBN1 ต่ำกว่า G1", g2 < g1 if EXPECT_DIR["G2"] == "lower" else None)
    for code, c in (("E1", "E1_lc4"), ("E2", "E2_lplc2")):
        v = np.mean([max(r["gf"]) for r in R if r["cond"] == c]); add(code, "หนีภัย", round(v, 1), "Giant Fiber ≥ 20 Hz", v >= 20)
    o = [r for r in R if r["cond"] == "O1_ornDM1"]
    dm1 = np.mean([max(r["dm1pn"]) for r in o])
    others = np.mean([np.mean([x for x in r["upn"] if x > 0] or [0]) for r in o])
    add("O1", "กลิ่น", f"DM1 PN {dm1:.0f} Hz vs uPN อื่น {others:.1f} Hz", "≥ 20 Hz และ ≥ 3×", dm1 >= 20 and dm1 >= 3 * max(others, 1e-9))
    n = [r for r in R if r["cond"] == "N1_epg6_pulse"]
    late = np.mean([np.sum(np.array(r["epg_late"])[np.array(r["epg6"])]) / 6 / 0.6 for r in n])  # Hz เฉลี่ยของ 6 ตัว ช่วง 400–1000 ms
    add("N1", "นำทาง", f"{late:.1f} Hz หลังหยุดกระตุ้น", "bump คงอยู่ ≥ 5 Hz (คาดว่าตก)", late >= 5)
    s = pd.DataFrame(out)
    json.dump(out, open(os.path.join(HERE, "scores.json" if variant == "D0" else f"scores_{variant}.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(s.to_string(index=False))
    sc = s.dropna(subset=["passed"]).groupby("system").passed.agg(["sum", "count"])
    print("\nคะแนนตามระบบ:\n", sc.to_string())


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "D0")

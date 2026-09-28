"""สำรวจ (post-hoc, ไม่ใช่ส่วนของเกณฑ์): ผลปิด Rattle/Usnea/Roundup ที่ความหวานต่ำ (ใกล้เกณฑ์ PER) vs สูง — โมเดลเดิม D0"""
import os, sys, json
os.environ["BM_VARIANT"] = "D0"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import benchmark as bm
for hz in (50, 75):
    bm.CONDS[f"X_sugar{hz}"] = ([("sugar", hz, 0, bm.T_MS)], [])
    for ko in ("rattle", "usnea", "roundup"):
        bm.CONDS[f"X_sugar{hz}_ko{ko}"] = ([("sugar", hz, 0, bm.T_MS)], [ko])
bm.CONDS = {k: v for k, v in bm.CONDS.items() if k.startswith("X_")}
bm.OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_lowrate.jsonl")
if __name__ == "__main__":  # รันในโปรเซสเดียว (worker แยกจะไม่เห็นเงื่อนไข X_ ที่เพิ่มตอนรัน)
    done = {(json.loads(l)["cond"], json.loads(l)["seed"]) for l in open(bm.OUT)} if os.path.exists(bm.OUT) else set()
    with open(bm.OUT, "a") as f:
        for c in bm.CONDS:
            for s in range(5):
                if (c, s) not in done:
                    f.write(json.dumps(bm.run_one((c, s))) + chr(10)); f.flush()
    print("done")

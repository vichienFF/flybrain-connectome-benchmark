"""
ขั้น A — รันชุดทดสอบความน่าเชื่อถือ (ดูเกณฑ์ใน BENCHMARK_PLAN.md — ตั้งไว้ก่อนรัน)
บันทึกผลดิบทีละการจำลอง (checkpoint) → results_raw.jsonl · วิเคราะห์ด้วย score.py
ใช้: python benchmark.py [workers=2]   (งดรันช่วง 07:10–07:40)
"""
import os, sys, json, time
from datetime import datetime
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from flybrain import load_shiu_783, Brain  # noqa: E402

VARIANT = os.environ.get("BM_VARIANT", "D0")  # ขั้น D: D0 (เดิม) / D1 / D2 / D3 — ดู D_PLAN.md
CONDSET = os.environ.get("BM_CONDSET", "main")  # main = ชุดทดสอบ 17 ข้อ · probe = ยืนยันใกล้เกณฑ์ (หวาน 50 Hz)
TAG = VARIANT + ("" if os.environ.get("BM_CONDSET", "main") == "main" else "_" + os.environ.get("BM_CONDSET", "main")) + ("" if os.environ.get("BM_SEEDS", "0,5") == "0,5" else "_s" + os.environ["BM_SEEDS"].replace(",", "-"))
OUT = os.path.join(HERE, "results_raw.jsonl" if TAG == "D0" else f"results_raw_{TAG}.jsonl")
SEEDS = range(*map(int, os.environ.get("BM_SEEDS", "0,5").split(",")))  # รอบยืนยัน D2 ใช้ 5,10
T_MS = 1000.0
B = {}


def setup():
    W, ids, ann = load_shiu_783()
    if VARIANT in ("D0", "D0U"):
        if VARIANT == "D0U":  # Usnea (CB0008) เป็นเซลล์กระตุ้น (NT ทาย GABA มั่นใจแค่ 49-66%) — สมมติฐานข้อ 3
            a0 = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids)
            Wc = W.tocsc(copy=True)
            for j in np.flatnonzero(a0.cell_type.values == "CB0008"):
                Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]] = np.abs(Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]])
            W = Wc.tocsr()
        b = Brain(W, ids)
    else:
        import calibrate_d
        k = json.load(open(os.path.join(HERE, f"calib_{VARIANT}.json")))["k"]
        b = calibrate_d.make(W, ids, VARIANT, k)
    ann = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids)
    L = json.load(open(os.path.join(HERE, "ids_shiu_figures.json")))
    idset = set(ids.tolist())
    ix = lambda lst: b.ix([x for x in lst if x in idset])
    ct = lambda t: np.flatnonzero(ann.cell_type.values == t)
    g = dict(sugar=ix(L["neu_sugar"]), water=ix(L["neu_water"]), bitter=ix(L["neu_bitter"]), ir94e=ix(L["neu_ir94e"]),
             jon_ce=ix(L["neu_JON_CE"]), jon_f=ix(L["neu_JON_F"]), rattle=ct("CB0499"), usnea=ct("CB0008"),
             fdg=ct("CB0038"), roundup=ct("CB0553"), lc4=ct("LC4"),
             sugar_other=ix(L["neu_sugar_left"]),
             dnge031_L=np.flatnonzero((ann.cell_type.values == "DNge031") & (ann.side.values == "left")),
             dnge031_R=np.flatnonzero((ann.cell_type.values == "DNge031") & (ann.side.values == "right")), lplc2=ct("LPLC2"), orn_dm1=ct("ORN_DM1"))
    # EPG 6 ตัวที่อยู่ใกล้กันที่สุด (ตามตำแหน่งจริง) = "ส่วนหนึ่งของเข็มทิศ"
    epg = ct("EPG"); P = ann[["pos_x", "pos_y", "pos_z"]].to_numpy(float)[epg]
    d = np.linalg.norm(P[:, None] - P[None], axis=2); seed_i = np.nanargmin(np.nansum(np.sort(d, 1)[:, :6], 1))
    g["epg6"] = epg[np.argsort(d[seed_i])[:6]]
    from phase1_v2_ids import MN9 as _MN9  # MN9 ทั้งสองข้าง (FlyWire v783)
    rd = dict(mn9b=b.ix(_MN9), mn9=b.ix([L["id_mn9"][0]]), abn1=b.ix([L["id_aBN1"][0]]), gf=ct("DNp01"), dm1pn=ct("DM1_lPN"),
              upn=np.flatnonzero(ann.cell_sub_class.values == "uniglomerular"), epg=epg)
    B.update(b=b, g=g, rd=rd, n=len(ids))


# เงื่อนไข: ชื่อ → (รายการ (กลุ่ม, Hz, t0, t1), กลุ่มที่ปิด)
CONDS = {
    "C1_none": ([], []),
    "F1_sugar": ([("sugar", 150, 0, T_MS)], []),
    "F2_water": ([("water", 150, 0, T_MS)], []),
    "F3_bitter": ([("bitter", 150, 0, T_MS)], []),
    "F4_sugar_bitter": ([("sugar", 150, 0, T_MS), ("bitter", 150, 0, T_MS)], []),
    "F5_sugar_ir94e": ([("sugar", 150, 0, T_MS), ("ir94e", 150, 0, T_MS)], []),
    **{f"F6_sugar_{hz}": ([("sugar", hz, 0, T_MS)], []) for hz in (25, 50, 100, 200)},
    "F7_sugar_koRattle": ([("sugar", 150, 0, T_MS)], ["rattle"]),
    "F8_sugar_koUsnea": ([("sugar", 150, 0, T_MS)], ["usnea"]),
    "F9_fdg": ([("fdg", 150, 0, T_MS)], []),
    "F10_sugar_koRoundup": ([("sugar", 150, 0, T_MS)], ["roundup"]),  # post-hoc (เปิดเผยใน BENCHMARK_PLAN.md)
    "G1_jonCE": ([("jon_ce", 150, 0, T_MS)], []),
    "G2_jonF": ([("jon_f", 150, 0, T_MS)], []),
    "E1_lc4": ([("lc4", 150, 0, T_MS)], []),
    "E2_lplc2": ([("lplc2", 150, 0, T_MS)], []),
    "O1_ornDM1": ([("orn_dm1", 150, 0, T_MS)], []),
    "N1_epg6_pulse": ([("epg6", 150, 0, 300)], []),
}


PHZ = float(os.environ.get("BM_PROBE_HZ", "50"))  # (ง) จุดทำงานเท่ากัน: ความแรงหวานของ probe
if CONDSET == "rates":  # (ง) เลือกความแรงหวานของ D0U ให้ MN9 ≈ D0 @50 Hz
    CONDS = {f"R_sugar{hz}": ([("sugar", hz, 0, T_MS)], []) for hz in (30, 35, 40, 45)}
if CONDSET == "opm":  # (ง) ทดสอบปิดเซลล์ที่จุดทำงานเท่ากัน (seed ใหม่)
    CONDS = {"O_sugar": ([("sugar", PHZ, 0, T_MS)], []),
             **{f"O_sugar_ko{k}": ([("sugar", PHZ, 0, T_MS)], [k]) for k in ("rattle", "usnea", "roundup")}}
if CONDSET == "sides":  # (ค) DNge031 ซ้าย/ขวา/สองข้าง × หวานสองด้าน
    CONDS = {f"S_{st}_{ko}": ([(st, 150, 0, T_MS)], {"none": [], "L": ["dnge031_L"], "R": ["dnge031_R"], "LR": ["dnge031_L", "dnge031_R"]}[ko])
             for st in ("sugar", "sugar_other") for ko in ("none", "L", "R", "LR")}
if CONDSET == "water":  # รอบ 4: ปิด Usnea ขณะกระตุ้นรสน้ำ
    CONDS = {"W_water150": ([("water", 150, 0, T_MS)], []), "W_water150_kousnea": ([("water", 150, 0, T_MS)], ["usnea"])}
if CONDSET == "probe":  # ยืนยันผลปิดเซลล์ใกล้เกณฑ์การกิน (ดู HANDOFF ข้อ 3)
    CONDS = {"P_sugar50": ([("sugar", 50, 0, T_MS)], []),
             **{f"P_sugar50_ko{k}": ([("sugar", 50, 0, T_MS)], [k]) for k in ("rattle", "usnea", "roundup")}}


def run_one(task):
    name, seed = task
    if not B:
        setup()
    b, g = B["b"], B["g"]
    stims, ko = CONDS[name]
    dt = b.p["dt"]
    def rate_fn(s):
        t = s * dt; idx, r = [], []
        for grp, hz, t0, t1 in stims:
            if t0 <= t < t1:
                idx.append(g[grp]); r.append(np.full(len(g[grp]), float(hz)))
        return (np.concatenate(idx), np.concatenate(r)) if idx else (np.empty(0, np.int64), np.empty(0))
    sil = np.concatenate([g[k] for k in ko]) if ko else None
    b.rng = np.random.default_rng(seed)
    t = time.time()
    c, bins = b.run(T_MS, rate_fn=rate_fn, silence=sil, bin_ms=100.0)  # bins: 10 ช่อง × 100 ms
    rd = B["rd"]
    late = bins[4:].sum(0)  # 400–1000 ms (สำหรับ N1)
    nz = np.flatnonzero(c)
    return dict(cond=name, seed=seed, sec=round(time.time() - t, 1), n_active=int(len(nz)),
                mn9=c[rd["mn9"]].tolist(), mn9b=c[rd["mn9b"]].tolist(), abn1=c[rd["abn1"]].tolist(), gf=c[rd["gf"]].tolist(),
                dm1pn=c[rd["dm1pn"]].tolist(), upn=c[rd["upn"]].tolist(),
                epg=c[rd["epg"]].tolist(), epg_late=late[rd["epg"]].tolist(), epg6=np.isin(rd["epg"], g["epg6"]).tolist())


def main():
    hm = datetime.now().hour * 60 + datetime.now().minute
    if 7 * 60 + 5 <= hm <= 7 * 60 + 45 and "KAGGLE_KERNEL_RUN_TYPE" not in os.environ:  # กฎนี้สำหรับเครื่อง user เท่านั้น
        sys.exit("ช่วงบอทเทรดทำงาน — งดรัน")
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT, encoding="utf-8"):
            try:
                r = json.loads(l); done.add((r["cond"], r["seed"]))
            except Exception:
                pass
    todo = [(c, s) for c in CONDS for s in SEEDS if (c, s) not in done]
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    print(f"todo {len(todo)} sims, workers {workers}", flush=True)
    t0 = time.time()
    with open(OUT, "a", encoding="utf-8") as f, ProcessPoolExecutor(workers) as ex:
        for i, fu in enumerate(as_completed([ex.submit(run_one, t) for t in todo]), 1):
            f.write(json.dumps(fu.result()) + "\n"); f.flush()
            if i % 10 == 0:
                print(f"{i}/{len(todo)} | {(time.time() - t0) / 60:.0f} min", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()

"""
ข้อ 1 (2026-09-26): กลับเครื่องหมายสารสื่อประสาททีละชนิดเซลล์ (ผู้ต้องสงสัย 44 ชนิดจาก nt_audit) แล้ววัดชุดทดสอบการกิน 8 ข้อ
หา "Usnea ตัวต่อไป" — ชนิดเซลล์ที่ถ้ากลับเครื่องหมายแล้วโมเดลตรงแมลงจริงขึ้น (เกณฑ์ใน NT_FLIP_PLAN.md ตั้งก่อนรัน)
flip = ขาออกทั้งหมดของชนิดนั้นเปลี่ยนเครื่องหมาย (กระตุ้น↔ยับยั้ง) · "none" = D0 เดิม · CB0008 (Usnea) = ตัวควบคุมบวก
ใช้: NF_SHARD=i/N NF_WORKERS=4 python nt_flip_screen.py → nt_flip_s<i>of<N>.jsonl (checkpoint ทุกการจำลอง)
"""
import os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from flybrain import load_shiu_783, Brain  # noqa: E402

T_MS = 1000.0
SEEDS = range(4)
PANEL = {  # ชื่อ: (สิ่งเร้า [(กลุ่ม, Hz)], กลุ่มที่ปิด)
    "P1_sugar150": ([("sugar", 150)], []), "P2_water150": ([("water", 150)], []), "P3_bitter150": ([("bitter", 150)], []),
    "P4_sugar_bitter": ([("sugar", 150), ("bitter", 150)], []), "P5_sugar_ir94e": ([("sugar", 150), ("ir94e", 150)], []),
    "P6_sugar50": ([("sugar", 50)], []), "P7_sugar50_koRattle": ([("sugar", 50)], ["rattle"]),
    "P8_sugar50_koUsnea": ([("sugar", 50)], ["usnea"]), "P9_sugar50_koRoundup": ([("sugar", 50)], ["roundup"]),
}
S = {}


def base():
    if "W" not in S:
        W, ids, ann = load_shiu_783()
        a = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids)
        L = json.load(open(os.path.join(HERE, "ids_shiu_figures.json")))
        idset = set(ids.tolist()); idx = {r: i for i, r in enumerate(ids.tolist())}
        ix = lambda lst: np.array([idx[x] for x in lst if x in idset])
        ct = lambda t: np.flatnonzero(a.cell_type.values == t)
        S.update(W=W.tocsc(), ids=ids, ct=a.cell_type.values, mn9=idx[L["id_mn9"][0]],
                 g=dict(sugar=ix(L["neu_sugar"]), water=ix(L["neu_water"]), bitter=ix(L["neu_bitter"]), ir94e=ix(L["neu_ir94e"]),
                        rattle=ct("CB0499"), usnea=ct("CB0008"), roundup=ct("CB0553")))
    return S


def brain_for(flip):
    s = base()
    if S.get("flip") != flip:
        Wc = s["W"].copy()
        if flip != "none":
            for j in np.flatnonzero(s["ct"] == flip):
                Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]] *= -1
        S["brain"] = Brain(Wc.tocsr(), s["ids"]); S["flip"] = flip
    return S["brain"]


def run_one(task):
    flip, cond, seed = task
    b = brain_for(flip); g = S["g"]
    stims, ko = PANEL[cond]
    idx = np.concatenate([g[k] for k, _ in stims]); r = np.concatenate([np.full(len(g[k]), float(hz)) for k, hz in stims])
    b.rng = np.random.default_rng(seed)
    t = time.time()
    c, _ = b.run(T_MS, rate_fn=lambda s: (idx, r), silence=np.concatenate([g[k] for k in ko]) if ko else None)
    return dict(flip=flip, cond=cond, seed=seed, mn9=int(c[S["mn9"]]), n_active=int((c > 0).sum()), sec=round(time.time() - t, 1))


def main():
    si, sn = map(int, os.environ.get("NF_SHARD", "0/1").split("/"))
    workers = int(os.environ.get("NF_WORKERS", "4"))
    deadline = float(os.environ.get("NF_DEADLINE_S", str(10.3 * 3600)))
    out = os.path.join(os.environ.get("NF_OUT", HERE), f"nt_flip_s{si}of{sn}.jsonl")
    flips = ["none"] + json.load(open(os.path.join(HERE, "results", "flip_types.json")))
    mine = flips[si::sn]
    done = set()
    if os.path.exists(out):
        for l in open(out, encoding="utf-8"):
            try:
                r = json.loads(l); done.add((r["flip"], r["cond"], r["seed"]))
            except Exception:
                pass
    todo = [(f, c, s) for f in mine for c in PANEL for s in SEEDS if (f, c, s) not in done]
    print(f"shard {si}/{sn}: flips {len(mine)} todo {len(todo)} workers {workers}", flush=True)
    # แบ่งงานเป็นก้อนตามชนิดเซลล์ต่อ worker → แต่ละ worker สร้าง W ที่กลับเครื่องหมายน้อยครั้ง
    t0 = time.time(); n = 0
    with open(out, "a", encoding="utf-8") as f, ProcessPoolExecutor(workers) as ex:
        it, pend = iter(todo), set()
        while True:
            while len(pend) < workers * 2 and time.time() - t0 < deadline:
                t = next(it, None)
                if t is None:
                    break
                pend.add(ex.submit(run_one, t))
            if not pend:
                break
            fin, pend = wait(pend, return_when=FIRST_COMPLETED)
            for fu in fin:
                f.write(json.dumps(fu.result()) + chr(10)); f.flush(); n += 1
                if n % 36 == 0:
                    print(f"{n}/{len(todo)} | {(time.time() - t0) / 60:.0f} min", flush=True)
    print(f"finished {n}/{len(todo)} in {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()

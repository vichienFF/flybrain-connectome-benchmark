"""
knockout_screen.py — กวาดปิดเซลล์ทีละตัว (in-silico knockout screen) หาเซลล์ที่คุม "กิน/ไม่กิน"
วัดผลที่ MN9 (สั่งยื่นงวง, ซ้าย+ขวา) · โมเดล: flybrain.py (Shiu et al. 2024, FlyWire v783) · แยกจากบอททุกตัว

หลักคิด (ทำให้ทั้งสมอง 138,639 เซลล์ เหลืองานไม่กี่ร้อย–พันครั้งแบบไม่ประมาณ):
- ปิดเซลล์ที่ "ไม่ยิงเลย" ใน baseline → ผลเหมือน baseline ทุกบิต (ตัวสุ่มใช้เฉพาะ Poisson อินพุต)
  ดังนั้นกวาดเฉพาะเซลล์ที่ยิงใน baseline (รวมทุก seed ของ baseline) ก็ครอบคลุมทั้งสมอง
- common random numbers: KO กับ baseline ใช้ seed เดียวกัน → เทียบแบบจับคู่ ลดสัญญาณรบกวน
- รอบ 1 (pass 1): ทุกเซลล์ 1 seed · รอบ 2 (pass 2): ผู้ต้องสงสัยอันดับต้น × หลาย seed (ยืนยันสถิติ)

ตั้งค่าผ่าน env: KS_COND (sugar|sugar_bitter) KS_PASS (1|2) KS_SHARD "i/N" KS_WORKERS KS_DEADLINE_S
                 KS_OUT (โฟลเดอร์ผล) KS_SEEDS (pass2, เช่น "0-9") KS_CANDIDATES (pass2, ไฟล์ json รายการ index)
ผล: KS_OUT/ks_<cond>_p<pass>_s<i>of<N>.jsonl — 1 บรรทัด/การจำลอง (checkpoint, รันซ้ำจะข้ามที่เสร็จแล้ว)
"""
import os, sys, json, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("KS_ROOT", os.path.dirname(HERE)))
from flybrain import load_shiu_783, Brain  # noqa: E402
from phase1_v2_ids import SUGAR, MN9  # noqa: E402

T_MS = float(os.environ.get("KS_TMS", "1000"))
SUGAR_HZ = 150.0
BITTER_HZ = float(os.environ.get("KS_BITTER_HZ", "100"))
B = {}


def setup():
    W, ids, ann = load_shiu_783()
    b = Brain(W, ids)
    idset = set(ids.tolist())
    sugar = b.ix([x for x in SUGAR if x in idset])
    bitter = b.ix([x for x in ann.root_id[(ann.cell_sub_class == "bitter") & (ann.side == "right")] if x in idset])
    B.update(b=b, ids=ids, ann=ann, mn9=b.ix(MN9), sugar=sugar, bitter=bitter)


def stim_of(cond):
    if cond == "sugar":
        return B["sugar"], np.full(len(B["sugar"]), SUGAR_HZ)
    if cond == "sugar_bitter":
        st = np.concatenate([B["sugar"], B["bitter"]])
        return st, np.r_[np.full(len(B["sugar"]), SUGAR_HZ), np.full(len(B["bitter"]), BITTER_HZ)]
    raise ValueError(cond)


def run_one(task):
    """task = (cond, ko_index or -1, seed) → dict ผล (นับ spike เฉพาะเซลล์ที่ยิง)"""
    cond, ko, seed = task
    b = B["b"]
    st, rates = stim_of(cond)
    b.rng = np.random.default_rng(seed)
    t = time.time()
    c, _ = b.run(T_MS, rate_fn=lambda s: (st, rates), silence=None if ko < 0 else np.array([ko]))
    nz = np.flatnonzero(c)
    return dict(cond=cond, ko=int(ko), seed=int(seed), mn9=c[B["mn9"]].tolist(), n_active=int(len(nz)),
                idx=nz.tolist(), cnt=c[nz].tolist(), sec=round(time.time() - t, 1))


def _init():  # worker: fork ใช้หน่วยความจำร่วม (Linux) · Windows โหลดใหม่ในแต่ละ worker
    if not B:
        setup()


def load_done(path):
    done = set()
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            try:
                r = json.loads(line); done.add((r["cond"], r["ko"], r["seed"]))
            except Exception:
                pass  # บรรทัดสุดท้ายที่เขียนไม่ครบ (ถูกตัดกลางคัน) → ทำใหม่
    return done


def parse_seeds(s):
    a, _, z = s.partition("-")
    return list(range(int(a), int(z or a) + 1))


def main():
    t0 = time.time()
    cond = os.environ.get("KS_COND", "sugar")
    ps = int(os.environ.get("KS_PASS", "1"))
    si, sn = map(int, os.environ.get("KS_SHARD", "0/1").split("/"))
    workers = int(os.environ.get("KS_WORKERS", str(os.cpu_count() or 1)))
    deadline = float(os.environ.get("KS_DEADLINE_S", str(10.5 * 3600)))
    out = os.environ.get("KS_OUT", HERE)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"ks_{cond}_p{ps}_s{si}of{sn}.jsonl")
    setup()
    base_seeds = [0, 1, 2]
    if ps == 1:
        # baseline 3 seed → รายชื่อเซลล์ที่ยิง (union) = รายการปิดเซลล์ทั้งสมอง
        bpath = os.path.join(out, f"ks_{cond}_baseline.json")
        if os.path.exists(bpath):
            base = json.load(open(bpath, encoding="utf-8"))
        else:
            base = [run_one((cond, -1, s)) for s in base_seeds]
            json.dump(base, open(bpath, "w", encoding="utf-8"))
        ko_list = sorted(set().union(*[set(r["idx"]) for r in base]))
        lim = int(os.environ.get("KS_LIMIT", "0"))  # ใช้ทดสอบเล็กเท่านั้น
        if lim:
            ko_list = ko_list[:lim]
        tasks = [(cond, k, 0) for k in ko_list]
        print(f"[{cond}] baseline MN9={[r['mn9'] for r in base]} active={[r['n_active'] for r in base]} "
              f"-> KO list {len(ko_list)} cells", flush=True)
    else:
        seeds = parse_seeds(os.environ.get("KS_SEEDS", "0-9"))
        cands = json.load(open(os.environ["KS_CANDIDATES"], encoding="utf-8"))
        tasks = [(cond, -1, s) for s in seeds] + [(cond, k, s) for k in cands for s in seeds]
    tasks = tasks[si::sn]
    done = load_done(path)
    todo = [t for t in tasks if t not in done]
    print(f"shard {si}/{sn}: {len(tasks)} tasks, done {len(done)}, todo {len(todo)}, workers {workers}", flush=True)
    n_ok, stopped = 0, False
    with open(path, "a", encoding="utf-8") as f, ProcessPoolExecutor(workers, initializer=_init) as ex:
        it, pend = iter(todo), set()
        while True:
            while not stopped and len(pend) < workers:
                if time.time() - t0 > deadline:
                    stopped = True
                    print(f"deadline {deadline / 3600:.1f} h reached - no new tasks", flush=True)
                    break
                t = next(it, None)
                if t is None:
                    break
                pend.add(ex.submit(run_one, t))
            if not pend:
                break
            fin, pend = wait(pend, return_when=FIRST_COMPLETED)
            for fu in fin:
                f.write(json.dumps(fu.result()) + "\n"); f.flush()
                n_ok += 1
                if n_ok % 20 == 0:
                    el = time.time() - t0
                    print(f"{n_ok}/{len(todo)} | {el / 60:.0f} min | ETA {el / n_ok * (len(todo) - n_ok) / 60:.0f} min", flush=True)
    left = len(todo) - n_ok
    print(f"finished {n_ok} tasks in {(time.time() - t0) / 60:.1f} min | remaining {left}", flush=True)
    open(os.path.join(out, f"DONE_{cond}_p{ps}_s{si}of{sn}.txt" if left == 0 else f"PARTIAL_{cond}_p{ps}_s{si}of{sn}.txt"),
         "w").write(f"{n_ok} ok, {left} left\n")


if __name__ == "__main__":
    main()

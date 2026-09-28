"""
flybrain.py — ตัวจำลองสมองแมลงหวี่ทั้งหัว (Leaky Integrate-and-Fire) จาก FlyWire v783
อ้างอิงพารามิเตอร์: Shiu et al., Nature 2024 (whole-brain LIF model)
ข้อมูล: D:\\Research_FlyBrain\\flywire\\ (CC-BY 4.0, FlyWire Consortium) — แยกจากบอททุกตัว

สมมติฐานของโมเดล (ต้องรายงานทุกครั้ง):
- acetylcholine = กระตุ้น (+), GABA/glutamate = ยับยั้ง (−)
- dopamine/serotonin/octopamine ถือเป็นกระตุ้น (+) แบบเร็ว (ของจริงเป็น neuromodulator — ข้อจำกัด)
- น้ำหนัก = จำนวน synapse × w_syn, ไม่มี plasticity
- v2 (2026-09-23): แก้ให้ตรง Brian2 ต้นฉบับ — reset g=0 ตอนยิง, v/g หยุดระหว่างช่วงพัก, Poisson เข้า v, เซลล์กระตุ้นไม่มีช่วงพัก
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
FW = os.path.join(HERE, "flywire")
MIN_SYN = int(os.environ.get("FLY_MIN_SYN", "5"))  # มาตรฐาน FlyWire: ตัดการเชื่อมต่อ < 5 synapse
CACHE = os.path.join(FW, f"W_783_min{MIN_SYN}.npz")

# Shiu et al. 2024
P = dict(v0=-52.0, v_reset=-52.0, v_th=-45.0, tau_m=20.0, tau_syn=5.0,
         t_ref=2.2, t_dly=1.8, w_syn=0.275, r_poi=150.0, f_poi=250.0, dt=0.1)
EXC = {"acetylcholine": 1, "gaba": -1, "glutamate": -1,
       "dopamine": 1, "serotonin": 1, "octopamine": 1}


def load_annotations():
    return pd.read_csv(os.path.join(FW, "neuron_annotations.tsv"), sep="\t", low_memory=False)


def load_connectome():
    """คืน (W_post_by_pre เป็น CSR ลายเซ็นแล้ว, root_ids, ann). cache ไว้ใน flywire/"""
    ann = load_annotations()
    ids = np.load(os.path.join(FW, f"root_ids_min{MIN_SYN}.npy")) if os.path.exists(CACHE) else None
    if ids is not None:
        W = sp.load_npz(CACHE).tocsr()
    else:
        c = pd.read_feather(os.path.join(FW, "proofread_connections_783.feather"))
        c = c.groupby(["pre_pt_root_id", "post_pt_root_id"], as_index=False)["syn_count"].sum()
        c = c[c.syn_count >= MIN_SYN]
        ids = np.union1d(ann.root_id.values, np.union1d(c.pre_pt_root_id, c.post_pt_root_id))
        idx = pd.Series(np.arange(len(ids)), index=ids)
        nt = ann.set_index("root_id").top_nt.map(EXC).reindex(ids).fillna(1).values
        pre, post = idx[c.pre_pt_root_id].values, idx[c.post_pt_root_id].values
        w = c.syn_count.values * nt[pre]
        W = sp.csr_matrix((w.astype(np.float32), (post, pre)), shape=(len(ids), len(ids)))
        sp.save_npz(CACHE, W)
        np.save(os.path.join(FW, f"root_ids_min{MIN_SYN}.npy"), ids)
    return W, ids, ann


def load_shiu_783():
    """โหลด Connectivity_783 ของ repo ทางการ (MIT) — ใช้ทุกเส้น ไม่ตัด threshold, เครื่องหมายตามต้นฉบับ"""
    base = os.path.join(HERE, "shiu_model")
    comp = pd.read_csv(os.path.join(base, "Completeness_783.csv"), index_col=0)
    con = pd.read_parquet(os.path.join(base, "Connectivity_783.parquet"))
    n = len(comp)
    W = sp.csr_matrix((con["Excitatory x Connectivity"].values.astype(np.float32),
                       (con.Postsynaptic_Index.values, con.Presynaptic_Index.values)), shape=(n, n))
    return W, comp.index.values.astype(np.int64), load_annotations()


class Brain:
    def __init__(self, W, ids, p=P, seed=0):
        self.W, self.ids, self.p = W, ids, p
        self.WT = W.T.tocsr()  # row j = ขาออกของ neuron j
        self.idx = pd.Series(np.arange(len(ids)), index=ids)
        self.rng = np.random.default_rng(seed)

    def ix(self, root_ids):
        return self.idx[np.asarray(root_ids)].values

    def run(self, t_ms=1000.0, poisson=None, rate_fn=None, record=None, silence=None, bin_ms=None, progress=None):
        """
        poisson: {index_array: rate_Hz} อินพุต Poisson คงที่
        rate_fn: f(step)->(index_array, rates_Hz) อินพุตแปรตามเวลา (สำหรับ task)
        silence: index_array ที่ถูกปิด (in-silico knockout)
        record: index_array ที่บันทึก spike ทุก step (คืน matrix bool)
        คืน spike counts ต่อ neuron (+ raster ของ record)
        """
        p, n = self.p, self.W.shape[0]
        dt = p["dt"]
        steps = int(t_ms / dt)
        dly = int(round(p["t_dly"] / dt))
        ref = int(round(p["t_ref"] / dt))
        v = np.full(n, p["v0"])
        g = np.zeros(n)
        ref_left = np.zeros(n, dtype=np.int32)
        counts = np.zeros(n, dtype=np.int32)
        buf = [np.empty(0, dtype=np.int64)] * dly  # ส่งถึงปลายทางหลัง dly steps พอดี
        a_m, a_s = dt / p["tau_m"], np.exp(-dt / p["tau_syn"])
        e_m = np.exp(-dt / p["tau_m"])
        k_vg = p["tau_syn"] / (p["tau_syn"] - p["tau_m"]) * (a_s - e_m)
        poi_w = p["w_syn"] * p["f_poi"]
        mask_ok = np.ones(n, bool)
        if silence is not None:
            mask_ok[silence] = False
        raster = np.zeros((steps, len(record)), bool) if record is not None else None
        # bin_ms: นับ spike เป็นช่วงเวลา (สำหรับภาพเคลื่อนไหว) → คืนเป็น raster แทน
        bsteps = int(round(bin_ms / dt)) if bin_ms else 0
        bins = np.zeros((int(np.ceil(steps / bsteps)), n), np.uint16) if bsteps else None
        no_ref = np.zeros(n, bool)  # เซลล์ที่รับ Poisson ไม่มีช่วงพัก (ตามต้นฉบับ)
        for s in range(steps):
            srcs = []
            if poisson:
                srcs += list(poisson.items())
            if rate_fn is not None:
                srcs.append(rate_fn(s))
            active = ref_left == 0
            # ODE (unless refractory): v, g หยุดนิ่งระหว่างช่วงพัก
            # exact linear integration (เหมือน Brian2 method='linear')
            g_a = g[active]
            v[active] = p["v0"] + (v[active] - p["v0"]) * e_m + g_a * k_vg
            g[active] = g_a * a_s
            ref_left[~active] -= 1
            # Poisson → v โดยตรง
            for pidx, r in srcs:
                no_ref[pidx] = True
                hit = self.rng.random(len(pidx)) < np.asarray(r) * dt / 1000.0
                if hit.any():
                    v[pidx[hit]] += poi_w
            spk = np.flatnonzero((v > p["v_th"]) & mask_ok & ((ref_left == 0) | no_ref))
            v[spk] = p["v_reset"]
            g[spk] = 0.0  # reset g = 0 (จุดที่เวอร์ชันแรกขาด)
            ref_left[spk] = np.where(no_ref[spk], 0, ref)
            counts[spk] += 1
            if bins is not None and spk.size:
                bins[s // bsteps, spk] += 1
            if progress is not None and s % 500 == 0:
                progress(s / steps)
            # synaptic events ที่ครบ delay → g (เกิดได้แม้อยู่ในช่วงพัก)
            arr = buf[s % dly]
            if arr.size:
                g += p["w_syn"] * np.asarray(self.WT[arr].sum(axis=0)).ravel()
            buf[s % dly] = spk
            if raster is not None:
                raster[s] = np.isin(record, spk)
        return counts, (bins if bins is not None else raster)

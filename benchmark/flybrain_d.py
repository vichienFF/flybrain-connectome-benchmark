"""
flybrain_d.py — ขั้น D: โมเดลอัปเกรด (ไม่แตะ flybrain.py ต้นฉบับ ซึ่ง Simulator/knockout screen ใช้อยู่)
เพิ่ม 2 กลไกที่มีหลักฐานทางชีววิทยา (ปิดได้ — เมื่อปิดทั้งคู่ ผลต้องเหมือน Brain เดิมทุกบิต):
1) STD  synapse ล้าระยะสั้น (Tsodyks–Markram แบบย่อ ต่อเซลล์ต้นทาง): ทุกครั้งที่ยิง ใช้ทรัพยากร x ไป U·x,
        แรงส่ง = x/1 ของครั้งแรก (normalize ให้ spike แรกเท่าเดิม), x ฟื้นกลับหา 1 ด้วยเวลา tau_rec
        หลักฐาน: Kazama & Wilson 2008 Neuron (ORN→PN depression)
2) SFA  เซลล์ล้า (spike-frequency adaptation แบบ threshold): ทุกครั้งที่ยิง เกณฑ์ยิงสูงขึ้น dth mV แล้วคืนตัวด้วย tau_th
"""
import numpy as np
from flybrain import Brain


class BrainD(Brain):
    def __init__(self, W, ids, std=None, sfa=None, **kw):
        super().__init__(W, ids, **kw)
        self.std, self.sfa = std, sfa  # std=dict(U, tau_rec[, mask]) · sfa=dict(dth, tau_th)
        # mask (bool ต่อเซลล์ต้นทาง): ใช้ synapse ล้าเฉพาะเซลล์ที่มีหลักฐาน (เช่น ORN→PN) · ไม่ระบุ = ทุกเซลล์

    def run(self, t_ms=1000.0, poisson=None, rate_fn=None, record=None, silence=None, bin_ms=None, progress=None):
        p, n = self.p, self.W.shape[0]
        dt = p["dt"]; steps = int(t_ms / dt)
        dly = int(round(p["t_dly"] / dt)); ref = int(round(p["t_ref"] / dt))
        v = np.full(n, p["v0"]); g = np.zeros(n)
        ref_left = np.zeros(n, dtype=np.int32); counts = np.zeros(n, dtype=np.int32)
        buf = [np.empty(0, dtype=np.int64)] * dly
        effbuf = [np.empty(0)] * dly  # แรงส่งของแต่ละ spike ณ เวลาที่ยิง (STD)
        a_s = np.exp(-dt / p["tau_syn"]); e_m = np.exp(-dt / p["tau_m"])
        k_vg = p["tau_syn"] / (p["tau_syn"] - p["tau_m"]) * (a_s - e_m)
        poi_w = p["w_syn"] * p["f_poi"]
        mask_ok = np.ones(n, bool)
        if silence is not None:
            mask_ok[silence] = False
        bsteps = int(round(bin_ms / dt)) if bin_ms else 0
        bins = np.zeros((int(np.ceil(steps / bsteps)), n), np.uint16) if bsteps else None
        no_ref = np.zeros(n, bool)
        x = np.ones(n) if self.std else None
        if self.std:
            U, rec = self.std["U"], dt / self.std["tau_rec"]
        th = np.zeros(n) if self.sfa else None
        if self.sfa:
            dth, e_th = self.sfa["dth"], np.exp(-dt / self.sfa["tau_th"])
        for s in range(steps):
            srcs = list(poisson.items()) if poisson else []
            if rate_fn is not None:
                srcs.append(rate_fn(s))
            active = ref_left == 0
            g_a = g[active]
            v[active] = p["v0"] + (v[active] - p["v0"]) * e_m + g_a * k_vg
            g[active] = g_a * a_s
            ref_left[~active] -= 1
            for pidx, r in srcs:
                no_ref[pidx] = True
                hit = self.rng.random(len(pidx)) < np.asarray(r) * dt / 1000.0
                if hit.any():
                    v[pidx[hit]] += poi_w
            vth = p["v_th"] + th if self.sfa else p["v_th"]
            spk = np.flatnonzero((v > vth) & mask_ok & ((ref_left == 0) | no_ref))
            v[spk] = p["v_reset"]; g[spk] = 0.0
            ref_left[spk] = np.where(no_ref[spk], 0, ref)
            counts[spk] += 1
            if self.std:
                x += (1.0 - x) * rec
                eff = x[spk].copy()          # แรงส่งของ spike นี้ (เทียบกับ spike แรก = 1)
                if "mask" in self.std:
                    m = self.std["mask"][spk]; eff[~m] = 1.0
                    x[spk[m]] -= U * x[spk[m]]
                else:
                    x[spk] -= U * x[spk]
            if self.sfa:
                th *= e_th
                th[spk] += dth
            if bins is not None and spk.size:
                bins[s // bsteps, spk] += 1
            arr = buf[s % dly]
            if arr.size:
                if self.std:
                    g += p["w_syn"] * np.asarray(self.WT[arr].T @ effbuf[s % dly]).ravel()
                else:
                    g += p["w_syn"] * np.asarray(self.WT[arr].sum(axis=0)).ravel()
            buf[s % dly] = spk
            if self.std:
                effbuf[s % dly] = eff
        return counts, bins

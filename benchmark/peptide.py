"""
peptide.py (2026-10-05) - a rule that MIMICS a slow neuropeptide co-transmitter (no chemistry is simulated).
Motivation: Amon (neuropeptide-processing enzyme) knockdown in Usnea reproduces the Usnea silencing phenotype
(G. Sterne, personal communication, 2026-10-01). Usnea keeps its GABA (fast, inhibitory) synapses unchanged.
Rule: every spike of a source neuron adds 1 to a shared "peptide level" c, which decays with time constant tau_ms.
Each target neuron j receives a slow depolarising drive  dv_j = gain * c * w_j * dt  (mV), w_j = synapses from the
source set to j divided by the maximum over targets (targets = direct postsynaptic partners of the source set).
With gain = 0 the model is bit-identical to Brain (flybrain.py).
"""
import numpy as np
from flybrain import Brain


class BrainP(Brain):
    def __init__(self, W, ids, src=None, gain=0.0, tau_ms=500.0, **kw):
        super().__init__(W, ids, **kw)
        self.src = np.asarray(src if src is not None else [], dtype=np.int64)
        self.gain, self.tau = float(gain), float(tau_ms)
        if len(self.src):
            syn = np.asarray(abs(self.WT[self.src]).sum(axis=0)).ravel()  # synapse counts (sign removed)
            syn[self.src] = 0.0
            self.tgt = np.flatnonzero(syn > 0)
            self.w_t = syn[self.tgt] / syn[self.tgt].max()
        else:
            self.tgt = np.empty(0, np.int64); self.w_t = np.empty(0)

    def run(self, t_ms=1000.0, poisson=None, rate_fn=None, record=None, silence=None, bin_ms=None, progress=None):
        p, n = self.p, self.W.shape[0]
        dt = p["dt"]; steps = int(t_ms / dt)
        dly = int(round(p["t_dly"] / dt)); ref = int(round(p["t_ref"] / dt))
        v = np.full(n, p["v0"]); g = np.zeros(n)
        ref_left = np.zeros(n, dtype=np.int32); counts = np.zeros(n, dtype=np.int32)
        buf = [np.empty(0, dtype=np.int64)] * dly
        a_s = np.exp(-dt / p["tau_syn"]); e_m = np.exp(-dt / p["tau_m"])
        k_vg = p["tau_syn"] / (p["tau_syn"] - p["tau_m"]) * (a_s - e_m)
        poi_w = p["w_syn"] * p["f_poi"]
        mask_ok = np.ones(n, bool)
        if silence is not None:
            mask_ok[silence] = False
        bsteps = int(round(bin_ms / dt)) if bin_ms else 0
        bins = np.zeros((int(np.ceil(steps / bsteps)), n), np.uint16) if bsteps else None
        raster = np.zeros((steps, len(record)), bool) if record is not None else None
        no_ref = np.zeros(n, bool)
        use_pep = self.gain > 0 and len(self.tgt) > 0
        c, e_c = 0.0, np.exp(-dt / self.tau)
        is_src = np.zeros(n, bool); is_src[self.src] = True
        for s in range(steps):
            srcs = list(poisson.items()) if poisson else []
            if rate_fn is not None:
                srcs.append(rate_fn(s))
            active = ref_left == 0
            g_a = g[active]
            v[active] = p["v0"] + (v[active] - p["v0"]) * e_m + g_a * k_vg
            g[active] = g_a * a_s
            if use_pep and c > 0:
                c *= e_c
                act_t = active[self.tgt]
                v[self.tgt[act_t]] += self.gain * c * self.w_t[act_t] * dt
            ref_left[~active] -= 1
            for pidx, r in srcs:
                no_ref[pidx] = True
                hit = self.rng.random(len(pidx)) < np.asarray(r) * dt / 1000.0
                if hit.any():
                    v[pidx[hit]] += poi_w
            spk = np.flatnonzero((v > p["v_th"]) & mask_ok & ((ref_left == 0) | no_ref))
            v[spk] = p["v_reset"]; g[spk] = 0.0
            ref_left[spk] = np.where(no_ref[spk], 0, ref)
            counts[spk] += 1
            if use_pep and spk.size:
                c += float(is_src[spk].sum())
            if bins is not None and spk.size:
                bins[s // bsteps, spk] += 1
            arr = buf[s % dly]
            if arr.size:
                g += p["w_syn"] * np.asarray(self.WT[arr].sum(axis=0)).ravel()
            buf[s % dly] = spk
            if raster is not None:
                raster[s] = np.isin(record, spk)
        return counts, (bins if bins is not None else raster)

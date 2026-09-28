r"""
สร้างรูปสำหรับ preprint จากไฟล์ผลจริงใน D:\Research_FlyBrain (ไม่มีตัวเลขพิมพ์มือ) → preprint/figures/*.png (300 dpi) + .pdf
"""
import os, json, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

R = r"D:\Research_FlyBrain"; B = os.path.join(R, "benchmark"); K = os.path.join(R, "knockout_screen", "results"); M = os.path.join(R, "malecns")
OUT = os.path.join(R, "preprint", "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
C_OK, C_BAD, C_GRAY, C_ACC = "#1D9E75", "#D85A30", "#888780", "#534AB7"
jl = lambda p: [json.loads(l) for l in open(p, encoding="utf-8")]


def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), dpi=300); fig.savefig(os.path.join(OUT, name + ".pdf")); plt.close(fig)
    print("saved", name)


def paired(path, base, ko, key=lambda r: r["mn9"][0]):
    rr = jl(path); b = {r["seed"]: key(r) for r in rr if r["cond"] == base}; k = {r["seed"]: key(r) for r in rr if r["cond"] == ko}
    s = sorted(set(b) & set(k)); d = np.array([k[i] - b[i] for i in s], float); bm = np.mean([b[i] for i in s])
    return 100 * d.mean() / bm, 100 * d.std(ddof=1) / np.sqrt(len(d)) / bm, stats.ttest_1samp(d, 0).pvalue, len(s)


# ---------- Fig 1: benchmark pass matrix ----------
def fig1():
    order = ["C1", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "G1", "G2", "E1", "E2", "O1", "N1"]
    names = {"D0": "D0 original (Shiu 2024)", "D1": "D1 +STD (global)", "D2": "D2 +SFA", "D3": "D3 +STD+SFA", "D4": "D4 +SFA +ORN-STD", "D0U": "D0U Usnea excitatory"}
    rows = []
    for v in names:
        f = os.path.join(B, "scores.json" if v == "D0" else f"scores_{v}.json")
        s = {x["test"]: x["passed"] for x in json.load(open(f, encoding="utf-8"))}
        rows.append([1 if s.get(t) else 0 for t in order])
    A = np.array(rows)
    fig, ax = plt.subplots(figsize=(7, 2.5))
    ax.imshow(A, cmap=matplotlib.colors.ListedColormap([C_BAD, C_OK]), aspect="auto")
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.set_yticks(range(len(names)))
    ax.set_yticklabels([f"{n}  ({A[i].sum()}/16)" for i, n in enumerate(names.values())])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title("Pre-registered benchmark vs. fly experiments (green = pass)", loc="left")
    save(fig, "fig1_benchmark")


# ---------- Fig 2: knockout screen confirmation ----------
def fig2():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6))
    for ax, cond, title in zip(axs, ("sugar", "sugar_bitter"), ("Sugar", "Sugar + bitter")):
        d = pd.read_csv(os.path.join(K, f"p2_rank_{cond}.csv"), dtype={"root_id": str})
        p1 = pd.read_csv(os.path.join(K, f"p1_rank_{cond}.csv"), dtype={"root_id": str})[["ko", "side", "lit_name"]]
        d = d.merge(p1, on="ko", how="left").sort_values("p2_d_pct", key=lambda x: -x.abs()).head(16).iloc[::-1]
        lab = [f"{t} ({s[0].upper() if isinstance(s, str) else '?'})" + (f" [{n}]" if isinstance(n, str) else "") for t, s, n in zip(d.cell_type, d.side, d.lit_name)]
        col = [C_OK if s else C_GRAY for s in d.p2_sig]
        ax.barh(range(len(d)), d.p2_d_pct, xerr=[d.p2_d_pct - d.p2_ci_lo, d.p2_ci_hi - d.p2_d_pct], color=col, ecolor="#444", capsize=1.5)
        ax.set_yticks(range(len(d))); ax.set_yticklabels(lab, fontsize=6.5); ax.axvline(0, color="#444", lw=.6)
        ax.set_title(title, loc="left")
    fig.supxlabel("ΔMN9 after single-cell knockout (%; mean ± 95% CI, 10 paired seeds; green = FDR q<0.05 and |Δ|≥5%)", fontsize=7.5)
    save(fig, "fig2_knockout_screen")


# ---------- Fig 3: why knockouts fail ----------
def fig3():
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6))
    ax = axs[0]; labs = ["Rattle", "Usnea", "Roundup"]; hi, lo = [], []
    for k, c in zip(("rattle", "usnea", "roundup"), ("F7_sugar_koRattle", "F8_sugar_koUsnea", "F10_sugar_koRoundup")):
        hi.append(paired(os.path.join(B, "results_raw.jsonl"), "F1_sugar", c)[0])
        lo.append(paired(os.path.join(B, "results_raw_D0_opm_s20-40.jsonl"), "O_sugar", f"O_sugar_ko{k}")[0])
    x = np.arange(3); ax.bar(x - .2, hi, .4, color=C_GRAY, label="sugar 150 Hz (5 seeds)"); ax.bar(x + .2, lo, .4, color=C_ACC, label="sugar 50 Hz (fresh seeds 20–39)")
    ax.axhline(0, color="#444", lw=.6); ax.set_xticks(x); ax.set_xticklabels(labs); ax.set_ylabel("ΔMN9 (%)"); ax.legend(fontsize=6, frameon=False, loc="lower left")
    ax.set_title("A  Operating point, model D0\n(fly: ↓, ↓, no change; 20 seeds at 50 Hz)", loc="left", fontsize=7.5)
    ax = axs[1]
    comp = pd.read_csv(os.path.join(R, "shiu_model", "Completeness_783.csv"), index_col=0); ids = comp.index.values.astype(np.int64)
    con = pd.read_parquet(os.path.join(R, "shiu_model", "Connectivity_783.parquet"))
    ann = pd.read_csv(os.path.join(R, "flywire", "neuron_annotations.tsv"), sep="\t", low_memory=False).drop_duplicates("root_id").set_index("root_id").reindex(ids)
    base = json.load(open(os.path.join(K, "p1", "ks_sugar_baseline.json")))[0]; rate = dict(zip(base["idx"], base["cnt"]))
    mn9 = int(np.flatnonzero(ids == 720575940660219265)[0]); inp = con[con.Postsynaptic_Index == mn9].copy()
    inp["drive"] = inp["Excitatory x Connectivity"] * inp.Presynaptic_Index.map(rate).fillna(0); e = inp[inp.drive > 0].sort_values("drive", ascending=False)
    share = 100 * e.drive / e.drive.sum(); names = ann.cell_type.values[e.Presynaptic_Index]
    top = list(zip(names[:3], share[:3])) + [("other\n(n=%d)" % (len(e) - 3), share[3:].sum())]
    ax.bar([t for t, _ in top], [s for _, s in top], color=[C_BAD, C_GRAY, C_GRAY, "#ccc"]); ax.set_ylabel("excitatory drive to MN9 (%)")
    ax.set_title("B  Bottleneck\n(CB0553 = Roundup)", loc="left", fontsize=7.5); ax.tick_params(axis="x", labelsize=6)
    ax = axs[2]; a = pd.read_csv(os.path.join(B, "results", "nt_audit.csv"))
    ax.scatter(a.top_nt_conf, a.male_conf, s=6, c=np.where(a.suspect, C_BAD, C_GRAY), alpha=.7)
    u = a[a.cell_type == "CB0008"]; ax.scatter(u.top_nt_conf, u.male_conf, s=30, facecolors="none", edgecolors=C_ACC, lw=1.2)
    ax.annotate("Usnea", (u.top_nt_conf.iloc[0], u.male_conf.iloc[0]), xytext=(.3, .45), fontsize=7, color=C_ACC, arrowprops=dict(arrowstyle="-", color=C_ACC, lw=.6))
    ax.axvline(.7, ls=":", color="#444", lw=.6); ax.axhline(.7, ls=":", color="#444", lw=.6)
    ax.set_xlabel("FlyWire NT confidence"); ax.set_ylabel("Male CNS NT confidence"); ax.set_title(f"C  NT uncertainty\n({int(a.suspect.sum())}/{len(a)} suspect)", loc="left", fontsize=7.5)
    save(fig, "fig3_why_knockouts_fail")


# ---------- Fig 4: Usnea ----------
def fig4():
    fig, axs = plt.subplots(1, 4, figsize=(7.6, 2.4))
    ax = axs[0]
    w0 = np.mean([r["mn9"][0] for r in jl(os.path.join(B, "results_raw.jsonl")) if r["cond"] == "F2_water"])
    w1 = np.mean([r["mn9"][0] for r in jl(os.path.join(B, "results_raw_D0U.jsonl")) if r["cond"] == "F2_water"])
    ax.bar(["Usnea\nGABA", "Usnea\nexcitatory"], [w0, w1], color=[C_GRAY, C_OK]); ax.axhline(20, ls=":", color="#444", lw=.6)
    ax.set_ylabel("MN9 (Hz)"); ax.set_title("A  Water → MN9\n(out-of-sample)", loc="left", fontsize=7.5)
    for ax, (path_a, path_b, base, ko, t) in zip(axs[1:3], [
            ("results_raw_D0_probe_s0-20.jsonl", "results_raw_D0U_probe_s0-20.jsonl", "P_sugar50", "P_sugar50_kousnea", "B  Knock out Usnea\n(sugar 50 Hz, 20 seeds)"),
            ("results_raw_D0_water_s0-10.jsonl", "results_raw_D0U_water_s0-10.jsonl", "W_water150", "W_water150_kousnea", "C  Knock out Usnea\n(water, 10 seeds)")]):
        v = [paired(os.path.join(B, p), base, ko) for p in (path_a, path_b)]
        ax.bar(["GABA", "excitatory"], [x[0] for x in v], yerr=[x[1] for x in v], color=[C_BAD, C_OK], capsize=2)
        ax.axhline(0, color="#444", lw=.6); ax.set_ylabel("ΔMN9 (%)"); ax.set_title(t, loc="left", fontsize=7.5)
    ax = axs[3]; t = pd.read_csv(os.path.join(B, "results", "ntflip_scores.csv"), index_col=0).drop(index="none")
    better = t.score > 5
    ax.bar(["better", "same", "worse"], [int(better.sum()), int((t.score == 5).sum()), int((t.score < 5).sum())], color=[C_OK, C_GRAY, C_BAD])
    ax.set_title("D  Flip 44 uncertain\nNT types (Usnea only)", loc="left", fontsize=7.5); ax.set_ylabel("cell types")
    save(fig, "fig4_usnea")


# ---------- Fig 5: male CNS ----------
def fig5():
    fig, axs = plt.subplots(1, 2, figsize=(5.6, 2.4))
    cal = json.load(open(os.path.join(M, "calib_male.json"))); lg = sorted(cal["log"])
    ax = axs[0]; ax.semilogy([x[0] for x in lg], [x[1] for x in lg], "o-", color=C_ACC, ms=3); ax.axhline(384, ls=":", color="#444", lw=.6)
    ax.set_xlabel("synaptic weight scale k"); ax.set_ylabel("active cells (sugar, 500 ms)"); ax.set_title("A  Male CNS calibration", loc="left", fontsize=7.5)
    rr = jl(os.path.join(M, "replicate_calib.jsonl")); b = {r["seed"]: sum(r["mn9"]) for r in rr if r["cond"] == "M_sugar"}; k = {r["seed"]: sum(r["mn9"]) for r in rr if r["cond"] == "M_sugar_koDNge031"}
    s = sorted(set(b) & set(k)); ax = axs[1]
    for i in s: ax.plot([0, 1], [b[i], k[i]], color=C_GRAY, lw=.6)
    ax.plot([0, 1], [np.mean([b[i] for i in s]), np.mean([k[i] for i in s])], color=C_BAD, lw=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["intact", "DNge031 KO"]); ax.set_ylabel("MN9 (spikes/s, L+R)")
    ax.set_title(f"B  DNge031 brake in male ({len(s)} seeds)", loc="left", fontsize=7.5)
    save(fig, "fig5_male_cns")




# ================= v1 (หลังตรวจค้านรอบ 1) =================
def ci95(d):
    d = np.asarray(d, float); return stats.t.ppf(.975, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))


def fig2v1():  # ตัด seed 0 (seed ที่ใช้คัดเลือก) + FDR เทียบ 400/402
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6))
    for ax, cond, title in zip(axs, ("sugar", "sugar_bitter"), ("Sugar", "Sugar + bitter")):
        rows = [r for f in glob.glob(os.path.join(K, "p2", f"ks_{cond}_p2_*.jsonl")) for r in jl(f) if r["seed"] != 0]
        base = {r["seed"]: sum(r["mn9"]) for r in rows if r["ko"] < 0}
        ko = {}
        for r in rows:
            if r["ko"] >= 0:
                ko.setdefault(r["ko"], []).append(100 * (sum(r["mn9"]) - base[r["seed"]]) / base[r["seed"]])
        t = pd.read_csv(os.path.join(K, f"p2_rank_{cond}_noseed0.csv")).head(16).iloc[::-1]
        p1 = pd.read_csv(os.path.join(K, f"p1_rank_{cond}.csv"))[["ko", "lit_name"]].drop_duplicates("ko")
        t = t.merge(p1, on="ko", how="left")
        m = [np.mean(ko[k]) for k in t.ko]; e = [ci95(ko[k]) for k in t.ko]
        lab = [f"{c} ({str(sd)[0].upper()})" + (f" [{n}]" if isinstance(n, str) else "") for c, sd, n in zip(t.cell_type, t.side, t.lit_name)]
        ax.barh(range(len(t)), m, xerr=e, color=[C_OK if x else C_GRAY for x in t.sig400], ecolor="#444", capsize=1.5)
        ax.set_yticks(range(len(t))); ax.set_yticklabels(lab, fontsize=6.5); ax.axvline(0, color="#444", lw=.6); ax.set_title(title, loc="left")
    fig.supxlabel("ΔMN9 after single-cell knockout (%; mean ± 95% CI, seeds 1–9)\ntop 16 of 64/63 re-tested cells shown; green = BH-FDR q<0.05 over all 400/402 screened cells and |Δ|≥5%", fontsize=6.8)
    save(fig, "fig2_knockout_screen")


def fig4v1():
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.6))
    ax = axs[0]; t = pd.read_csv(os.path.join(B, "results", "ntflip_scores.csv"), index_col=0)
    w = t.P2_water150.sort_values(); col = [C_OK if i == "CB0008" else (C_ACC if i == "none" else C_GRAY) for i in w.index]
    ax.bar(range(len(w)), w.values, color=col); ax.axhline(20, ls=":", color="#444", lw=.6)
    ax.set_xticks([]); ax.set_ylabel("water → MN9 (Hz, 4 seeds)")
    nb = w.drop(["CB0008", "none"]).max()
    ax.set_title(f"A  Re-sign each of 44 uncertain NT types:\nonly Usnea rescues water → MN9 (72.5 vs next {nb:.1f} Hz)", loc="left", fontsize=7.2)
    ax.text(len(w) - 1, w["CB0008"], "Usnea", ha="right", va="bottom", fontsize=7, color=C_OK)
    ax = axs[1]; labs = ["Rattle", "Usnea", "Roundup"]; x = np.arange(3); vals = {}
    for tag, base in (("D0_opm_s20-40", "O_sugar"), ("D0U_opm_s20-40", "O_sugar")):
        vals[tag] = [paired(os.path.join(B, f"results_raw_{tag}.jsonl"), base, f"O_sugar_ko{k}") for k in ("rattle", "usnea", "roundup")]
    ax.bar(x - .2, [v[0] for v in vals["D0_opm_s20-40"]], .4, yerr=[v[1] for v in vals["D0_opm_s20-40"]], color=C_GRAY, capsize=2, label="Usnea GABA (50 Hz)")
    ax.bar(x + .2, [v[0] for v in vals["D0U_opm_s20-40"]], .4, yerr=[v[1] for v in vals["D0U_opm_s20-40"]], color=C_OK, capsize=2, label="Usnea excitatory (45 Hz)")
    ax.axhline(0, color="#444", lw=.6); ax.set_xticks(x); ax.set_xticklabels(labs); ax.set_ylabel("ΔMN9 (%, mean ± SEM)")
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    ax.set_title("B  Matched operating point, fresh seeds 20–39\n(fly: ↓, ↓, no change)", loc="left", fontsize=7.2)
    save(fig, "fig4_usnea")


def fig5v1():
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.5))
    ax = axs[0]; path = os.path.join(B, "results_raw_D0_sides_s20-30.jsonl")
    kk = lambda r: sum(r["mn9b"]); x = np.arange(3)
    for off, st, c, lab in ((-.2, "sugar", C_BAD, "sugar GRNs, side A (20)"), (.2, "sugar_other", "#F0997B", "sugar GRNs, side B (9)")):
        v = [paired(path, f"S_{st}_none", f"S_{st}_{k}", key=kk) for k in ("L", "R", "LR")]
        ax.bar(x + off, [a[0] for a in v], .4, yerr=[a[1] for a in v], color=c, capsize=2, label=lab)
    ax.set_xticks(x); ax.set_xticklabels(["KO left", "KO right", "KO both"]); ax.set_ylabel("ΔMN9 L+R (%, mean ± SEM)")
    ax.legend(fontsize=5.8, frameon=False); ax.set_title("A  Female: DNge031 knockout\n(10 fresh seeds)", loc="left", fontsize=7.2)
    ax = axs[1]; ks, eff, err, base = [], [], [], []
    for k, f, sd in ((0.50, "replicate_k0.500.jsonl", range(20, 30)), (0.52, "replicate_k0.520.jsonl", range(20, 30)), (0.547, "replicate_calib.jsonl", range(20, 30))):
        rr = [r for r in jl(os.path.join(M, f)) if r["seed"] in sd]
        b = {r["seed"]: sum(r["mn9"]) for r in rr if r["cond"] == "M_sugar"}; q = {r["seed"]: sum(r["mn9"]) for r in rr if r["cond"] == "M_sugar_koDNge031"}
        s = sorted(set(b) & set(q)); d = np.array([q[i] - b[i] for i in s], float); bm = np.mean([b[i] for i in s])
        ks.append(k); eff.append(100 * d.mean() / bm); err.append(100 * d.std(ddof=1) / np.sqrt(len(d)) / bm); base.append(bm)
    ax.errorbar(ks, eff, yerr=err, fmt="o-", color=C_BAD, capsize=2); ax.axhline(0, color="#444", lw=.6)
    for k, e, bm in zip(ks, eff, base): ax.annotate(f"MN9 {bm:.0f}", (k, e), textcoords="offset points", xytext=(4, -10), fontsize=6)
    ax.set_xlabel("synaptic scale k (male CNS)"); ax.set_ylabel("ΔMN9 after DNge031 KO (%, ± SEM)"); ax.set_title("B  Male: brake depends on\nnetwork drive (seeds 20–29)", loc="left", fontsize=7.2)
    ax = axs[2]; rr = [r for r in jl(os.path.join(M, "replicate_calib.jsonl")) if r["seed"] >= 20]
    b = {r["seed"]: sum(r["mn9"]) for r in rr if r["cond"] == "M_sugar"}; q = {r["seed"]: sum(r["mn9"]) for r in rr if r["cond"] == "M_sugar_koRoundup"}
    s = sorted(set(b) & set(q))
    for i in s: ax.plot([0, 1], [b[i], q[i]], color=C_GRAY, lw=.6)
    ax.plot([0, 1], [np.mean([b[i] for i in s]), np.mean([q[i] for i in s])], color=C_BAD, lw=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["intact", "Roundup KO"]); ax.set_ylabel("MN9 (spikes/s, L+R)")
    ax.set_title("C  Male (with VNC): Roundup\nremains a bottleneck", loc="left", fontsize=7.2)
    save(fig, "fig5_dnge031_male")


if __name__ == "__main__":
    for f in (fig1, fig2v1, fig3, fig4v1, fig5v1):
        f()

"""สำรวจ (post-hoc): ถ้า Usnea (CB0008) จริงๆ เป็นเซลล์กระตุ้น (NT ทาย GABA มั่นใจแค่ 49–66%) ผลปิด Usnea จะตรงแมลงจริงไหม
แมลงจริง (Shiu 2022): hyperpolarize Usnea → PER ลด · โมเดลเดิม (GABA): ปิดแล้ว MN9 เพิ่ม +10–24%"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from flybrain import load_shiu_783, Brain

W, ids, ann = load_shiu_783()
a = ann.drop_duplicates("root_id").set_index("root_id").reindex(ids)
usnea = np.flatnonzero(a.cell_type.values == "CB0008")
L = json.load(open(os.path.join(HERE, "ids_shiu_figures.json")))
idset = set(ids.tolist())
W2 = W.tocsc(copy=True)  # คอลัมน์ = เซลล์ต้นทาง → กลับเครื่องหมายขาออกของ Usnea
for j in usnea:
    W2.data[W2.indptr[j]:W2.indptr[j + 1]] = np.abs(W2.data[W2.indptr[j]:W2.indptr[j + 1]])
W2 = W2.tocsr()
out = open(os.path.join(HERE, "probe_usnea_sign.jsonl"), "w")
for tag, Wm in (("orig", W), ("usnea_exc", W2)):
    b = Brain(Wm, ids); sug = b.ix([x for x in L["neu_sugar"] if x in idset]); mn9 = b.ix([L["id_mn9"][0]])[0]
    for hz in (50, 150):
        for ko in (False, True):
            for seed in range(3):
                b.rng = np.random.default_rng(seed)
                c, _ = b.run(1000.0, rate_fn=lambda s: (sug, float(hz)), silence=usnea if ko else None)
                out.write(json.dumps(dict(model=tag, hz=hz, ko=ko, seed=seed, mn9=int(c[mn9]))) + chr(10)); out.flush()
print("done")

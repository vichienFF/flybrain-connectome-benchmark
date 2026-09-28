"""
แปลง Male CNS v1.0 (Janelia, CC-BY) → เมทริกซ์สำหรับตัวจำลอง LIF เดียวกับ flybrain.py (Shiu 2024)
- เซลล์: ตัดออก Glia / Orphan (ชิ้นส่วน) / Unimportant · เก็บทั้งสมอง + optic lobe + เส้นประสาทลำตัว (VNC)
- เครื่องหมาย: consensus_nt (สำรอง celltype_predicted_nt) · ACh/DA/5HT/OA = +, GABA/Glu = − (สมมติฐานเดียวกับโมเดลตัวเมีย) · ไม่ทราบ = +
- น้ำหนัก = จำนวน synapse × เครื่องหมาย (Shiu ใช้ทุกเส้น ไม่ตัด threshold)
- อ่านไฟล์ weights 1.05 GB ทีละ record batch (ประหยัด RAM)
→ malecns/model/W_male.npz (CSR post×pre), ids.npy, ann.parquet
"""
import os, time
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.ipc as ipc, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "model")
SIGN = {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "dopamine": 1, "serotonin": 1, "octopamine": 1}


def main():
    t0 = time.time(); os.makedirs(OUT, exist_ok=True)
    a = pd.read_feather(os.path.join(HERE, "body-annotations-male-cns-v1.0-minconf-0.5.feather"))
    a = a[~a.statusLabel.isin(["Glia", "Orphan", "Unimportant"]) & a.superclass.notna()].copy()
    nt = pd.read_feather(os.path.join(HERE, "body-neurotransmitters-male-cns-v1.0.feather"))
    nt = nt.drop_duplicates("body").set_index("body")
    ntc = nt.consensus_nt.where(nt.consensus_nt.isin(SIGN.keys()), nt.celltype_predicted_nt)
    a["nt"] = a.bodyId.map(ntc)
    a["sign"] = a.nt.map(SIGN).fillna(1).astype(np.int8)
    a = a.sort_values("bodyId").reset_index(drop=True)
    ids = a.bodyId.values.astype(np.int64)
    idx = pd.Index(ids)
    sign = a.sign.values
    print(f"neurons {len(ids):,} · signs + {int((sign > 0).sum()):,} / − {int((sign < 0).sum()):,} · {time.time() - t0:.0f}s", flush=True)
    rows, cols, vals = [], [], []
    with pa.memory_map(os.path.join(HERE, "connectome-weights-male-cns-v1.0-minconf-0.5.feather")) as src:
        rd = ipc.open_file(src)
        for bi in range(rd.num_record_batches):
            b = rd.get_batch(bi)
            pre = idx.get_indexer(b.column("body_pre").to_numpy()); post = idx.get_indexer(b.column("body_post").to_numpy())
            ok = (pre >= 0) & (post >= 0) & (pre != post)
            w = b.column("weight").to_numpy()[ok].astype(np.float32)
            rows.append(post[ok].astype(np.int32)); cols.append(pre[ok].astype(np.int32)); vals.append(w * sign[pre[ok]])
            if bi % 20 == 0:
                print(f"batch {bi + 1}/{rd.num_record_batches} kept {sum(len(r) for r in rows):,} · {time.time() - t0:.0f}s", flush=True)
    r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
    W = sp.csr_matrix((v, (r, c)), shape=(len(ids), len(ids)))
    W.sum_duplicates()
    sp.save_npz(os.path.join(OUT, "W_male.npz"), W)
    np.save(os.path.join(OUT, "ids.npy"), ids)
    a.drop(columns=[c for c in a.columns if a[c].dtype == object and c not in
                    ("type", "flywireType", "somaSide", "superclass", "class", "subclass", "synonyms", "nt", "statusLabel")]).to_parquet(os.path.join(OUT, "ann.parquet"))
    print(f"W {W.shape} nnz {W.nnz:,} synapses {int(np.abs(W.data).sum()):,} · done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

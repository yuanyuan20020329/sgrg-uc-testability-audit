# -*- coding: utf-8 -*-
"""数据集1：X 归一化状态 + SGRG 基因覆盖 + 细胞组成"""
import collections
import json
import numpy as np
import h5py

PATH = "data/raw/sc/CxG_IBD_HC_10x.h5ad"
SGRG = [l.split("\t")[0] for l in
        open("results/sgrg_gene_set_clean.tsv", encoding="utf-8").read().splitlines()[1:] if l.strip()]
print("SGRG 面板基因数：", len(SGRG))


def decode(arr):
    return [x.decode() if isinstance(x, bytes) else x for x in arr]


f = h5py.File(PATH, "r")
X = f["X"]
n_cells, n_genes = [int(v) for v in X.attrs["shape"]]
d = X["data"][:]
print("X.data: min=%.4f max=%.4f mean=%.4f  整数占比=%.4f  nnz/cell=%.0f" % (
    d.min(), d.max(), d.mean(), float(np.mean(d == np.floor(d))), d.size / n_cells))

sym = np.array(decode(f["var"]["gene_symbols"][:]))
print("var 基因数：", len(sym), " 唯一 symbol：", len(set(sym)))
present = [g for g in SGRG if g in set(sym)]
missing = [g for g in SGRG if g not in set(sym)]
print("SGRG 可检出 symbol：%d / %d" % (len(present), len(SGRG)))
print("缺失：", missing)

lut = np.full(n_genes, -1, dtype=np.int32)
gene_idx = np.array([np.where(sym == g)[0][0] for g in present], dtype=np.int64)
lut[gene_idx] = np.arange(len(gene_idx), dtype=np.int32)
G = len(gene_idx)

counts = np.zeros(G, dtype=np.int64)
indptr, indices, data = X["indptr"], X["indices"], X["data"]
BS = 4000
for s in range(0, n_cells, BS):
    e = min(s + BS, n_cells)
    lo, hi = int(indptr[s]), int(indptr[e])
    cols = indices[lo:hi]
    vals = data[lo:hi]
    pos = lut[cols]
    sel = (pos >= 0) & (vals > 0)
    if sel.any():
        counts += np.bincount(pos[sel], minlength=G)

det = counts / n_cells
order = np.argsort(-det)
rows = [(present[r], int(counts[r]), float(det[r])) for r in order]
print("\n top15:", ", ".join("%s=%.1f%%" % (g, 100 * v) for g, _, v in rows[:15]))
print(" bot15:", ", ".join("%s=%.1f%%" % (g, 100 * v) for g, _, v in rows[-15:]))
print(" 检测率>=1%% 的 SGRG：%d    >=5%%：%d" % (int((det >= .01).sum()), int((det >= .05).sum())))

obs = f["obs"]
ct_cat, ct_code = decode(obs["cell_type"]["categories"][:]), obs["cell_type"]["codes"][:]
dis_cat, dis_code = decode(obs["disease"]["categories"][:]), obs["disease"]["codes"][:]
don_cat, don_code = decode(obs["donor_id"]["categories"][:]), obs["donor_id"]["codes"][:]

print("\n疾病 × 细胞类型 细胞数：")
tab = collections.Counter((dis_cat[a], ct_cat[b]) for a, b in zip(dis_code, ct_code))
allct = sorted(set(ct_cat))
print("  %-16s" % "disease", "".join("%13s" % c[:12] for c in allct))
for d0 in sorted(set(dis_cat)):
    print("  %-16s" % d0, "".join("%13d" % tab[(d0, c)] for c in allct))

print("\n供体 × 疾病：")
dc = collections.Counter((don_cat[a], dis_cat[b]) for a, b in zip(don_code, dis_code))
for k in sorted(dc):
    print("   ", k, dc[k])

f.close()
json.dump({"n_cells": n_cells, "n_genes": n_genes, "sgrg_panel": len(SGRG),
           "sgrg_present": present, "sgrg_missing": missing,
           "sgrg_detection": {g: v for g, _, v in rows},
           "x_stats": {"min": float(d.min()), "max": float(d.max()), "mean": float(d.mean()),
                       "integer_fraction": float(np.mean(d == np.floor(d)))}},
          open("results/partD_ds1_overview.json", "w", encoding="utf-8"), indent=2, ensure_ascii=False)
print("\n[saved] results/partD_ds1_overview.json")

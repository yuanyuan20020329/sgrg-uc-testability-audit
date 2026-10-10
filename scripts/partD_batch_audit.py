# -*- coding: utf-8 -*-
"""Part D / 批次-疾病嵌套审计

为什么做
--------
ds4（Tearle2026）的「UC vs 正常」比较在控制测序深度后全灭。但核查 obs 元数据发现
一个**更上游**的问题：`batch_condition` 与疾病标签几乎完全嵌套
（35 个批次里只有 1 个同时含 UC 与正常；UC 只在 2022_3~6 + 2025，正常只在 A3x/CRC 切除批次），
而且批次间的细胞中位深度相差 10~45 倍。

如果疾病与处理批次不可分离，那么「疾病效应」与「技术效应」在数学上不可识别——
深度校正只能处理深度这一条通道，处理不掉批次带来的其他差异（化学版本、样本类型、
组织处理、解离效率、ambient RNA 水平……）。

本脚本对**四套数据**统一做同一件事：
  1) 找出可用的「批次类」变量（study / batch_condition / assay / donor）；
  2) 检验 disease × 批次 是否嵌套（是否存在同时含病例与对照的批次）；
  3) 算每个 donor / 批次的**细胞级中位深度**，看深度是随疾病分层还是随批次分层；
  4) 给出「在批次内可识别的病例-对照对比」还剩多少 donor。

输出
----
  results/partD_batch_audit_donor.tsv     每个 donor 一行（疾病/批次/深度/细胞数）
  results/partD_batch_audit_summary.tsv   每套数据一行（嵌套指标）
  results/partD_batch_audit_summary.json
用法
----
  python partD_batch_audit.py
"""
import json
import os
import time

import h5py
import numpy as np
import pandas as pd

T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)
OUT = "results"

# (tag, 路径, 疾病匹配 UC, 疾病匹配 对照, 批次类变量候选)
DATASETS = [
    ("ds1_garrido", "data/raw/sc/CxG_IBD_HC_10x.h5ad",
     "ulcer", "normal", ["assay"]),
    ("ds2_myeloid", "data/raw/sc/CxG_Oliver2024_Myeloid.h5ad",
     "ulcer", "normal", ["study", "assay", "sample_retrieval"]),
    ("ds3_largeintestine", "data/raw/sc/CxG_Oliver2024_LargeIntestine.h5ad",
     "ulcer", "normal", ["study", "assay", "sample_retrieval"]),
    ("ds4_tearle", "data/raw/sc/CxG_Tearle2026_PSC_UC.h5ad",
     "ulcer", "normal", ["batch_condition", "assay"]),
]


# ------------------------------------------------------------------ 读 obs 列
def obs_col(f, col):
    if col not in f["obs"]:
        return None
    d = f["obs"][col]
    if isinstance(d, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x)
                         for x in d["categories"][:]], dtype=object)
        codes = d["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    arr = d[:]
    if arr.dtype.kind in "OS":
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr], dtype=object)
    return arr


def pick_matrix(f):
    """返回 (名称, h5py 组)：优先整数 counts"""
    cand = []
    if "X" in f:
        cand.append(("X", f["X"]))
    if "raw" in f and "X" in f["raw"]:
        cand.append(("raw/X", f["raw"]["X"]))
    if "layers" in f:
        for k in f["layers"].keys():
            cand.append(("layers/" + k, f["layers"][k]))
    for nm, g in cand:
        try:
            d = g["data"][:200000]
            if float(np.mean(d == np.floor(d))) > 0.999:
                return nm, g
        except Exception:
            pass
    return cand[0] if cand else (None, None)


def cell_totals(grp, block=100000):
    """流式求每细胞总 counts（只读 data 段）"""
    ip = grp["indptr"][:]
    n = len(ip) - 1
    tot = np.zeros(n, dtype=np.int64)
    for i0 in range(0, n, block):
        b = np.arange(i0, min(i0 + block, n))
        s = int(ip[b[0]]); e = int(ip[b[-1] + 1])
        if e <= s:
            continue
        off = np.append(ip[b] - s, e - s)
        seg = np.diff(off)
        d = grp["data"][s:e].astype(np.int64)
        cs = np.add.reduceat(d, np.minimum(off[:-1], max(len(d) - 1, 0)))
        cs[seg == 0] = 0
        tot[i0:i0 + len(b)] = cs
    return tot


donor_rows, summ_rows = [], []

for tag, path, uc_pat, ctl_pat, batchvars in DATASETS:
    if not os.path.exists(path):
        log("跳过（缺失）", tag); continue
    log("=" * 70)
    log("审计", tag, "|", path)
    f = h5py.File(path, "r")
    disease = obs_col(f, "disease")
    donor = obs_col(f, "donor_id")
    if donor is None:
        donor = np.array(["d%d" % i for i in range(len(disease))], dtype=object)
    if disease is None:
        disease = np.array(["unknown"] * len(donor), dtype=object)

    batchvars = [b for b in batchvars if b in f["obs"]]
    log("  可用批次类变量:", batchvars)
    bv = {b: obs_col(f, b) for b in batchvars}

    mname, grp = pick_matrix(f)
    log("  矩阵 =", mname, grp.attrs.get("shape"))
    tot = cell_totals(grp)

    df = pd.DataFrame({"donor": donor, "disease": disease, "nUMI": tot})
    for b in batchvars:
        df[b] = bv[b]
    df["_grp"] = np.where(pd.Series(disease).str.lower().str.contains(uc_pat), "UC",
                  np.where(pd.Series(disease).str.lower().str.contains(ctl_pat), "CTRL", "OTHER"))

    # donor 级汇总（批次取该 donor 的主批次）
    g = df.groupby(["donor", "_grp"]).agg(n_cells=("nUMI", "size"),
                                          med_nUMI=("nUMI", "median"),
                                          mean_nUMI=("nUMI", "mean")).reset_index()
    for b in batchvars:
        dom = df.groupby(["donor", b]).size().reset_index(name="n")
        dom = dom.sort_values("n").groupby("donor").tail(1).set_index("donor")[b]
        g[b] = g.donor.map(dom)
    g["dataset"] = tag
    donor_rows.append(g)

    uc = g[g._grp == "UC"]; ct = g[g._grp == "CTRL"]
    rec = {"dataset": tag, "path": path, "matrix": mname,
           "n_cells": int(len(df)), "n_donors_total": int(g.donor.nunique()),
           "n_uc_donors": int(len(uc)), "n_ctrl_donors": int(len(ct)),
           "med_nUMI_uc": float(np.median(df[df._grp == "UC"].nUMI)) if len(uc) else np.nan,
           "med_nUMI_ctrl": float(np.median(df[df._grp == "CTRL"].nUMI)) if len(ct) else np.nan,
           "depth_fold_uc_over_ctrl": (float(np.median(df[df._grp == "UC"].nUMI) /
                                             max(np.median(df[df._grp == "CTRL"].nUMI), 1e-9))
                                       if len(uc) and len(ct) else np.nan)}
    # is_primary_data（若存在）：区分「本研究原产」与「复用他人数据」
    prim = obs_col(f, "is_primary_data")
    if prim is not None:
        pf = pd.Series(np.asarray(prim).astype(bool), index=df.index)
        ctrl_mask = df._grp == "CTRL"
        rec["n_ctrl_cells_primary"] = int((ctrl_mask & pf).sum())
        rec["n_ctrl_cells_total"] = int(ctrl_mask.sum())
        rec["frac_ctrl_cells_primary"] = (float((ctrl_mask & pf).sum() / max(ctrl_mask.sum(), 1)))
        rec["n_ctrl_donors_primary"] = int(g[(g._grp == "CTRL") & g.donor.isin(
            df.loc[pf, "donor"].unique())].donor.nunique())
    # 嵌套指标（对每个批次类变量）
    for b in batchvars:
        dsg = df[df._grp != "OTHER"]
        both = dsg.groupby(b)._grp.nunique()
        rec["n_batch_" + b] = int(dsg[b].nunique())
        rec["n_batch_mixed_" + b] = int((both > 1).sum())
        # donor 级：UC 中有多少 donor 的批次里同时存在对照 donor
        db = g[g._grp != "OTHER"].groupby(b)._grp.nunique()
        mixed = set(db[db > 1].index)
        rec["n_uc_donor_in_mixed_batch_" + b] = int(uc[b].isin(mixed).sum())
        rec["n_ctrl_donor_in_mixed_batch_" + b] = int(ct[b].isin(mixed).sum())
        # ★ 同源对照比例：对照细胞/供体中，其批次同时含病例的比例
        cmask = df._grp == "CTRL"
        rec["frac_ctrl_cells_samebatch_" + b] = (
            float(df.loc[cmask, b].isin(mixed).mean()) if int(cmask.sum()) else np.nan)
        rec["frac_ctrl_donors_samebatch_" + b] = (
            float(ct[b].isin(mixed).mean()) if len(ct) else np.nan)
    summ_rows.append(rec)

    log("  UC donor %d | 对照 donor %d | 细胞级中位深度 UC %.0f vs 对照 %.0f（%.1f×）" % (
        rec["n_uc_donors"], rec["n_ctrl_donors"], rec["med_nUMI_uc"], rec["med_nUMI_ctrl"],
        rec["depth_fold_uc_over_ctrl"]))
    for b in batchvars:
        log("    [%s] 批次数 %d | 同时含病例与对照的批次 %d | 落在混合批次里的 UC donor %d / 对照 %d" % (
            b, rec["n_batch_" + b], rec["n_batch_mixed_" + b],
            rec["n_uc_donor_in_mixed_batch_" + b], rec["n_ctrl_donor_in_mixed_batch_" + b]))
    f.close()

D = pd.concat(donor_rows, ignore_index=True)
S = pd.DataFrame(summ_rows)
D.to_csv(os.path.join(OUT, "partD_batch_audit_donor.tsv"), sep="\t", index=False,
         float_format="%.6g")
S.to_csv(os.path.join(OUT, "partD_batch_audit_summary.tsv"), sep="\t", index=False,
         float_format="%.6g")
json.dump(json.loads(S.to_json(orient="records")),
          open(os.path.join(OUT, "partD_batch_audit_summary.json"), "w", encoding="utf-8"),
          indent=2, ensure_ascii=False, default=str)

print("\n" + "=" * 96)
print("批次-疾病嵌套审计汇总")
print("=" * 96)
print(S.to_string(index=False, float_format=lambda v: "%.0f" % v))
log("完成 %.1f s" % (time.time() - T0))

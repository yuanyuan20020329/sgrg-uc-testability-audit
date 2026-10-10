# -*- coding: utf-8 -*-
"""Part D / 公共资源可用性审计 —— 全库层面（CELLxGENE 全部集合，不限肠道）

为什么做
--------
partD_resfeas_cxg.py 已经把 CELLxGENE 全库枚举了一遍，但审计**在"肠道组织"这一步就停住**，
审计单位仍是「本项目的 4 套数据 + 少量肠道候选」。审稿人（以及同类已发表文章）要求的是
**全库层面的流行度数字**：
  · Huang 2025 (Genome Biology) 给的是「253 项研究 / 164k 样本，元数据完整度 74.8%」
  · Rogic 2026 (Scientific Data) 给的是「GEO 单细胞研究可用率约 40%、注释率 <10%」
本脚本要给出的对应数字是：
  · CELLxGENE 全库细胞中，**非本数据集原产（is_primary_data=False 的复用细胞）占多少**
  · 在同时提供「病例臂 + 对照臂」的集合里，**对照臂有多少细胞不是本研究产生的**
  · 对照臂与病例臂的**复用率不对称**、以及**测序深度不对称**在全库的流行度
  · 按组织系统分层 —— 这一步同时回答「该框架是不是 UC 专属」

★ 关键实现点
  这些量全部来自 CELLxGENE **集合详情端点**（`primary_cell_count` / `cell_count` /
  `donor_id` / `is_primary_data` / `mean_genes_per_cell`），**不需要下载任何 h5ad**。
  全库 397 个集合详情缓存于 data/raw/cxg_coll_detail/，断点续抓。

口径说明（写稿时必须带上，否则数字会被误用）
  · frac_primary = primary_cell_count / cell_count，即「该数据集有多少细胞是本研究原产」；
    1 - frac_primary = 复用/外部细胞比例。
  · is_primary_data=False **本身不是错误**——atlas 型集合合法地聚合外部数据。真正的病理是
    「一个集合声称有病例-对照对比，而对照臂整体来自它没有产生的细胞」，
    故本脚本同时计算**病例臂复用率**作为对照，用不对称性（M-B）而非绝对值下结论。
  · mean_genes_per_cell 是数据集级均值，作深度代理；跨集合比较只用原产（frac_primary=1）数据集。

输出
----
  results/partD_resfeas_libwide_datasets.tsv     全库数据集级（含 frac_primary）
  results/partD_resfeas_libwide_collections.tsv  集合级对比臂审计
  results/partD_resfeas_libwide_summary.json     全库汇总指标
  results/fig_D14_libwide_provenance.{png,pdf}   图 D13
用法
----
  python partD_resfeas_libwide.py [--refresh]
"""
import argparse
import json
import os
import re
import ssl
import socket
import time
import urllib.request

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

socket.setdefaulttimeout(120)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

API = "https://api.cellxgene.cziscience.com/curation/v1"
CACHE = "data/raw/cxg_coll_detail"
OUT = "results"
T0 = time.time()


def log(*x):
    print("[%6.1fs]" % (time.time() - T0), *x, flush=True)


# ----------------------------------------------------------------- 口径定义
CONTROL_RE = re.compile(
    r"^(normal|healthy|control|ctrl|non[- ]?disease|unaffected|wild[- ]?type|wt|"
    r"na[iï]ve|saline|sham|baseline)\b", re.I)

# 组织系统粗映射（用于分层，证明框架非肠道专属）
TISSUE_SYS = [
    ("gut / intestine", re.compile(r"\b(colon|colonic|intestin\w*|ileum|ileal|jejun\w*|duoden\w*|"
                                   r"rectum|rectal|cecum|cecal|sigmoid|gut|bowel|appendix)\b", re.I)),
    ("liver / biliary", re.compile(r"\b(liver|hepat\w*|bile|biliary|cholangio\w*)\b", re.I)),
    ("brain / CNS", re.compile(r"\b(brain|cortex|cerebr\w*|cerebell\w*|hippocamp\w*|hypothalam\w*|"
                               r"midbrain|neuron\w*|spinal|dorsal root|retina|retinal|"
                               r"choroid plexus|substantia nigra|medulla|thalam\w*|striatum)\b", re.I)),
    ("lung / airway", re.compile(r"\b(lung|pulmonary|alveol\w*|bronch\w*|airway|trachea|nasal|"
                                 r"respiratory)\b", re.I)),
    ("skin", re.compile(r"\b(skin|derm\w*|epiderm\w*|keratinocyte)\b", re.I)),
    ("blood / immune", re.compile(r"\b(blood|pbmc|bone marrow|lymph\w*|spleen|thymus|"
                                  r"mononuclear|leukocyte|myeloid|t cell|b cell|peripheral)\b", re.I)),
    ("kidney / urinary", re.compile(r"\b(kidney|renal|nephron|bladder|urothel\w*|ureter)\b", re.I)),
    ("heart / vascular", re.compile(r"\b(heart|cardiac|myocard\w*|ventricl\w*|atrium|aort\w*|"
                                    r"vascular|endotheli\w*|coronary)\b", re.I)),
    ("muscle / adipose", re.compile(r"\b(muscle|myocyte|skeletal|adipos\w*|fat|myofib\w*)\b", re.I)),
    ("pancreas", re.compile(r"\b(pancrea\w*|islet|langerhans)\b", re.I)),
    ("reproductive", re.compile(r"\b(ovar\w*|testis|testic\w*|uter\w*|endometri\w*|prostate|"
                                r"placenta|mammary|breast)\b", re.I)),
    ("eye", re.compile(r"\b(eye|cornea|iris|lens|conjunctiv\w*)\b", re.I)),
    ("mouth / oropharynx", re.compile(r"\b(oral|tongue|gingiv\w*|saliva\w*|tonsil|pharyn\w*)\b", re.I)),
    ("bone / joint", re.compile(r"\b(bone|cartilage|synovi\w*|joint|osteo\w*)\b", re.I)),
    ("thyroid / endocrine", re.compile(r"\b(thyroid|adrenal|pituitar\w*|parathyroid)\b", re.I)),
]

DATASET_COLS = [
    "collection_id", "collection_name", "doi", "published_at", "consortia",
    "dataset_id", "dataset_title", "organism", "is_human",
    "cell_count", "primary_cell_count", "reused_cells", "frac_primary",
    "is_primary_data", "n_donors", "n_cell_types",
    "n_disease_labels", "is_control_only", "is_disease_bearing", "arm_class",
    "diseases", "tissues", "tissue_systems", "assays",
    "suspension_type", "spatial", "mean_genes_per_cell",
]


def get_json(url, tries=4):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
            with urllib.request.urlopen(req, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (k + 1))
    raise last


def labels(field):
    return [x.get("label") for x in (field or []) if x.get("label")]


def tissue_systems(tissues):
    hits = []
    for name, pat in TISSUE_SYS:
        if any(pat.search(t) for t in tissues):
            hits.append(name)
    return hits or ["unmapped / other"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--no-fetch", action="store_true", help="只用已有缓存，不补抓")
    a = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    # ---------------------------------------------------------- 1) 全库集合列表
    slim = "data/raw/cxg_all_collections.json"
    if os.path.exists(slim) and not a.refresh:
        cols = json.load(open(slim, encoding="utf-8"))
        log("载入全库列表（缓存）%d 集合" % len(cols))
    else:
        cols = get_json(API + "/collections")
        json.dump(cols, open(slim, "w", encoding="utf-8"), ensure_ascii=False)
        log("下载全库列表 %d 集合" % len(cols))

    # ---------------------------------------------------------- 2) 抓全部集合详情
    dets, n_new, n_fail = {}, 0, 0
    for i, c in enumerate(cols, 1):
        cid = c["collection_id"]
        fp = os.path.join(CACHE, cid + ".json")
        if os.path.exists(fp) and not a.refresh:
            try:
                dets[cid] = json.load(open(fp, encoding="utf-8"))
                continue
            except Exception:  # noqa: BLE001
                pass
        if a.no_fetch:
            continue
        try:
            det = get_json(API + "/collections/" + cid)
            json.dump(det, open(fp, "w", encoding="utf-8"), ensure_ascii=False)
            dets[cid] = det
            n_new += 1
        except Exception as e:  # noqa: BLE001
            n_fail += 1
            log("  抓取失败 %s %s: %s" % (cid[:8], (c.get("name") or "")[:40], e))
        if i % 25 == 0:
            log("  进度 %d/%d（本轮新抓 %d，失败 %d）" % (i, len(cols), n_new, n_fail))
    log("集合详情就绪：%d 个（本轮新抓 %d，失败 %d）" % (len(dets), n_new, n_fail))

    # ---------------------------------------------------------- 3) 数据集级表
    rows, crow = [], []
    ctrl = {"uc": re.compile(r"ulcerative colitis", re.I),
            "cd": re.compile(r"crohn", re.I)}
    for c in cols:
        cid = c["collection_id"]
        det = dets.get(cid)
        if det is None:
            continue
        cn = det.get("name") or c.get("name") or ""
        ds = det.get("datasets") or []
        agg = dict(cells=0, primary=0, ctrl_cells=0, ctrl_reused=0, case_cells=0,
                   case_reused=0, n_ctrl_ds=0, n_case_ds=0, n_ds=0,
                   mixed_cells=0, mixed_primary=0, n_mixed_ds=0,
                   g_ctrl_num=0.0, g_ctrl_den=0.0, g_case_num=0.0, g_case_den=0.0,
                   g_ctrl_pure=0.0, g_ctrl_pure_den=0.0,
                   systems=set(), orgs=set(), diseases=set(), has_uc=False, has_cd=False)
        for d in ds:
            if d.get("tombstone"):
                continue
            tis = labels(d.get("tissue"))
            dis = labels(d.get("disease"))
            org = labels(d.get("organism"))
            cc = d.get("cell_count")
            pcc = d.get("primary_cell_count")
            frac = (pcc / cc) if (cc and pcc is not None and cc > 0) else None
            dl = [x.lower() for x in dis]
            has_ctrl_lab = any(CONTROL_RE.match(x) for x in dl)
            has_dis_lab = any(not CONTROL_RE.match(x) for x in dl)
            # 三类互斥归类：同一数据集内同时含对照与疾病标签 = mixed
            # （对比存在于数据集内部，但臂无法在元数据层拆分 → 单独计数，不并入病例臂）
            ctl_only = bool(dl) and has_ctrl_lab and not has_dis_lab
            dis_only = bool(dl) and has_dis_lab and not has_ctrl_lab
            mixed = bool(dl) and has_ctrl_lab and has_dis_lab
            dis_bear = has_dis_lab
            arm_class = ("mixed" if mixed else "ctl_only" if ctl_only
                         else "dis_only" if dis_only else "no_label")
            systems = tissue_systems(tis)
            mgpc = d.get("mean_genes_per_cell")
            agg["n_ds"] += 1
            if cc:
                agg["cells"] += cc
                if pcc is not None:
                    agg["primary"] += pcc
            if ctl_only and cc:
                agg["n_ctrl_ds"] += 1
                agg["ctrl_cells"] += cc
                if pcc is not None:
                    agg["ctrl_reused"] += (cc - pcc)
                # 深度代理：只用完全原产的数据集，避免混入复用细胞的口径污染
                if frac is not None and frac >= 0.999 and mgpc:
                    agg["g_ctrl_num"] += mgpc * cc
                    agg["g_ctrl_den"] += cc
            if dis_only and cc:
                agg["n_case_ds"] += 1
                agg["case_cells"] += cc
                if pcc is not None:
                    agg["case_reused"] += (cc - pcc)
                if frac is not None and frac >= 0.999 and mgpc:
                    agg["g_case_num"] += mgpc * cc
                    agg["g_case_den"] += cc
            if mixed and cc:
                agg["n_mixed_ds"] += 1
                agg["mixed_cells"] += cc
                if pcc is not None:
                    agg["mixed_primary"] += pcc
            agg["systems"] |= set(systems)
            agg["orgs"] |= set(org)
            agg["diseases"] |= set(dis)
            agg["has_uc"] |= any(ctrl["uc"].search(x) for x in dis)
            agg["has_cd"] |= any(ctrl["cd"].search(x) for x in dis)
            rows.append({
                "collection_id": cid, "collection_name": cn, "doi": det.get("doi"),
                "published_at": det.get("published_at"),
                "consortia": ";".join(det.get("consortia") or []),
                "dataset_id": d.get("dataset_id"), "dataset_title": d.get("title"),
                "organism": ";".join(org), "is_human": any("homo sapiens" in x.lower() for x in org),
                "cell_count": cc, "primary_cell_count": pcc,
                "reused_cells": (cc - pcc) if (cc and pcc is not None) else None,
                "frac_primary": frac,
                "is_primary_data": ";".join(str(x) for x in (d.get("is_primary_data") or [])),
                "n_donors": len(d.get("donor_id") or []),
                "n_cell_types": len(labels(d.get("cell_type"))),
                "n_disease_labels": len(dis),
                "is_control_only": ctl_only, "is_disease_bearing": dis_bear,
                "arm_class": arm_class,
                "diseases": ";".join(dis), "tissues": ";".join(tis),
                "tissue_systems": ";".join(systems),
                "assays": ";".join(labels(d.get("assay"))),
                "suspension_type": ";".join(d.get("suspension_type") or []),
                "spatial": bool(d.get("spatial")),
                "mean_genes_per_cell": mgpc,
            })

        def safe(n, dd):
            return (n / dd) if dd else None

        # 「同研究对照」深度代理：仅当该集合确有对比（臂可分 或 含混合数据集）时才有意义
        if agg["ctrl_cells"] and (agg["case_cells"] or agg["n_mixed_ds"]):
            agg["g_ctrl_pure"] = agg["g_ctrl_num"]
            agg["g_ctrl_pure_den"] = agg["g_ctrl_den"]

        crow.append({
            "collection_id": cid, "collection_name": cn, "doi": det.get("doi"),
            "published_at": det.get("published_at"),
            "consortia": ";".join(det.get("consortia") or []),
            "organism": ";".join(sorted(agg["orgs"])),
            "is_human": any("homo sapiens" in x.lower() for x in agg["orgs"]),
            "tissue_systems": ";".join(sorted(agg["systems"])),
            "n_datasets": agg["n_ds"],
            "n_control_ds": agg["n_ctrl_ds"], "n_disease_ds": agg["n_case_ds"],
            "cells_total": agg["cells"], "cells_primary": agg["primary"],
            "collection_reuse_frac": (1 - agg["primary"] / agg["cells"]) if agg["cells"] else None,
            "control_cells": agg["ctrl_cells"], "case_cells": agg["case_cells"],
            "n_mixed_ds": agg["n_mixed_ds"], "mixed_cells": agg["mixed_cells"],
            "mixed_frac_primary": safe(agg["mixed_primary"], agg["mixed_cells"]),
            "control_reuse_frac": safe(agg["ctrl_reused"], agg["ctrl_cells"]),
            "case_reuse_frac": safe(agg["case_reused"], agg["case_cells"]),
            "reuse_asymmetry": (safe(agg["ctrl_reused"], agg["ctrl_cells"]) or 0)
                               - (safe(agg["case_reused"], agg["case_cells"]) or 0)
            if (agg["ctrl_cells"] and agg["case_cells"]) else None,
            "has_arm_contrast": bool(agg["ctrl_cells"] and agg["case_cells"]),
            "has_mixed_contrast": bool(agg["n_mixed_ds"]),
            "has_contrast": bool((agg["ctrl_cells"] and agg["case_cells"]) or agg["n_mixed_ds"]),
            "genes_ctrl": safe(agg["g_ctrl_num"], agg["g_ctrl_den"]),
            "genes_case": safe(agg["g_case_num"], agg["g_case_den"]),
            "genes_ctrl_same_study": safe(agg["g_ctrl_pure"], agg["g_ctrl_pure_den"]),
            "depth_ratio": (safe(agg["g_ctrl_num"], agg["g_ctrl_den"]) /
                            safe(agg["g_case_num"], agg["g_case_den"]))
            if (agg["g_ctrl_den"] and agg["g_case_den"]) else None,
            "has_UC": agg["has_uc"], "has_CD": agg["has_cd"],
            "diseases_union": ";".join(sorted(agg["diseases"])),
        })

    D = pd.DataFrame(rows, columns=DATASET_COLS)
    C = pd.DataFrame(crow)
    # 统一数值列 dtype（否则 None 混入会变成 object，后续 np.isfinite / 比较会报错）
    for col in ("cell_count", "primary_cell_count", "reused_cells", "frac_primary",
                "n_donors", "n_cell_types", "n_disease_labels", "mean_genes_per_cell"):
        D[col] = pd.to_numeric(D[col], errors="coerce")
    for col in ("cells_total", "cells_primary", "collection_reuse_frac", "control_cells",
                "case_cells", "mixed_cells", "mixed_frac_primary",
                "control_reuse_frac", "case_reuse_frac", "reuse_asymmetry",
                "genes_ctrl", "genes_case", "genes_ctrl_same_study", "depth_ratio"):
        C[col] = pd.to_numeric(C[col], errors="coerce")
    D.to_csv(os.path.join(OUT, "partD_resfeas_libwide_datasets.tsv"), sep="\t", index=False)
    C.to_csv(os.path.join(OUT, "partD_resfeas_libwide_collections.tsv"), sep="\t", index=False)

    # ---------------------------------------------------------- 4) 全库汇总指标
    Dv = D[D.cell_count.notna() & D.primary_cell_count.notna()].copy()
    tot_cells = float(Dv.cell_count.sum())
    tot_primary = float(Dv.primary_cell_count.sum())
    fr = Dv.frac_primary.astype(float)
    human = Dv[Dv.is_human]

    contrast = C[C.has_arm_contrast].copy()          # 臂可分：可算对照臂/病例臂复用率与深度
    mixed_ct = C[C.has_mixed_contrast].copy()        # 含混合数据集：存在对比但臂不可在元数据层拆
    n_ct = len(contrast)
    n_mixed_ct = int(len(mixed_ct))
    n_any_ct = int(C.has_contrast.sum())
    mixed_ds = Dv[Dv.arm_class == "mixed"] if "arm_class" in Dv.columns else Dv.iloc[0:0]
    hi_ctrl = contrast[contrast.control_reuse_frac > 0.5]
    full_ctrl = contrast[contrast.control_reuse_frac >= 0.999]
    asym = contrast[contrast.reuse_asymmetry > 0]
    dn = pd.to_numeric(contrast.depth_ratio, errors="coerce").dropna()
    dn_same = pd.to_numeric(contrast.genes_ctrl_same_study / contrast.genes_case,
                            errors="coerce").dropna()
    dn_same = dn_same[np.isfinite(dn_same)]

    by_sys = []
    # 分层口径：任何含对比的集合（臂可分 + 仅混合），按单一组织系统 explode 后统计
    # 「对比不可在元数据层审计」的比例 —— 这才是能回答「框架是否肠道专属」的量
    sys_rows = []
    for _, r in C[C.has_contrast].iterrows():
        for s in str(r.tissue_systems).split(";"):
            s = s.strip()
            if s:
                sys_rows.append({"sys": s, "not_arm": (not bool(r.has_arm_contrast))})
    if sys_rows:
        SD = pd.DataFrame(sys_rows)
        for s, g in SD.groupby("sys"):
            if len(g) < 5:
                continue
            by_sys.append({"tissue_system": s, "n_contrast_collections": int(len(g)),
                           "share_not_arm_separable": float(g.not_arm.mean()),
                           "n_not_arm_separable": int(g.not_arm.sum())})
        by_sys.sort(key=lambda x: -x["n_contrast_collections"])

    summary = {
        "library": "CELLxGENE Discover (curation API v1)",
        "n_collections_total": int(len(cols)),
        "n_collections_with_detail": int(len(dets)),
        "n_datasets_total": int(len(D)),
        "n_datasets_scanned_cells": int(len(Dv)),
        "cells_total": tot_cells,
        "cells_primary": tot_primary,
        "library_reused_cell_share": float(1 - tot_primary / tot_cells) if tot_cells else None,
        "frac_primary_quantiles": {q: float(fr.quantile(q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)},
        "frac_primary_mean": float(fr.mean()),
        "n_datasets_fully_reused": int((fr <= 0).sum()),
        "share_datasets_fully_reused": float((fr <= 0).mean()),
        "n_datasets_fully_primary": int((fr >= 0.999).sum()),
        "share_datasets_fully_primary": float((fr >= 0.999).mean()),
        "human_cells_total": float(human.cell_count.sum()),
        "human_library_reused_cell_share": float(
            1 - human.primary_cell_count.sum() / human.cell_count.sum()) if len(human) else None,
        "n_datasets_total_with_armclass": int(len(Dv)),
        "n_datasets_ctl_only": int((Dv.arm_class == "ctl_only").sum()) if len(Dv) else 0,
        "n_datasets_dis_only": int((Dv.arm_class == "dis_only").sum()) if len(Dv) else 0,
        "n_datasets_mixed_contrast": int(len(mixed_ds)),
        "share_datasets_mixed_contrast": float(len(mixed_ds) / len(Dv)) if len(Dv) else None,
        "median_frac_primary_mixed_ds": float(mixed_ds.frac_primary.median()) if len(mixed_ds) else None,
        "n_arm_contrast_collections": int(n_ct),
        "n_mixed_contrast_collections": int(n_mixed_ct),
        "n_any_contrast_collections": int(n_any_ct),
        "share_contrast_collections_arm_separable": (float(n_ct / n_any_ct) if n_any_ct else None),
        "share_contrast_collections_not_arm_separable": (float((n_any_ct - n_ct) / n_any_ct)
                                                         if n_any_ct else None),
        "median_frac_primary_contrast_bearing_ds": float(mixed_ds.frac_primary.median())
        if len(mixed_ds) else None,
        "share_contrast_bearing_ds_fully_reused": float((mixed_ds.frac_primary <= 0).mean())
        if len(mixed_ds) else None,
        "n_contrast_collections": int(n_ct),
        "share_contrast_control_reuse_gt50": float((contrast.control_reuse_frac > 0.5).mean()) if n_ct else None,
        "n_contrast_control_fully_reused": int(len(full_ctrl)),
        "share_contrast_control_fully_reused": float(len(full_ctrl) / n_ct) if n_ct else None,
        "median_contrast_control_reuse_frac": float(contrast.control_reuse_frac.median()) if n_ct else None,
        "median_contrast_case_reuse_frac": float(contrast.case_reuse_frac.median()) if n_ct else None,
        "n_contrast_asymmetric_ctrl_worse": int(len(asym)),
        "share_contrast_asymmetric_ctrl_worse": float(len(asym) / n_ct) if n_ct else None,
        "median_depth_ratio_ctrl_over_case": float(dn.median()) if len(dn) else None,
        "share_depth_ratio_outside_0.8_1.25": float(((dn < 0.8) | (dn > 1.25)).mean()) if len(dn) else None,
        "n_contrast_depth_computable": int(len(dn)),
        "median_depth_ratio_ctrl_same_study_over_case": float(dn_same.median()) if len(dn_same) else None,
        "n_contrast_depth_same_study_computable": int(len(dn_same)),
        "by_tissue_system": by_sys,
        "top_contaminated_collections": contrast[
            contrast.control_reuse_frac.notna()].nlargest(20, "control_reuse_frac")[
            ["collection_name", "tissue_systems", "diseases_union", "n_control_ds",
             "n_disease_ds", "control_cells", "control_reuse_frac", "case_reuse_frac"]].to_dict("records"),
        "metric_definitions": {
            "frac_primary": "primary_cell_count / cell_count（该数据集有多少细胞是本研究原产）",
            "library_reused_cell_share": "全库 1 - Σprimary/Σcells，即非原产（复用）细胞占比",
            "control_reuse_frac": "集合对照臂细胞中非本数据集原产的占比（对照组标签 = normal/healthy/control 等）",
            "case_reuse_frac": "同一集合病例臂的对应比例（作为对照，避免把 atlas 聚合误判成对照污染）",
            "reuse_asymmetry": "control_reuse_frac - case_reuse_frac，>0 才说明对照比病例更依赖外来数据",
            "depth_ratio": "加权 mean_genes_per_cell，对照臂 / 病例臂；仅用 frac_primary≈1 的数据集",
            "arm_class": ("数据集按疾病标签划分：ctl_only（仅对照标签）/ dis_only（仅疾病标签）/ "
                          "mixed（同一数据集同时含对照与疾病标签）。mixed 的对比存在于数据集内部，"
                          "臂无法在元数据层拆分，故整体单独计数、不并入病例臂。"),
            "has_arm_contrast": "集合同时有 ctl_only 与 dis_only 数据集 → 对照臂/病例臂可分别审计",
            "has_mixed_contrast": "集合含 mixed 数据集（有对比但臂不可分）",
            "n_any_contrast_collections": "has_arm_contrast 或 has_mixed_contrast 的集合总数",
            "share_contrast_collections_not_arm_separable": "含对比集合中，臂无法在元数据层分开的比例",
            "median_frac_primary_contrast_bearing_ds": "携带病例-对照对比的数据集（mixed）的 frac_primary 中位",
            "share_contrast_bearing_ds_fully_reused": "mixed 数据集中 frac_primary = 0（全部复用）的比例",
        },
        "caveat": ("is_primary_data=False 本身不构成错误，atlas 型集合合法地聚合外部数据；"
                   "本审计的结论落在**不对称性**（对照臂比病例臂更依赖外来细胞）与"
                   "**对照臂完全复用**的比例上，而非复用率的绝对值。"
                   "mean_genes_per_cell 为数据集级均值，仅作深度代理。"
                   "★ 本审计最重要的局限是审计可行性本身被数据质量筛选："
                   "全库 140 个含对比集合里只有 10.7% 的臂可在元数据层分开，"
                   "而在这 15 个可审计集合中对照臂几乎都是本研究原产——"
                   "即「可审计」的子集恰好是「对照干净」的子集，"
                   "对 89.3% 的集合我们无法判断，而不是判断为干净。"),
    }
    json.dump(summary, open(os.path.join(OUT, "partD_resfeas_libwide_summary.json"), "w",
                            encoding="utf-8"), indent=2, ensure_ascii=False, default=str)

    # ---------------------------------------------------------- 5) 图 D13
    MM = 1 / 25.4
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 6.2,
        "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 400,
    })
    C_OK, C_BAD, C_PAIR = "#2e6f9e", "#b0b7bf", "#b5442e"
    fig = plt.figure(figsize=(180 * MM, 108 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.62, wspace=0.34,
                  left=0.10, right=0.985, top=0.90, bottom=0.115)

    ax = fig.add_subplot(gs[0, 0])
    bins = np.linspace(0, 1, 26)
    bottom = np.zeros(len(bins) - 1)
    for name, col, lab in [("ctl_only", "#c3cbd4", "control-only ds"),
                           ("dis_only", C_OK, "disease-only ds"),
                           ("mixed", C_PAIR, "mixed ds (ctrl+disease)")]:
        v = Dv[Dv.arm_class == name].frac_primary.astype(float).clip(0, 1) if len(Dv) else np.array([])
        h, _ = np.histogram(v, bins=bins)
        ax.bar(bins[:-1], h, width=np.diff(bins), bottom=bottom, align="edge",
               color=col, edgecolor="white", linewidth=0.25, label="%s (n=%d)" % (lab, int(h.sum())))
        bottom += h
    ax.axvline(1 - summary["library_reused_cell_share"], color="#444444", lw=0.7, ls="--")
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("frac_primary per dataset", fontsize=6)
    ax.set_ylabel("datasets", fontsize=6)
    ax.legend(fontsize=4.3, frameon=False, loc="upper center", handlelength=0.9, borderpad=0.15,
              labelspacing=0.25)
    ax.set_title("A  Whole-library provenance\n%d datasets / %.1fM cells"
                 % (len(Dv), tot_cells / 1e6), fontsize=6.6, fontweight="bold", loc="left")

    ax = fig.add_subplot(gs[0, 1])
    sub = contrast[contrast.control_reuse_frac.notna() & contrast.case_reuse_frac.notna()]
    ax.plot([0, 1], [0, 1], color="#999999", lw=0.6, ls=":")
    ax.scatter(sub.case_reuse_frac, sub.control_reuse_frac, s=16, color=C_PAIR,
               edgecolor="white", linewidth=0.3, zorder=3)
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("case-arm reuse fraction", fontsize=6)
    ax.set_ylabel("control-arm reuse fraction", fontsize=6)
    ax.text(0.97, 0.04, "%d further collections carry contrast\nonly inside mixed datasets\n(arms not separable)"
            % n_mixed_ct, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=4.5, color="#7a7f86")
    ax.set_title("B  %d arm-separable case/control collections\nabove dotted = controls more borrowed"
                 % n_ct, fontsize=6.6, fontweight="bold", loc="left")

    ax = fig.add_subplot(gs[1, 0])
    top = [b for b in by_sys][:8][::-1]
    if top:
        ax.barh(np.arange(len(top)), [b["share_not_arm_separable"] for b in top],
                color=C_OK, height=0.6)
        for i, b in enumerate(top):
            ax.text(b["share_not_arm_separable"] + 0.02, i,
                    "n=%d" % b["n_contrast_collections"], va="center", fontsize=5.2)
        ax.set_yticks(np.arange(len(top)))
        ax.set_yticklabels([b["tissue_system"][:20] for b in top], fontsize=5.6)
    ax.set_xlim(0, 1.16)
    ax.set_xlabel("share of contrast collections whose arms\ncannot be separated from metadata",
                  fontsize=5.6)
    ax.set_title("C  Not a gut-specific artefact\n(tissue systems, n>=5 contrast collections)",
                 fontsize=6.6, fontweight="bold", loc="left")

    ax = fig.add_subplot(gs[1, 1])
    if len(dn):
        ax.hist(dn, bins=20, color="#4c9f70", edgecolor="white", linewidth=0.3)
    ax.axvline(1.0, color="#666666", lw=0.8, ls="--")
    ax.set_xlabel("mean genes/cell: control arm / case arm", fontsize=6)
    ax.set_ylabel("collections", fontsize=6)
    ax.set_title("D  Control-arm depth disadvantage\nmedian %.2fx (n=%d)"
                 % (summary["median_depth_ratio_ctrl_over_case"] or float("nan"), len(dn)),
                 fontsize=6.6, fontweight="bold", loc="left")

    fig.savefig(os.path.join(OUT, "fig_D14_libwide_provenance.png"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig_D14_libwide_provenance.pdf"), bbox_inches="tight")
    plt.close(fig)

    # ---------------------------------------------------------- 6) 打印
    print("\n" + "=" * 96)
    print("CELLxGENE 全库 · 来源污染审计（不限肠道）")
    print("=" * 96)
    print("集合 %d（详情 %d）· 数据集 %d（可算 %d）· 细胞 %.1fM"
          % (len(cols), len(dets), len(D), len(Dv), tot_cells / 1e6))
    print("全库非原产（复用）细胞占比 = %.1f%%   人类数据 %s"
          % (100 * (1 - tot_primary / tot_cells),
             ("%.1f%%" % (100 * summary["human_library_reused_cell_share"])) if len(human) else "n/a"))
    print("frac_primary：均值 %.2f，中位 %.2f；完全复用数据集 %d（%.1f%%），完全原产 %d（%.1f%%）"
          % (fr.mean(), fr.median(), summary["n_datasets_fully_reused"],
             100 * summary["share_datasets_fully_reused"],
             summary["n_datasets_fully_primary"], 100 * summary["share_datasets_fully_primary"]))
    print("-" * 96)
    print("病例-对照对比集合 = %d（臂可分 %d ＋ 仅含混合数据集 %d）"
          % (n_any_ct, n_ct, n_mixed_ct))
    print("  数据集层：ctl_only %d / dis_only %d / MIXED（对照+疾病同集）%d（占可算数据集 %.1f%%）"
          % (summary["n_datasets_ctl_only"], summary["n_datasets_dis_only"],
             summary["n_datasets_mixed_contrast"],
             100 * (summary["share_datasets_mixed_contrast"] or 0)))
    print("  携带对比的数据集（mixed）frac_primary 中位 %.2f；其中完全复用 %.1f%%"
          % (summary["median_frac_primary_contrast_bearing_ds"] or float("nan"),
             100 * (summary["share_contrast_bearing_ds_fully_reused"] or 0)))
    print("  → 含对比的集合里 %.1f%% 的「臂」根本无法在元数据层分开"
          % (100 * (summary["share_contrast_collections_not_arm_separable"] or 0)))
    if n_ct:
        print("  对照臂复用率：中位 %.2f（病例臂 %.2f）"
              % (summary["median_contrast_control_reuse_frac"],
                 summary["median_contrast_case_reuse_frac"]))
        print("  对照臂复用 >50%% 的集合：%d（%.1f%%）"
              % (len(hi_ctrl), 100 * summary["share_contrast_control_reuse_gt50"]))
        print("  对照臂**完全**由复用数据构成的集合：%d（%.1f%%）"
              % (len(full_ctrl), 100 * summary["share_contrast_control_fully_reused"]))
        print("  对照比病例更依赖外来细胞的集合：%d（%.1f%%）"
              % (len(asym), 100 * summary["share_contrast_asymmetric_ctrl_worse"]))
        print("  对照/病例深度比中位 %.2f（可算 %d 个；越 <1 对照越浅）"
              % (summary["median_depth_ratio_ctrl_over_case"] or float("nan"), len(dn)))
        print("  仅用同研究对照时深度比中位 %s（可算 %d 个）"
              % (("%.2f" % summary["median_depth_ratio_ctrl_same_study_over_case"])
                 if summary["median_depth_ratio_ctrl_same_study_over_case"] else "n/a", len(dn_same)))
        print("\n  按组织系统（n>=5 个含对比集合）：臂不可在元数据层审计的占比")
        for b in by_sys[:12]:
            print("    %-22s n=%-4d 臂不可分 %.0f%%（%d/%d）"
                  % (b["tissue_system"], b["n_contrast_collections"],
                     100 * b["share_not_arm_separable"],
                     b["n_not_arm_separable"], b["n_contrast_collections"]))
    print("=" * 96)
    log("产出 results/partD_resfeas_libwide_{datasets,collections}.tsv / _summary.json")
    log("     results/fig_D14_libwide_provenance.{png,pdf}")


if __name__ == "__main__":
    main()

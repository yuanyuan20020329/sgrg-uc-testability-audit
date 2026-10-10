# -*- coding: utf-8 -*-
"""Part D / 公共资源可用性审计 —— CELLxGENE 侧（全库）

为什么做
--------
Part D 的终判定是「现有公共单细胞资源在结构上无法做出有检验力、批次受控的检验」。
审稿人必问的第一句是「你真的找遍了吗」。本脚本把这句话变成**可审计的枚举**：
遍历 CELLxGENE Discover 全库（397 集合 / 2238 数据集），筛出全部**肠道组织**数据集，
逐一记录其设计参数（细胞数、原产细胞数、供体清单、疾病标签、assay、是否空间），
并判定每个集合/数据集**能否构成一个「病例 + 同源对照」的受控比较**。

★ 关键实现点（决定这步能不能成立）
  CELLxGENE 的 `primary_cell_count` 与 `donor_id` / `is_primary_data` 在**集合详情**端点里就有，
  不需要下载 h5ad。也就是说来源审计（ds4 那次发现「对照 94.4% 是复用外部数据」）可以
  **在元数据层面对全库零成本复现**。这是本脚本相对 partD_batch_audit.py（读 h5ad obs）的根本升级：
  从「四套数据」扩到「全库」。

判定口径
--------
  · 组内供体数：从 donor_id 名单按命名前缀归类（HC-/UC-/CD-/PSC- 等），无法归类时记为 unclassified；
  · 同源对照：同一**集合**内同时存在病例与对照，且对照细胞的主体是 is_primary_data=True 的
    —— 注意 is_primary_data 是**数据集级**的 True/False 列表，需配合 primary_cell_count 解读；
  · 检验力门槛：同 donor 配对设计需 ≥6 对（双侧符号翻转最小 p=2/2^n 才 <0.05），
    独立两组设计需每组 ≥10 例才有名义检验力（df≈18）。

输出
----
  results/partD_resfeas_cxg_datasets.tsv     每个肠道数据集一行
  results/partD_resfeas_cxg_collections.tsv  每个含肠道数据的集合一行
  results/partD_resfeas_cxg_summary.json
  data/raw/cxg_coll_detail/<cid>.json        集合详情缓存（重跑不重复请求）
用法
----
  python partD_resfeas_cxg.py [--refresh]
"""
import argparse
import json
import os
import re
import ssl
import socket
import time
import urllib.request

import pandas as pd

socket.setdefaulttimeout(120)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

API = "https://api.cellxgene.cziscience.com/curation/v1"
CACHE = "data/raw/cxg_coll_detail"
OUT = "results"
T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)

GUT = re.compile(r"\b(colon|colonic|intestin\w*|rectum|rectal|ileum|ileal|cecum|cecal|"
                 r"duoden\w*|jejun\w*|anal|sigmoid|ascending colon|descending colon|gut)\b", re.I)
IBD = re.compile(r"ulcerative colitis|crohn|inflammatory bowel|colitis|IBD", re.I)

# donor_id 前缀归类（顺序敏感：先具体后笼统）
DONOR_PATTERNS = [
    ("PSC", re.compile(r"^(psc|psc[-_ ]?uc)\b", re.I)),
    ("UC",  re.compile(r"^(uc|ulcer\w*|u)[-_ .]?\d", re.I)),
    ("CD",  re.compile(r"^(cd|crohn\w*)[-_ .]?\d", re.I)),
    ("CTL", re.compile(r"^(hc|ctrl|ctl|con|cont|hv|healthy|normal|norm|n|ld|non[-_ ]?ibd|ctl)[-_ .]?\d", re.I)),
    ("IBD", re.compile(r"^(ibd)[-_ .]?\d", re.I)),
]


def get_json(url, tries=4):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0",
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            last = e
            time.sleep(1.5 * (k + 1))
    raise last


def classify_donor(d):
    s = str(d).strip()
    for grp, pat in DONOR_PATTERNS:
        if pat.search(s):
            return grp
    return "UNK"


def labels(field):
    return [x.get("label") for x in (field or []) if x.get("label")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="忽略集合详情缓存")
    a = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)

    # ---- 1) 全库列表（slim）
    slim_path = "data/raw/cxg_all_collections.json"
    if os.path.exists(slim_path) and not a.refresh:
        cols = json.load(open(slim_path, encoding="utf-8"))
        log("载入缓存全库列表 %d 集合" % len(cols))
    else:
        cols = get_json(API + "/collections")
        json.dump(cols, open(slim_path, "w", encoding="utf-8"), ensure_ascii=False)
        log("下载全库列表 %d 集合" % len(cols))

    # ---- 2) 候选集合：含肠道组织数据集，或疾病标签含 colitis
    cand = []
    for c in cols:
        hit_gut = hit_ibd = False
        for d in c.get("datasets", []):
            if GUT.search(" ".join(labels(d.get("tissue")))):
                hit_gut = True
            if IBD.search(" ".join(labels(d.get("disease")))) or IBD.search(c.get("name") or ""):
                hit_ibd = True
        if hit_gut or hit_ibd:
            cand.append(c)
    log("候选集合（含肠道组织 或 含 colitis 疾病/名称）= %d / %d" % (len(cand), len(cols)))

    # ---- 3) 抓集合详情（带缓存）
    rows, coll_rows, n_fetch = [], [], 0
    for i, c in enumerate(cand, 1):
        cid = c["collection_id"]
        fp = os.path.join(CACHE, cid + ".json")
        if os.path.exists(fp) and not a.refresh:
            det = json.load(open(fp, encoding="utf-8"))
        else:
            try:
                det = get_json(API + "/collections/" + cid)
                json.dump(det, open(fp, "w", encoding="utf-8"), ensure_ascii=False)
                n_fetch += 1
            except Exception as e:
                log("  抓取失败 %s %s: %s" % (cid[:8], c.get("name", "")[:40], e))
                continue
        cn = det.get("name") or c.get("name")
        doi = det.get("doi")
        ds = det.get("datasets") or []
        coll_gut = 0
        coll_dis, coll_donors = set(), set()
        for d in ds:
            tis = labels(d.get("tissue"))
            dis = labels(d.get("disease"))
            gut = [t for t in tis if GUT.search(t)]
            if not gut:
                continue
            coll_gut += 1
            donors = [str(x) for x in (d.get("donor_id") or [])]
            cls = {g: 0 for g in ("UC", "CD", "CTL", "PSC", "IBD", "UNK")}
            for x in donors:
                cls[classify_donor(x)] += 1
            cc = d.get("cell_count")
            pcc = d.get("primary_cell_count")
            row = {
                "collection_id": cid, "collection_name": cn, "doi": doi,
                "dataset_id": d.get("dataset_id"), "dataset_title": d.get("title"),
                "cell_count": cc, "primary_cell_count": pcc,
                "frac_primary": (pcc / cc) if (cc and pcc is not None) else None,
                "n_donors": len(donors),
                "n_UC_donors": cls["UC"], "n_CD_donors": cls["CD"], "n_CTL_donors": cls["CTL"],
                "n_PSC_donors": cls["PSC"], "n_IBD_donors": cls["IBD"], "n_UNK_donors": cls["UNK"],
                "donors": ";".join(donors),
                "diseases": ";".join(dis), "tissues": ";".join(tis),
                "gut_tissues": ";".join(sorted(set(gut))),
                "assays": ";".join(labels(d.get("assay"))),
                "suspension_type": ";".join(d.get("suspension_type") or []),
                "spatial": bool(d.get("spatial")),
                "is_primary_data": ";".join(str(x) for x in (d.get("is_primary_data") or [])),
                "organism": ";".join(labels(d.get("organism"))),
            }
            dl = set(x.lower() for x in dis)
            row["has_UC"] = any("ulcerative colitis" in x for x in dl)
            row["has_CD"] = any("crohn" in x for x in dl)
            row["has_normal"] = any(x in ("normal", "healthy", "control") or "normal" in x for x in dl)
            row["has_PSC"] = any("sclerosing cholangitis" in x for x in dl)
            rows.append(row)
            coll_dis |= set(dis)
            coll_donors |= set(donors)

        cd = set(x.lower() for x in coll_dis)
        has_uc = any("ulcerative colitis" in x for x in cd)
        has_norm = any("normal" in x or "healthy" in x or "control" in x for x in cd)
        coll_rows.append({
            "collection_id": cid, "collection_name": cn, "doi": doi,
            "published_at": det.get("published_at"),
            "n_datasets_gut": coll_gut, "n_datasets_total": len(ds),
            "n_donors_union": len(coll_donors),
            "collection_has_UC": has_uc, "collection_has_normal": has_norm,
            "collection_has_CD": any("crohn" in x for x in cd),
            "collection_has_PSC": any("sclerosing cholangitis" in x for x in cd),
            "diseases_union": ";".join(sorted(set(labels(
                [{"label": x} for x in coll_dis])))),
            "supports_controlled_UC_comparison": bool(has_uc and has_norm),
        })
        if i % 10 == 0:
            log("  进度 %d/%d（本轮新抓 %d）" % (i, len(cand), n_fetch))
    log("集合详情抓取完成：本轮新抓 %d 个" % n_fetch)

    D = pd.DataFrame(rows)
    C = pd.DataFrame(coll_rows)
    D = D.sort_values(["has_UC", "cell_count"], ascending=[False, False])
    C = C.sort_values(["supports_controlled_UC_comparison", "n_datasets_gut"],
                      ascending=[False, False])
    D.to_csv(os.path.join(OUT, "partD_resfeas_cxg_datasets.tsv"), sep="\t", index=False)
    C.to_csv(os.path.join(OUT, "partD_resfeas_cxg_collections.tsv"), sep="\t", index=False)

    summ = {
        "n_collections_total": len(cols),
        "n_collections_scan_detail": len(cand),
        "n_collections_with_gut_data": int(len(C)),
        "n_datasets_gut": int(len(D)),
        "n_datasets_gut_has_UC": int(D.has_UC.sum()),
        "n_datasets_gut_has_CD": int(D.has_CD.sum()),
        "n_datasets_gut_has_normal": int(D.has_normal.sum()),
        "n_collections_with_UC": int(C.collection_has_UC.sum()),
        "n_collections_UC_and_normal": int(C.supports_controlled_UC_comparison.sum()),
        "n_collections_UC_and_normal_and_spatial": int(
            C[C.supports_controlled_UC_comparison].merge(
                D[D.spatial][["collection_id"]].drop_duplicates(), on="collection_id").shape[0]),
        "n_datasets_spatial_gut": int(D.spatial.sum()),
        "note_primary_cell_count": (
            "primary_cell_count / donor_id / is_primary_data 均来自 CELLxGENE 集合详情端点，"
            "可在元数据层面复现来源审计，无需下载 h5ad"),
        "power_floors": {"paired_signed_flip_min_p_at_n": {str(n): 2.0 / (2 ** n) for n in (4, 6, 7, 8, 10, 15)},
                         "independent_two_group_rule_of_thumb": 10},
    }
    json.dump(summ, open(os.path.join(OUT, "partD_resfeas_cxg_summary.json"), "w",
                         encoding="utf-8"), indent=2, ensure_ascii=False, default=str)

    print("\n" + "=" * 92)
    print("CELLxGENE 全库 · 肠道数据集可用性")
    print("=" * 92)
    print("集合 %d → 扫描详情 %d → 含肠道数据 %d；肠道数据集 %d（UC %d / CD %d / 含正常 %d / 空间 %d）"
          % (summ["n_collections_total"], summ["n_collections_scan_detail"],
             summ["n_collections_with_gut_data"], summ["n_datasets_gut"],
             summ["n_datasets_gut_has_UC"], summ["n_datasets_gut_has_CD"],
             summ["n_datasets_gut_has_normal"], summ["n_datasets_spatial_gut"]))
    print("含 UC 的集合 %d；同时含 UC 与正常（可做受控比较）的集合 %d"
          % (summ["n_collections_with_UC"], summ["n_collections_UC_and_normal"]))
    if len(C):
        cols_show = ["collection_name", "n_datasets_gut", "n_donors_union",
                     "collection_has_UC", "collection_has_normal"]
        print("\n支持受控 UC 比较的集合：")
        sub = C[C.supports_controlled_UC_comparison]
        if len(sub):
            print(sub[cols_show].to_string(index=False, max_colwidth=58))
        else:
            print("  （无）")
        print("\n全部含 UC 的集合：")
        print(C[C.collection_has_UC][cols_show].to_string(index=False, max_colwidth=58))
    log("完成 %.1f s" % (time.time() - T0))
    log("产出 results/partD_resfeas_cxg_datasets.tsv / _collections.tsv / _summary.json")


if __name__ == "__main__":
    main()

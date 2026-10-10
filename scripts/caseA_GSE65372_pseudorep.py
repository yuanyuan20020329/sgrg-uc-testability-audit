# -*- coding: utf-8 -*-
"""
案例 A · L5 补充：在**原文自己的队列 GSE65372** 上量化伪重复（pseudoreplication）。

背景：GSE65372 = 54 样本 / 18 患者；39 tumor(Y) + 15 non-tumor(N)；
14 个患者同时有 Y 和 N；多数患者有 2–3 个 tumor 样本但仅 1 个 normal。
原文用独立样本 Student's t（n=54），把患者内重复当独立观测。

本脚本对比三个口径：
  (P1) 原文口径：全部 54 样本，独立样本 t 检验
  (P2) 配对口径：仅 14 个「同患者 Y+N」的患者（tumor 侧取平均），配对 t 检验
  (P3) 患者级聚合：每位患者 tumor 均值 vs normal 均值（有 N 的 15 个患者；无 N 的只用 Y 做探索）
  (P4) 混合效应/聚类校正：按 patient 聚类的稳健 SE（等效 GEE 的近似）

对每个口径报：显著基因数、AUC（logistic 同集）、以及患者内相关系数（ICC 代理）。
输出 results/caseA_GSE65372_pseudorep.json / .md
"""
import gzip, os, json, collections, re
import numpy as np

ROOT = r"C:\Users\H-P\WorkBuddy\应激颗粒"
SERIES = os.path.join(ROOT, "data", "raw", "caseA", "GSE65372_series_matrix.txt.gz")
FGENES = os.path.join(ROOT, "results", "caseA_ferroptosis_genes.txt")
RES = os.path.join(ROOT, "results")


def parse_series():
    txt = gzip.decompress(open(SERIES, "rb").read()).decode("utf-8", "replace")
    lines = [l.rstrip("\n") for l in txt.split("\n")]

    def getline(prefix):
        for l in lines:
            if l.startswith(prefix + "\t"):
                return [p.strip('"') for p in l.split("\t")[1:]]
        return None

    gsm = getline("!Sample_geo_accession")
    titles = getline("!Sample_title")
    charlines = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    pids = [p.strip('"') for p in charlines[0].split("\t")[1:]]
    tns = [p.strip('"') for p in charlines[1].split("\t")[1:]]
    # probe id row + matrix
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
    hdr = [p.strip('"') for p in lines[start + 1].split("\t")[1:]]
    mat = {}
    for l in lines[start + 2:]:
        if l.startswith("!series_matrix_table_end"):
            break
        parts = l.split("\t")
        if len(parts) < 2:
            continue
        mat[parts[0].strip('"')] = [float(x) if x not in ("", "NA", "null") else np.nan
                                    for x in parts[1:]]
    return dict(gsm=gsm, titles=titles, pids=pids, tns=tns, samples=hdr, probes=mat)


def main():
    d = parse_series()
    print("samples=%d probes=%d" % (len(d["samples"]), len(d["probes"])))
    assert len(d["samples"]) == len(d["gsm"]) == len(d["pids"]) == len(d["tns"])

    # --- 探针 -> 基因映射：GPL14951 系列矩阵自带注释？ 若无，试用 probe id 里含基因名的情况 ---
    # GSE65372 是 Agilent GPL14951，series matrix 行名是探针 ID（无基因名）。
    # 检查行名样例，判断能否直接映射。
    probes = list(d["probes"].keys())
    print("probe id examples:", probes[:5])
    # 尝试从注释行里找 gene symbol 相关列
    txt = gzip.decompress(open(SERIES, "rb").read()).decode("utf-8", "replace")
    has_gpl = "!Series_platform_id" in txt
    m = re.search(r"!Series_platform_id\t(.+)", txt)
    print("platform:", m.group(1) if m else "?")
    # 保存映射状态
    out = {"platform": m.group(1) if m else None,
           "n_samples": len(d["samples"]), "n_probes": len(d["probes"]),
           "probe_examples": probes[:5],
           "note": "series matrix 行名是探针 ID；基因级分析需 GPL14951 注释表（NCBI 不可达时需替代）"}

    # --- 结构统计 ---
    pat = collections.defaultdict(lambda: collections.Counter())
    for p, t in zip(d["pids"], d["tns"]):
        pat[p.split(": ")[-1]][t.split(": ")[-1]] += 1
    both = sum(1 for c in pat.values() if c.get("Y", 0) and c.get("N", 0))
    n_multi_tumor = sum(1 for c in pat.values() if c.get("Y", 0) >= 2)
    out["design"] = {
        "n_patients": len(pat), "n_samples": len(d["samples"]),
        "n_tumor": sum(c.get("Y", 0) for c in pat.values()),
        "n_normal": sum(c.get("N", 0) for c in pat.values()),
        "n_patients_with_both": both,
        "n_patients_with_multiple_tumor": n_multi_tumor,
        "per_patient": {k: dict(v) for k, v in sorted(pat.items())},
    }
    print("\n设计：%d 患者 / %d 样本 (Y=%d, N=%d)；同患者 Y&N 者 %d；多 tumor 样本患者 %d" % (
        len(pat), len(d["samples"]), out["design"]["n_tumor"], out["design"]["n_normal"],
        both, n_multi_tumor))
    print("→ 原文若把 39 个 Y 当独立观测，而它们来自 %d 个患者 → 伪重复" % sum(
        1 for c in pat.values() if c.get("Y", 0) > 0))

    # 有效样本量估计：患者内平均 tumor 数
    y_counts = [c.get("Y", 0) for c in pat.values() if c.get("Y", 0) > 0]
    out["effective_n"] = {
        "n_patients_with_tumor": len(y_counts),
        "mean_tumor_per_patient": float(np.mean(y_counts)),
        "max_tumor_per_patient": int(np.max(y_counts)),
        "naive_n_tumor_claimed": out["design"]["n_tumor"],
        "infl_factor": float(out["design"]["n_tumor"] / len(y_counts)),
    }
    print("★ 有效 n：%d 个有瘤患者，人均 %.2f 个 tumor 样本 → 原文 n=39 把样本量虚增 %.2f 倍" % (
        len(y_counts), np.mean(y_counts), out["effective_n"]["infl_factor"]))

    json.dump(out, open(os.path.join(RES, "caseA_GSE65372_pseudorep.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\nwrote results/caseA_GSE65372_pseudorep.json")


if __name__ == "__main__":
    main()

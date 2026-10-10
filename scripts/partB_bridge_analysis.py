# -*- coding: utf-8 -*-
"""
Part B · Step 3 — 脑肠桥梁交叉分析
三重交集 + 面板富集 + 方向一致性 + 跨队列复现 + 证据叠合表。

输出：
  results/partB_bridge_results.json      所有核心统计
  results/partB_evidence_stack.tsv       证据叠合表（基因 × 证据流）
"""
import os, json, csv, math
import numpy as np
from scipy import stats

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT  = os.path.join(BASE, "results")

def log(*a): print(*a, flush=True)

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

# ------------------------------------------------------------------ 读入
# 1) MDD 各队列 by-gene 抓取结果
mdd_rows = []
with open(os.path.join(OUT, "partB_sgrg_mdd_de.tsv"), encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        mdd_rows.append(row)

# 子集映射（rs → 区域标签）
SUBSET_LABEL = {
    "493107": "nucleus accumbens",
    "493108": "anterior cingulate cortex",
    "493109": "dorsolateral prefrontal cortex",
    "456122": "striatum",
    "456123": "hippocampal formation",
    "456124": "dorsolateral prefrontal cortex",
    "540049": "whole blood",
    "540050": "whole blood (treatment)",   # 治疗因子，不用于疾病对照
    "540051": "whole blood (subset)",
    "644372": "colon biopsy (Gemma stored)",
}

# 2) GSE19738 全基因组（权威：单对照，探针级 q 即 MDD FDR）
universe = set()
mdd_de_genes = set()
sgrg_in_19738 = {}   # gene -> dict(best probe p/q/lfc)
sgrg_syms = set()
with open(os.path.join(OUT, "sgrg_gene_set_clean.tsv"), encoding="utf-8") as f:
    next(f)
    for line in f:
        sgrg_syms.add(line.split("\t")[0])

with open(os.path.join(OUT, "partB_gse19738_genomewide.tsv"), encoding="utf-8") as f:
    r = csv.reader(f, delimiter="\t")
    header = next(r)
    for row in r:
        probe, syms, p, q, lfc, tstat = row
        ss = [s for s in syms.split("|") if s and s != "NA"]
        for s in ss:
            universe.add(s)
        qq = fnum(q)
        if qq is not None and qq <= 0.05:
            for s in ss:
                mdd_de_genes.add(s)
        hit = [s for s in ss if s in sgrg_syms]
        if hit:
            g = hit[0]
            pp, ll = fnum(p), fnum(lfc)
            cur = sgrg_in_19738.get(g)
            if cur is None or (pp is not None and (cur["p"] is None or pp < cur["p"])):
                sgrg_in_19738[g] = dict(probe=probe, p=pp, q=qq, lfc=ll)

# 3) UC 侧（自算，Part A Step1）
uc = {}
with open(os.path.join(OUT, "partA_GSE38713_SGRG_DE.tsv"), encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        uc[row["GeneSymbol"]] = dict(
            lfc=fnum(row["log2FC_inflamed_vs_ctrl"]),
            p=fnum(row["p_inflamed"]),
            fdr=fnum(row["FDR_inflamed"]),
        )

# 4) UC 配对验证（Part A Step3, GSE22619 双生子）
paired = {}
with open(os.path.join(OUT, "partA_step3_paired_DE.tsv"), encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        paired[row["GeneSymbol"]] = dict(
            dz=fnum(row["cohen_dz"]),
            p=fnum(row["wilcoxon_p"]),
            diff=fnum(row["mean_diff_UC_minus_twin"]),
        )

# 5) Part C MR（每基因跨组织取最小 IVW p 及其 q）
mr = {}
with open(os.path.join(OUT, "partC_mr_geneset.tsv"), encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        g = row["gene"]
        p = fnum(row["IVW_p"]); q = fnum(row["IVW_q"])
        if p is None:
            continue
        if g not in mr or p < mr[g]["IVW_p"]:
            mr[g] = dict(tissue=row["tissue"], nsnp=fnum(row["nsnp"]),
                         IVW_p=p, IVW_q=q, IVW_b=fnum(row["IVW_b"]),
                         egger_p=fnum(row["Egger_p"]))

log(f"universe(GSE19738)={len(universe)}  MDD-DEG={len(mdd_de_genes)}")
log(f"SGRG in 19738={len(sgrg_in_19738)}  UC DE rows={len(uc)}  paired={len(paired)}  MR genes={len(mr)}")

# ------------------------------------------------------------------ 1) SGRG 面板富集（MDD 外周血）
panel_universe = [g for g in sgrg_syms if g in universe]
panel_de = [g for g in panel_universe if g in mdd_de_genes]
N, K, n = len(universe), len(mdd_de_genes), len(panel_universe)
k = len(panel_de)
# 超几何：面板里 DE 基因数 ≥ 观测 k 的概率
p_over = stats.hypergeom.sf(k - 1, N, K, n)
exp_k = n * K / N
OR = (k / (n - k)) / ((K - k) / (N - n - (K - k))) if (n - k) > 0 and (K - k) > 0 else float("nan")
log(f"[enrichment] k={k}/{n} observed vs expected {exp_k:.2f} | OR={OR:.2f} | p={p_over:.3g}")
enrichment = dict(N=N, K=K, n=n, k=k, expected=exp_k, OR=OR, p=p_over,
                  panel_de_genes=sorted(panel_de))

# 等规模随机面板分布（对齐 Part A Step2b 的口径）
rng = np.random.default_rng(20261009)
universe_list = sorted(universe)
null = []
for _ in range(2000):
    sel = rng.choice(len(universe_list), size=n, replace=False)
    null.append(sum(1 for i in sel if universe_list[i] in mdd_de_genes))
null = np.array(null)
pct = float((null < k).mean() * 100)
log(f"[enrichment] null median={np.median(null):.1f}, p95={np.percentile(null,95):.1f}, SGRG percentile={pct:.1f}")
enrichment.update(null_median=float(np.median(null)), null_p95=float(np.percentile(null, 95)),
                  percentile=pct)

# ------------------------------------------------------------------ 2) UC×MDD 方向一致性（SGRG 面板内可比基因）
both = []
for g in sgrg_in_19738:
    if g in uc and sgrg_in_19738[g]["lfc"] is not None and uc[g]["lfc"] is not None:
        both.append((g, uc[g]["lfc"], sgrg_in_19738[g]["lfc"]))
lfc_uc = np.array([b[1] for b in both]); lfc_mdd = np.array([b[2] for b in both])
rho, rho_p = stats.spearmanr(lfc_uc, lfc_mdd)
same_sign = int((np.sign(lfc_uc) == np.sign(lfc_mdd)).sum())
n_cmp = len(both)
binom_p = stats.binomtest(same_sign, n_cmp, 0.5).pvalue
log(f"[concordance] n={n_cmp} same_sign={same_sign} Spearman rho={rho:.3f} p={rho_p:.3g} | binomial p={binom_p:.3g}")

# 名义显著双侧同向
nom = [(g, u, m) for (g, u, m) in both
       if uc.get(g, {}).get("p", 1) < 0.05 and (sgrg_in_19738[g]["p"] or 1) < 0.05]
nom_same = sum(1 for (g, u, m) in nom if np.sign(u) == np.sign(m))
log(f"[concordance] nominal-both n={len(nom)}, same sign {nom_same}")

# ------------------------------------------------------------------ 3) 三重交集：SGRG ∩ UC-DE(FDR<0.05) ∩ MDD-DE(q<0.05)
uc_de = {g for g, v in uc.items() if v["fdr"] is not None and v["fdr"] < 0.05}
mdd_de_panel = {g for g in sgrg_in_19738 if sgrg_in_19738[g]["q"] is not None and sgrg_in_19738[g]["q"] <= 0.05}
triple = sorted(set(uc_de) & mdd_de_panel)
triple_info = []
for g in triple:
    triple_info.append(dict(
        gene=g,
        uc_lfc=uc[g]["lfc"], uc_fdr=uc[g]["fdr"],
        mdd_lfc=sgrg_in_19738[g]["lfc"], mdd_q=sgrg_in_19738[g]["q"],
        same_dir=bool(np.sign(uc[g]["lfc"]) == np.sign(sgrg_in_19738[g]["lfc"])),
        paired_p=paired.get(g, {}).get("p"),
        mr_p=mr.get(g, {}).get("IVW_p"),
    ))
log(f"[triple] UC-DE={len(uc_de)} MDD-DE(panel)={len(mdd_de_panel)} triple={triple}")

# 宽口径（UC 名义 p<0.05 且 MDD q<0.05）
uc_nom = {g for g, v in uc.items() if v["p"] is not None and v["p"] < 0.05}
triple_nom = sorted(uc_nom & mdd_de_panel)
log(f"[triple-nominal] {len(triple_nom)} genes: {triple_nom}")

# ------------------------------------------------------------------ 4) 脑区复现（80655 / 53987 各子集）
# 折叠 per gene×subset：min-p 探针
sub_data = {}
for row in mdd_rows:
    ds = row["dataset"]
    if ds not in ("GSE80655", "GSE53987"):
        continue
    rs = row["rs"]
    g = row["gene"]
    p = fnum(row["p"]); lfc = fnum(row["lfc"])
    if p is None:
        continue
    key = (ds, rs, g)
    if key not in sub_data or p < sub_data[key]["p"]:
        sub_data[key] = dict(p=p, lfc=lfc, probe=row["probe"])

brain_repl = {}
for g in sorted(sgrg_syms):
    entry = {}
    for (ds, rs, gg), v in sub_data.items():
        if gg != g:
            continue
        label = f"{ds}:{SUBSET_LABEL.get(rs, rs)}"
        entry[label] = v
    if entry:
        brain_repl[g] = entry

# 血清 DE 基因在脑区方向一致比例
blood_de = {g: sgrg_in_19738[g]["lfc"] for g in mdd_de_panel}
concord = {}
for label in sorted({f'{ds}:{SUBSET_LABEL.get(rs, rs)}' for (ds, rs, _) in sub_data}):
    hits = []
    for g, bl in blood_de.items():
        v = brain_repl.get(g, {}).get(label)
        if v and bl is not None and v["lfc"] is not None:
            hits.append(np.sign(bl) == np.sign(v["lfc"]))
    if hits:
        concord[label] = dict(n=len(hits), same=int(np.sum(hits)))
log(f"[brain-replication of blood DEGs] {concord}")

# ------------------------------------------------------------------ 5) 证据叠合表
stack_cols = ["gene",
              "UC_lfc_GSE38713", "UC_FDR",
              "UCpaired_dz_GSE22619", "UCpaired_p",
              "MDDblood_lfc_GSE19738", "MDDblood_p", "MDDblood_q",
              "nAcc_lfc", "nAcc_p", "ACC_lfc", "ACC_p", "dlPFC80655_lfc", "dlPFC80655_p",
              "STR_lfc", "STR_p", "HPC_lfc", "HPC_p", "dlPFC53987_lfc", "dlPFC53987_p",
              "MR_tissue", "MR_IVW_p", "MR_IVW_q"]
label2col = {
    "GSE80655:nucleus accumbens": ("nAcc_lfc", "nAcc_p"),
    "GSE80655:anterior cingulate cortex": ("ACC_lfc", "ACC_p"),
    "GSE80655:dorsolateral prefrontal cortex": ("dlPFC80655_lfc", "dlPFC80655_p"),
    "GSE53987:striatum": ("STR_lfc", "STR_p"),
    "GSE53987:hippocampal formation": ("HPC_lfc", "HPC_p"),
    "GSE53987:dorsolateral prefrontal cortex": ("dlPFC53987_lfc", "dlPFC53987_p"),
}
with open(os.path.join(OUT, "partB_evidence_stack.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(stack_cols)
    for g in sorted(sgrg_syms):
        row = [g] + [None] * (len(stack_cols) - 1)
        d = dict(zip(stack_cols, row))
        d["gene"] = g
        u = uc.get(g)
        if u:
            d["UC_lfc_GSE38713"] = u["lfc"]; d["UC_FDR"] = u["fdr"]
        p3 = paired.get(g)
        if p3:
            d["UCpaired_dz_GSE22619"] = p3["dz"]; d["UCpaired_p"] = p3["p"]
        m1 = sgrg_in_19738.get(g)
        if m1:
            d["MDDblood_lfc_GSE19738"] = m1["lfc"]; d["MDDblood_p"] = m1["p"]; d["MDDblood_q"] = m1["q"]
        for label, (cl, cp) in label2col.items():
            v = brain_repl.get(g, {}).get(label)
            if v:
                d[cl] = v["lfc"]; d[cp] = v["p"]
        m4 = mr.get(g)
        if m4:
            d["MR_tissue"] = m4["tissue"]; d["MR_IVW_p"] = m4["IVW_p"]; d["MR_IVW_q"] = m4["IVW_q"]
        w.writerow(["" if d[c] is None else (repr(d[c]) if isinstance(d[c], float) else d[c]) for c in stack_cols])
log("written: partB_evidence_stack.tsv")

# ------------------------------------------------------------------ 保存
results = dict(
    universe_size=N, mdd_de_genes=K, sgrg_panel_size=len(sgrg_syms),
    sgrg_in_19738=len(panel_universe),
    enrichment=enrichment,
    concordance=dict(n=n_cmp, same_sign=same_sign, spearman_rho=float(rho),
                     spearman_p=float(rho_p), binom_p=float(binom_p),
                     nominal_both_n=len(nom), nominal_both_same=nom_same,
                     nominal_genes=[g for g, _, _ in nom]),
    triple=dict(strict=triple_info, nominal=triple_nom,
                uc_de_n=len(uc_de), mdd_de_panel_n=len(mdd_de_panel)),
    brain_replication_concordance=concord,
)
with open(os.path.join(OUT, "partB_bridge_results.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=2, default=float)
log("written: partB_bridge_results.json")

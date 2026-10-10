# -*- coding: utf-8 -*-
"""GSE133060 结果诊断：随机集零分布形态 + 检验校准。"""
import numpy as np
import pandas as pd

T = pd.read_csv("results/partD_gse133060_tests.tsv", sep="\t")
G = pd.read_csv("results/partD_gse133060_gene_scores.tsv", sep="\t")

x = G.log2FC_inflamed_vs_noninflamed.values
print("=== 基因级 Δz 分布（全部 %d 基因）===" % len(x))
print("  mean=%+.4f median=%+.4f sd=%.3f q05=%+.3f q95=%+.3f max=%+.2f min=%+.2f"
      % (x.mean(), np.median(x), x.std(), np.quantile(x, .05), np.quantile(x, .95),
         x.max(), x.min()))
print("  偏度=%.2f" % float(pd.Series(x).skew()))
print()

rd = T[T.kind == "random"].mean_paired_diff.values
sgrg = float(T[T.geneset == "SGRG"].mean_paired_diff.iloc[0])
print("=== 随机集 Δz 零分布 (n=%d) ===" % len(rd))
print("  mean=%+.4f median=%+.4f sd=%.3f" % (rd.mean(), np.median(rd), rd.std()))
for q in [1, 5, 25, 50, 75, 95, 99]:
    print("   q%02d = %+.4f" % (q, np.quantile(rd, q / 100)))
print("  SGRG Δz = %+.4f" % sgrg)
print()

n_up = int((rd >= sgrg).sum())
n_abs = int((np.abs(rd) >= abs(sgrg)).sum())
n_dn = int((rd <= -sgrg).sum())
print("  单侧(向上)经验 p = (%d+1)/%d = %.4f" % (n_up, len(rd), (n_up + 1) / (len(rd) + 1)))
print("  双侧(绝对值)经验 p = (%d+1)/%d = %.4f" % (n_abs, len(rd), (n_abs + 1) / (len(rd) + 1)))
print("  负向侧 n(<= -SGRG) = %d" % n_dn)
print("  SGRG 向上排名 = %d / %d" % (n_up + 1, len(rd)))
print()

print("=== 名义显著性（各集合）===")
for k in ["random", "positive", "target"]:
    s = T[T.kind == k]
    if not len(s):
        continue
    print("  %-9s n=%3d  名义p<0.05=%3d (%.1f%%)  n_up=8/8 的=%d"
          % (k, len(s), int((s.wilcoxon_p < 0.05).sum()), 100 * (s.wilcoxon_p < 0.05).mean(),
             int((s.n_up == 8).sum())))
print()

print("=== 随机集 Wilcoxon p 的取值分布（8 对 → 离散）===")
vc = T[T.kind == "random"].wilcoxon_p.round(4).value_counts().sort_index()
for p, c in vc.items():
    print("   p=%.4f  n=%d" % (p, c))

# -*- coding: utf-8 -*-
"""
Part B · Step 1 — 抓取 SGRG 基因在 MDD 三个队列（Gemma 存算差异分析）中的
MDD-对照统计量，外加 GSE38713 存算 UC 分析（批次混杂，仅作次要参考）。

数据通道（重要，可复用）：
  GET /rest/v2/datasets/analyses/differential/results/taxa/human/genes/{symbol}?query=<dataset>&limit=100
  - query 按数据集 shortName 全文检索，返回该数据集所有结果集中该基因的探针行
  - 每行 contrasts[] 含各疾病水平的 logFoldChange / pvalue（对照=基线）
  - 行级 pValue/correctedPvalue 为 Gemma 存的探针级统计（GSE19738 单对照时即 MDD 检验）

输出：results/partB_sgrg_mdd_de.tsv
"""
import os, sys, json, time, gzip, io, urllib.request, urllib.parse, urllib.error

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
OUT  = os.path.join(BASE, "results")
CACHE = os.path.join(BASE, "data", "partB_cache")
os.makedirs(CACHE, exist_ok=True)

SGRG_CLEAN = os.path.join(OUT, "sgrg_gene_set_clean.tsv")

# 数据集 → (查询名, 目标 resultSets, MDD 对照的 factorValue summary)
DATASETS = {
    "GSE19738": dict(label="MDD 外周血 (66 MDD / 66 healthy)",
                     rs_wanted=None,           # 全部结果集都要（540049 疾病 / 540050 治疗 / 540051 子集）
                     mdd_fv="major depressive disorder",
                     keep_fv=["major depressive disorder"]),
    "GSE80655": dict(label="MDD 伏隔核 RNA-seq (postmortem)",
                     rs_wanted=None,
                     mdd_fv="major depressive disorder",
                     keep_fv=["major depressive disorder"]),
    "GSE53987": dict(label="MDD 脑区 (PFC/HPC/STR)",
                     rs_wanted=None,
                     mdd_fv="major depressive disorder",
                     keep_fv=["major depressive disorder"]),
    "GSE38713": dict(label="UC 结肠活检 (存算分析, 批次混杂, 次要参考)",
                     rs_wanted=None,
                     mdd_fv="ulcerative colitis",
                     keep_fv=["ulcerative colitis", "inactive  ulcerative colitis",
                              "ulcerative colitis with non-involved"]),
}

def log(*a): print(*a, flush=True)

def load_sgrg():
    syms = []
    with open(SGRG_CLEAN, encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p and p[0]:
                syms.append(p[0])
    return syms

def http_get_json(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept-Encoding": "gzip",
                "User-Agent": "curl/8.0",
            })
            with urllib.request.urlopen(req, timeout=90) as r:
                raw = r.read()
                if r.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.GzipFile(fileobj=io.BytesIO(raw)).read()
                return json.loads(raw)
        except Exception as e:
            last = e
            time.sleep(2 + 3 * i)
    raise RuntimeError(f"GET failed after {tries} tries: {url} :: {last}")

def fetch_gene_dataset(gene, ds):
    """返回该基因在该数据集所有结果集的行"""
    cache = os.path.join(CACHE, f"{ds}__{gene}.json")
    if os.path.exists(cache):
        with open(cache, encoding="utf-8") as f:
            return json.load(f)
    url = ("https://gemma.msl.ubc.ca/rest/v2/datasets/analyses/differential/"
           f"results/taxa/human/genes/{urllib.parse.quote(gene)}"
           f"?query={urllib.parse.quote(ds)}&limit=100")
    d = http_get_json(url)
    rows = d.get("data", [])
    with open(cache, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False)
    return rows

def main():
    genes = load_sgrg()
    log(f"SGRG genes: {len(genes)}")
    records = []
    n_empty = 0
    for ds, meta in DATASETS.items():
        log(f"== {ds} — {meta['label']}")
        for gi, gene in enumerate(genes):
            try:
                rows = fetch_gene_dataset(gene, ds)
            except Exception as e:
                log(f"  !! {gene}: {e}")
                continue
            if not rows:
                n_empty += 1
            for r in rows:
                syms = [g.get("officialSymbol") for g in r.get("genes", [])]
                if gene not in syms and syms:
                    # 基因符号映射差异，跳过不相关行
                    continue
                for c in r.get("contrasts", []):
                    fv = c.get("factorValue") or {}
                    fvs = fv.get("summary") or fv.get("factorValue") or ""
                    if fvs not in meta["keep_fv"]:
                        continue
                    rec = dict(
                        dataset=ds,
                        gene=gene,
                        rs=r.get("resultSetId"),
                        exp=r.get("experimentAnalyzedId"),
                        probe=r.get("probeName"),
                        fv=fvs,
                        lfc=c.get("logFoldChange"),
                        p=c.get("pvalue"),
                        probe_p=r.get("pValue"),
                        probe_q=r.get("correctedPvalue"),
                    )
                    records.append(rec)
            if (gi + 1) % 20 == 0:
                log(f"  {gi+1}/{len(genes)} done (rows so far {len(records)})")
    # 输出
    path = os.path.join(OUT, "partB_sgrg_mdd_de.tsv")
    cols = ["dataset","gene","rs","exp","probe","fv","lfc","p","probe_p","probe_q"]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for rec in records:
            f.write("\t".join(
                "NA" if rec[c] is None else
                (rec[c] if isinstance(rec[c], str) else repr(rec[c]))
                for c in cols) + "\n")
    log(f"total records: {len(records)}; empty gene-dataset pairs: {n_empty}")
    log(f"written: {path}")

if __name__ == "__main__":
    main()

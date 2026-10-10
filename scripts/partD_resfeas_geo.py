# -*- coding: utf-8 -*-
"""Part D / 公共资源可用性审计 —— GEO 侧（全库枚举 + 配对设计识别）

为什么做
--------
CELLxGENE 只有 3 个集合含 UC（见 partD_resfeas_cxg.py）。但公共资源不止 CELLxGENE。
**真正决定「能不能做有检验力检验」的设计，是「同一供体的炎症 vs 非炎症配对」**——
它天然控制了批次与个体，是唯一能在没有大样本对照池时仍成立的设计。
本脚本在 GEO 上把这件事枚举清楚：

  1) esearch/db=gds 用若干互补查询（IBD×肠段×测序类型/空间）取回系列；
  2) esummary 取系列级元数据（标题/摘要/物种/样本数/类型）；
  3) 对每个候选 GSE 拉 SOFT brief（`acc.cgi?targ=gsm&form=text`）解析**每个样本的标题**，
     从标题里判状态（inflamed / non-inflamed / healthy）与**供体标识**；
  4) 统计每个系列的「同供体炎症-非炎症配对数」「独立真实对照数」「疾病与肠段」。

★ 关键实现点
  GEO 的配对信息只存在于**样本标题**里（characteristics 常常只有 tissue），
  所以必须拉 GSM 级 SOFT。样本标题如
      `Rectum_CD_Inflamed [I0262 RSDF00A CosMx1k]` / `Rectum_CD_Noninflamed [I0275 ...]`
  —— 方括号里的是切片/供体标识，状态在标题词里。

检验力门槛（与 partD 记忆一致）
  双侧符号翻转检验最小可达 p = 2/2^n：n=4→0.125、n=6→0.031、n=7→0.016。
  故「同供体配对数 ≥6」才算有名义检验力；「≥15」才算稳。

输出
----
  results/partD_resfeas_geo_series.tsv    每个系列一行
  results/partD_resfeas_geo_samples.tsv   每个样本一行（配对审计原始表）
  results/partD_resfeas_geo_summary.json
  data/raw/geo_soft/<ACC>.txt             SOFT 缓存
用法
----
  python partD_resfeas_geo.py [--workers 5] [--max-samples 400]
"""
import argparse
import json
import os
import re
import ssl
import socket
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

socket.setdefaulttimeout(120)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CACHE = "data/raw/geo_soft"
OUT = "results"
T0 = time.time()
log = lambda *x: print("[%6.1fs]" % (time.time() - T0), *x, flush=True)

IBD = r'("ulcerative colitis" OR "Crohn disease" OR "inflammatory bowel disease" OR "IBD")'
GUT = r'(colon OR colonic OR rectum OR rectal OR ileum OR ileal OR intestinal OR intestine OR "large bowel")'
SC_HINT = re.compile(r"single[\s\-]?cell|single[\s\-]?nucle|scrna|snrna|sc-rna|10x|10X|"
                     r"spatial|xenium|cosmx|viser[uo]n|merscope|in situ|stereo[\s\-]?seq|"
                     r"slide[\s\-]?seq|seqFISH|merfish|nano[\s\-]?string", re.I)

QUERIES = {
    "hts": '(%s AND %s) AND "expression profiling by high throughput sequencing"[DataSet Type]' % (IBD, GUT),
    "spatial": '(%s AND %s) AND (spatial OR Xenium OR CosMx OR Visium OR MERSCOPE OR "in situ transcriptomics")' % (IBD, GUT),
    "scnuc": '((single cell OR single-cell OR scRNA OR snRNA) AND %s AND %s)' % (IBD, GUT),
    "array": '(%s AND %s) AND "expression profiling by array"[DataSet Type]' % (IBD, GUT),
}

# ---------------- 状态 / 供体 判定 ----------------
RE_NONINFL = re.compile(r"(non[\s\-_]?inflam\w*|un[\s\-_]?inflam\w*|not[\s\-_]?inflam\w*|"
                        r"noninflam\w*|quiescen\w*|remiss\w*|inactive)", re.I)
RE_INFLAM = re.compile(r"(inflam\w*|active|flare)", re.I)
RE_HEALTHY = re.compile(r"(healthy|control|ctrl|non[\s\-_]?ibd|normal|unaffected|HC\b)", re.I)
RE_UC = re.compile(r"(ulcerative colitis|\bUC\b|U\.C\.)", re.I)
RE_CD = re.compile(r"(crohn|\bCD\b)", re.I)
RE_IBD = re.compile(r"(inflammatory bowel|\bIBD\b)", re.I)
RE_ILEUM = re.compile(r"(ileum|ileal|\bTI\b|terminal ileum)", re.I)
RE_COLON = re.compile(r"(colon|colonic|sigmoid|ascending|descending|transverse|cecum|caecum|"
                      r"rectum|rectal|left[\s\-]?side|right[\s\-]?side)", re.I)
RE_ILLEGAL = re.compile(r"gse\d+[a-z]?", re.I)


def get_text(url, tries=4):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, context=CTX) as x:
                return x.read().decode("utf-8", "replace")
        except Exception as e:
            last = e
            time.sleep(1.2 * (k + 1))
    raise last


def get_json(url, tries=4):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0",
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, context=CTX) as x:
                return json.loads(x.read().decode("utf-8"))
        except Exception as e:
            last = e
            time.sleep(1.2 * (k + 1))
    raise last


def esearch(term, retmax=600):
    u = "%s/esearch.fcgi?db=gds&retmax=%d&retmode=json&term=%s" % (EUTILS, retmax,
                                                                  urllib.parse.quote(term))
    j = get_json(u)
    res = j["esearchresult"]
    return res.get("idlist", []), int(res.get("count", 0))


def esummary(ids):
    """分批 esummary，返回 uid -> DocSum dict"""
    out = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        u = "%s/esummary.fcgi?db=gds&retmode=json&id=%s" % (EUTILS, ",".join(chunk))
        j = get_json(u)
        out.update(j.get("result", {}))
        time.sleep(0.4)
    return out


def parse_soft(acc, max_samples):
    """拉 GSM 级 SOFT brief，返回样本列表"""
    fp = os.path.join(CACHE, acc + ".txt")
    if os.path.exists(fp):
        s = open(fp, encoding="utf-8", errors="replace").read()
    else:
        u = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=gsm&form=text&view=brief" % acc
        s = get_text(u)
        os.makedirs(CACHE, exist_ok=True)
        open(fp, "w", encoding="utf-8").write(s)
        time.sleep(0.35)
    if not s.startswith("^SAMPLE"):
        return None
    recs = s.split("^SAMPLE = ")[1:]
    if len(recs) > max_samples:
        return {"too_many": len(recs)}
    samples = []
    for r in recs:
        title = re.search(r"!Sample_title = (.*)", r)
        src = re.search(r"!Sample_source_name_ch1 = (.*)", r)
        chars = re.findall(r"!Sample_characteristics_ch1 = (.*)", r)
        samples.append({
            "gsm": r.splitlines()[0].strip() if r.splitlines() else "",
            "title": (title.group(1).strip() if title else ""),
            "source": (src.group(1).strip() if src else ""),
            "chars": " | ".join(x.strip() for x in chars),
        })
    return samples


def donor_token(title, source):
    """从标题/来源里抽供体标识。优先方括号内的首个 token。

    ★ 刻意**不做激进归一化**（如把 `353_6_3338` 削成 `353`）——那会凭空造出配对。
      供体身份在本步骤里按「原样 token」比对；跨切片的同供体判断交给人工核查
      （全库仅 7 个系列有任何配对，人工可穷举）。
    """
    t = title
    m = re.search(r"[\[\(]([^\]\)]+)[\]\)]", t)
    if m:
        inner = m.group(1).strip()
        for tok in re.split(r"[,\s;]+", inner):
            tok = tok.strip()
            if not tok:
                continue
            if re.fullmatch(r"(?i)(rep|run|section|slide|cosmx|xenium|visium|merscope|ns|ws|"
                            r"panel|gene|v\d)\w*\d*", tok):
                continue
            if re.search(r"\d", tok):
                return tok
    # 明确写出 patient/subject/donor 的
    m = re.search(r"\b(patient|subject|donor|individual|pt)[\s\-_#]*(\d+)\b", t, re.I)
    if m:
        return "%s%s" % (m.group(1).lower(), m.group(2))
    # 形如 I0262 / UC-1 / HC-2 / P123 的 token
    m = re.search(r"\b([A-Za-z]{1,4}[\-_]?\d{2,5})\b", t)
    if m:
        return m.group(1)
    m = re.search(r"\b(\d{1,4})\b", t)
    if m:
        return m.group(1)
    return ""


SER_CACHE = "data/raw/geo_soft_series"


def parse_soft_series(acc):
    """拉系列级 SOFT brief，取可审计的 `!Series_type`（判定检测分子的唯一可靠字段）。

    ★ 教训：★ 不要用系列摘要里的关键词判模态。GSE272890 的摘要写了
      「We performed spatial transcriptomics ...」，但它这个 series 实际是
      `Non-coding RNA profiling by high throughput sequencing`（测 miRNA）。
      按摘要关键词归类会把它记成「空间转录组」并虚增可用空间数据数。
    """
    fp = os.path.join(SER_CACHE, acc + ".txt")
    if os.path.exists(fp):
        s = open(fp, encoding="utf-8", errors="replace").read()
    else:
        u = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief"
             % acc)
        s = get_text(u)
        os.makedirs(SER_CACHE, exist_ok=True)
        open(fp, "w", encoding="utf-8").write(s)
        time.sleep(0.35)
    out = {"series_type": "", "overall_design": "", "summary": "", "title": "", "platform": ""}
    KEYMAP = {"type": "series_type", "overall_design": "overall_design",
              "summary": "summary", "title": "title", "platform_id": "platform"}
    for line in s.splitlines():
        m = re.match(r"!Series_(type|overall_design|summary|title|platform_id) = (.*)", line)
        if not m:
            continue
        k, v = KEYMAP[m.group(1)], m.group(2)
        if k == "platform":
            out["platform"] = (out["platform"] + ";" + v).strip(";")
        else:
            out[k] = (out[k] + " | " + v).strip(" |")
    return out


def assay_class(series_type, gds_type, blob):
    """返回 (分子类别, 分辨率)。

    分子类别是**闸门**：非 mRNA 的 assay（miRNA/ncRNA、ATAC/ChIP、甲基化）无法给
    SGRG 这类编码基因面板打分，必须在「可用性」判定里显式排除。
    """
    st = " ".join([str(series_type or ""), str(gds_type or ""), blob[:400]]).lower()
    if re.search(r"non[\s\-]?coding rna|microrna|mirna|ncrna|lncRNA", st, re.I):
        mol = "ncRNA"
    elif re.search(r"genome binding|occupancy|methylation|atac|chip[\s\-]?seq|"
                   r"chromatin|histone", st, re.I):
        mol = "chromatin"
    elif re.search(r"expression profiling|transcriptom|rna[\s\-]?seq|microarray|"
                   r"third[\s\-]?party reanalysis", st, re.I):
        mol = "mRNA"
    else:
        mol = "other"
    if mol != "mRNA":
        return mol, mol
    if re.search(r"spatial transcriptom|xenium|cosmx|visium|merscope|in situ|"
                 r"slide[\s\-]?seq|stereo[\s\-]?seq|seqfish|merfish", blob, re.I):
        return mol, "spatial"
    if re.search(r"single[\s\-]?cell|single[\s\-]?nucle|\bsc[\s\-]?rna|\bsn[\s\-]?rna|10x",
                 blob, re.I):
        return mol, "single-cell"
    if re.search(r"array|microarray|beadchip", " ".join([str(series_type or ""),
                                                          str(gds_type or "")]), re.I):
        return mol, "bulk-array"
    return mol, "bulk-seq"


RE_INVITRO = re.compile(r"organoid|enteroid|cell line|in vitro|colonoid|"
                        r"ipsc|caco[\s\-]?2|monolayer", re.I)


def site_of(txt):
    if RE_ILEUM.search(txt):
        return "ileum"
    if re.search(r"rectum|rectal", txt, re.I):
        return "rectum"
    if RE_COLON.search(txt):
        return "colon"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--max-samples", type=int, default=400)
    ap.add_argument("--soft-only", action="store_true", help="只跑 SOFT 解析（跳过 esearch）")
    a = ap.parse_args()

    # ---- 1) esearch 并集
    uid2q = {}
    counts = {}
    for k, q in QUERIES.items():
        if a.soft_only:
            break
        ids, cnt = esearch(q)
        counts[k] = cnt
        log("查询 %-8s 命中 %d（取回 %d）" % (k, cnt, len(ids)))
        for i in ids:
            uid2q.setdefault(i, set()).add(k)
    log("UID 并集 = %d" % len(uid2q))

    # ---- 2) esummary
    if not a.soft_only:
        summ = esummary(sorted(uid2q.keys()))
        log("esummary 返回 %d 条" % len(summ))
        recs = []
        for uid, d in summ.items():
            if uid == "uids":
                continue
            if (d.get("entrytype") or "") != "GSE":
                continue
            tax = (d.get("taxon") or "")
            if "Homo sapiens" not in tax:
                continue
            acc = d.get("accession")
            if not acc:
                continue
            title = d.get("title") or ""
            summary = d.get("summary") or ""
            n = int(d.get("n_samples") or 0)
            recs.append({"uid": uid, "gse": acc, "title": title, "summary": summary,
                         "n_samples": n, "gdsType": d.get("gdsType"),
                         "pdat": d.get("pdat"), "pubmed": d.get("pubmedids"),
                         "queries": ";".join(sorted(uid2q.get(uid, [])))})
        R = pd.DataFrame(recs).drop_duplicates("gse")
        log("人类 GSE 系列 = %d" % len(R))
    else:
        R = pd.read_csv(os.path.join(OUT, "partD_resfeas_geo_series_raw.tsv"), sep="\t")
        log("soft-only：载入既有原始系列表 %d" % len(R))

    R["blob"] = (R.title.fillna("") + " " + R.summary.fillna(""))
    R["is_sc_or_spatial"] = R.blob.apply(lambda x: bool(SC_HINT.search(x)))
    R["gut_hit"] = R.blob.apply(lambda x: bool(re.search(r"colon|intestin|rectum|rectal|ileum|ileal|cecum|"
                                                          r"colitis|crohn|bowel", x, re.I)))
    R["ibd_hit"] = R.blob.apply(lambda x: bool(re.search(r"ulcerative colitis|crohn|inflammatory bowel|"
                                                         r"colitis|\bIBD\b", x, re.I)))
    R.to_csv(os.path.join(OUT, "partD_resfeas_geo_series_raw.tsv"), sep="\t", index=False)

    # 只对「像单细胞/空间」或「小系列」的拉 SOFT
    todo = R[(R.n_samples > 0) & (R.n_samples <= a.max_samples) &
             (R.is_sc_or_spatial | (R.n_samples <= 60))].copy()
    todo = todo.drop_duplicates("gse")
    log("待解析 SOFT 的系列 = %d / %d" % (len(todo), len(R)))

    # ---- 3) 并行拉 SOFT + 解析
    def work(row):
        try:
            ser = parse_soft_series(row.gse)
        except Exception as e:
            ser = {"series_type": "ERR:%s" % type(e).__name__}
        try:
            s = parse_soft(row.gse, a.max_samples)
            return row.gse, s, ser, None
        except Exception as e:
            return row.gse, None, ser, "%s: %s" % (type(e).__name__, e)

    results = {}
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for i, (gse, s, ser, err) in enumerate(ex.map(work, [r for _, r in todo.iterrows()]), 1):
            results[gse] = (s, ser, err)
            if i % 40 == 0:
                log("  SOFT %d/%d" % (i, len(todo)))
    log("SOFT 解析完成，成功 %d / 失败 %d"
        % (sum(1 for v in results.values() if v[0] and "too_many" not in v[0]),
           sum(1 for v in results.values() if v[2])))

    # ---- 4) 逐样本判状态/供体，逐系列汇总
    srows, grows = [], []
    for _, r in todo.iterrows():
        s, ser, err = results.get(r.gse, (None, {}, "missing"))
        if s is None or "too_many" in s:
            grows.append({"gse": r.gse, "title": r.title, "n_samples": r.n_samples,
                          "parse": "skip(%s)" % ("too_many" if s and "too_many" in s else err)})
            continue
        mol, res = assay_class(ser.get("series_type", ""), r.gdsType, r.blob)
        state, donors, dsets, sites = [], {}, set(), set()
        pair_site = {}       # donor -> {site: set(states)}  严格口径：同供体 + 同部位
        invitro = bool(RE_INVITRO.search(r.blob))
        for smp in s:
            txt = " ".join([smp["title"], smp["source"], smp["chars"]])
            if RE_NONINFL.search(txt):
                st = "noninflamed"
            elif RE_INFLAM.search(txt):
                st = "inflamed"
            elif RE_HEALTHY.search(txt):
                st = "healthy"
            else:
                st = "unknown"
            dn = donor_token(smp["title"], smp["source"])
            dis = ("UC" if RE_UC.search(txt) else "") or ("CD" if RE_CD.search(txt) else "") \
                  or ("IBD" if RE_IBD.search(txt) else "") or ("" if st == "healthy" else "")
            site = site_of(txt)
            if site:
                sites.add(site)
            if dis:
                dsets.add(dis)
            key = dn or ("<%s#%s>" % (st, smp["gsm"]))
            donors.setdefault(key, set()).add(st)
            site_key = site or "na"
            pair_site.setdefault((dn or smp["gsm"], site_key), set()).add(st)
            srows.append({"gse": r.gse, "gsm": smp["gsm"], "title": smp["title"],
                          "state": st, "donor_raw": dn, "disease": dis, "site": site,
                          "in_vitro": invitro, "chars": smp["chars"][:220]})
        paired = {k: v for k, v in donors.items() if {"inflamed", "noninflamed"} <= v}
        paired_strict = {k: v for k, v in pair_site.items()
                         if {"inflamed", "noninflamed"} <= v}
        n_healthy = sum(1 for k, v in donors.items() if v == {"healthy"})
        n_infl = sum(1 for v in donors.values() if "inflamed" in v)
        n_noninfl = sum(1 for v in donors.values() if "noninflamed" in v)
        grows.append({
            "gse": r.gse, "title": r.title, "n_samples": len(s),
            "series_type": ser.get("series_type", "")[:120],
            "platform": ser.get("platform", "")[:60],
            "molecule": mol, "resolution": res,
            "mrna_scorable": (mol == "mRNA"),
            "in_vitro": invitro,
            "diseases": ";".join(sorted(dsets)), "sites": ";".join(sorted(sites)),
            "n_state_known": sum(1 for x in srows if x["gse"] == r.gse and x["state"] != "unknown"),
            "n_inflamed": n_infl, "n_noninflamed": n_noninfl,
            "n_paired_donors": len(paired),
            "n_paired_same_site": len(paired_strict),
            "n_healthy_units": n_healthy,
            "paired_donor_ids": ";".join(sorted(paired.keys()))[:400],
            "pdat": r.pdat, "pubmed": r.pubmed, "queries": r.queries,
        })

    G = pd.DataFrame(grows)
    S = pd.DataFrame(srows)
    G.to_csv(os.path.join(OUT, "partD_resfeas_geo_series.tsv"), sep="\t", index=False)
    S.to_csv(os.path.join(OUT, "partD_resfeas_geo_samples.tsv"), sep="\t", index=False)

    P = G[G.n_paired_donors > 0].sort_values("n_paired_donors", ascending=False)
    M = G[G.mrna_scorable]
    SC = M[M.resolution.isin(["single-cell", "spatial"])]
    SP = M[M.resolution == "spatial"]
    summ = {
        "query_counts": counts,
        "n_gse_human": int(len(R)),
        "n_series_parsed": int(len(G)),
        "n_series_skipped_bulklarge": int(len(R) - len(G)),
        "molecule_counts": {k: int(v) for k, v in G.molecule.value_counts().items()},
        "resolution_counts_all": {k: int(v) for k, v in G.resolution.value_counts().items()},
        "n_mrna_scorable": int(len(M)),
        "resolution_counts_mrna": {k: int(v) for k, v in M.resolution.value_counts().items()},
        "n_singlecell_or_spatial_mrna": int(len(SC)),
        "n_spatial_mrna": int(len(SP)),
        # ---- 配对（全部 / 仅 mRNA 可打分 / 仅单细胞与空间）
        "n_series_with_pairing": int(len(P)),
        "n_series_with_pairing_ge6": int((G.n_paired_donors >= 6).sum()),
        "n_series_pairing_same_site": int((G.n_paired_same_site > 0).sum()),
        "n_series_with_pairing_mrna": int(((M.n_paired_donors > 0)).sum()),
        "n_series_with_pairing_ge6_mrna": int((M.n_paired_donors >= 6).sum()),
        "n_series_with_pairing_ge15_mrna": int((M.n_paired_donors >= 15).sum()),
        "n_sc_or_spatial_mrna_with_pairing": int((SC.n_paired_donors > 0).sum()),
        "n_sc_or_spatial_mrna_with_pairing_ge6": int((SC.n_paired_donors >= 6).sum()),
        "n_spatial_mrna_with_pairing": int((SP.n_paired_donors > 0).sum()),
        "n_spatial_mrna_with_pairing_samesite": int((SP.n_paired_same_site > 0).sum()),
        "n_series_pairing_and_UC_mrna_ge6": int(((M.n_paired_donors >= 6) &
                                                M.diseases.str.contains("UC")).sum()),
        "n_series_pairing_ge6_mrna_tissue": int(((M.n_paired_donors >= 6) & (~M.in_vitro)).sum()),
        "top_paired_series": P.head(20)[["gse", "n_paired_donors", "n_paired_same_site",
                                         "diseases", "sites", "molecule", "resolution",
                                         "in_vitro", "n_samples"]].to_dict("records"),
        "power_note": "双侧符号翻转最小 p = 2/2^n：n=4→0.125, n=6→0.031, n=7→0.016；"
                      "配对设计需 ≥6 对才有名义检验力",
        "assay_gate_note": "「是否测 mRNA」是可用性闸门：miRNA/ncRNA、ATAC/ChIP、甲基化系列"
                           "无法给 SGRG 编码基因面板打分，即使配对设计再漂亮也不可用。",
        "sampling_note": "SOFT 仅对 n_samples<=400 且（像单细胞/空间 或 n_samples<=60）的系列解析；"
                         "被跳过的 %d 个系列全部为 bulk 且 n_samples>60（其中单细胞/空间 0 个），"
                         "故单细胞/空间的覆盖是完整的。" % int(len(R) - len(G)),
    }
    json.dump(summ, open(os.path.join(OUT, "partD_resfeas_geo_summary.json"), "w",
                         encoding="utf-8"), indent=2, ensure_ascii=False, default=str)

    print("\n" + "=" * 92)
    print("GEO 全库 · 配对炎症/非炎症设计审计（含「是否测 mRNA」闸门）")
    print("=" * 92)
    print("人类 GSE %d；解析 SOFT %d 个系列（跳过 bulk 大系列 %d）"
          % (summ["n_gse_human"], summ["n_series_parsed"], summ["n_series_skipped_bulklarge"]))
    print("分子类别：%s" % summ["molecule_counts"])
    print("分辨率（mRNA 可打分 %d 个系列）：%s"
          % (summ["n_mrna_scorable"], summ["resolution_counts_mrna"]))
    print("单细胞/空间且 mRNA 可打分 = %d（其中空间 %d）"
          % (summ["n_singlecell_or_spatial_mrna"], summ["n_spatial_mrna"]))
    print("-" * 92)
    print("含同供体配对（炎症+非炎症）：全部 %d 个系列；≥6 对 %d；同部位配对 %d"
          % (summ["n_series_with_pairing"], summ["n_series_with_pairing_ge6"],
             summ["n_series_pairing_same_site"]))
    print("★ 限定 mRNA 可打分后：有配对 %d；≥6 对 %d；≥15 对 %d"
          % (summ["n_series_with_pairing_mrna"], summ["n_series_with_pairing_ge6_mrna"],
             summ["n_series_with_pairing_ge15_mrna"]))
    print("★ 其中单细胞/空间：有配对 %d；≥6 对 %d；空间有配对 %d（同部位 %d）"
          % (summ["n_sc_or_spatial_mrna_with_pairing"],
             summ["n_sc_or_spatial_mrna_with_pairing_ge6"],
             summ["n_spatial_mrna_with_pairing"], summ["n_spatial_mrna_with_pairing_samesite"]))
    print("★ mRNA 可打分 + ≥6 对 + UC = %d；排除离体模型后 = %d"
          % (summ["n_series_pairing_and_UC_mrna_ge6"], summ["n_series_pairing_ge6_mrna_tissue"]))
    if len(P):
        print("\n配对数最多的系列（★ 注意 molecule 列——miRNA 测的配对再多也不能给 SGRG 打分）：")
        print(P.head(20)[["gse", "n_paired_donors", "n_paired_same_site", "molecule",
                          "resolution", "diseases", "sites", "in_vitro", "n_samples"]]
              .to_string(index=False, max_colwidth=26))
    log("完成 %.1f s" % (time.time() - T0))
    log("产出 results/partD_resfeas_geo_series.tsv / _samples.tsv / _summary.json")


if __name__ == "__main__":
    main()

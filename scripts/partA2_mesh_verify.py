# -*- coding: utf-8 -*-
"""A2 步 0：对 160 个「关键词命中」的 UC 队列做独立疾病归因核验 + 供体级去重

为什么要做
----------
A2 的队列池来自 `query=ulcerative colitis`（Gemma 的关键词/文本相似检索）。
这类归因本身有量级的失效率（已发表审计：36.4% 关键词式疾病归因被出版物 MeSH 推翻）。
而「160 → 15」这个数字一旦进入引言就是承重结构，必须用**独立于 Gemma 的策展**去核验。

核验路线（全部走 NCBI，独立于 Gemma 的 curation）
-------------------------------------------------
  GSE accession --esearch(db=gds)--> GDS uid --esummary(db=gds)--> PubMedID(s)
      --efetch(db=pubmed, XML)--> MeSH DescriptorName + MajorTopicYN

判据
----
  uc_major : 至少一篇关联文献的 MeSH 中 "Colitis, Ulcerative" 标为 MajorTopicYN=Y
  uc_any   : 至少一篇含 "Colitis, Ulcerative"（任意主/次）
  ibd_only : 含 IBD 相关词（Crohn / Inflammatory Bowel / Colitis）但不含 UC

供体级去重
----------
  · 同一 PubMedID 的多个 GSE = 同一研究 → 供体只应计一次
  · 同一数据集不重复（1 GSE = 1 Gemma 数据集）
  另：GEO 的 SuperSeries/SubSeries 会以不同 GSE 出现但指向同一篇文献 → 由 PMID 捕获

用法
----
  python partA2_mesh_verify.py [--refresh]
"""
import argparse
import collections
import json
import os
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CACHE = os.path.join(ROOT, "data", "raw", "mesh")
OUT = os.path.join(ROOT, "results")
os.makedirs(CACHE, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
HDR = {"User-Agent": "Mozilla/5.0 (academic research; a2mesh)", "Accept": "*/*"}
TOOL = "a2_mesh_verify"
MAIL = "a2mesh@example.org"
PAUSE = 0.36           # 无 API key 时限 3 req/s
AP = argparse.ArgumentParser()
AP.add_argument("--refresh", action="store_true")
ARGS = AP.parse_args()

_t0 = time.time()
def log(*x):
    print("[%6.1fs]" % (time.time() - _t0), *x, flush=True)


def http(u, retries=4):
    for k in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=HDR),
                                        context=CTX, timeout=60) as r:
                return r.read()
        except Exception as e:
            if k == retries - 1:
                raise
            time.sleep(1.5 * (k + 1))


def cache_path(name):
    return os.path.join(CACHE, name)


def fetch_json(name, url, refresh=False):
    p = cache_path(name)
    if os.path.exists(p) and os.path.getsize(p) > 2 and not refresh:
        try:
            return json.load(open(p, encoding="utf-8"))
        except Exception:
            pass
    raw = http(url)
    d = json.loads(raw)
    json.dump(d, open(p, "w", encoding="utf-8"))
    return d


def fetch_text(name, url, refresh=False):
    p = cache_path(name)
    if os.path.exists(p) and os.path.getsize(p) > 2 and not refresh:
        return open(p, encoding="utf-8").read()
    raw = http(url).decode("utf-8", "replace")
    open(p, "w", encoding="utf-8").write(raw)
    return raw


def eurl(endpoint, **kw):
    kw.setdefault("tool", TOOL)
    kw.setdefault("email", MAIL)
    return EUTILS + endpoint + "?" + urllib.parse.urlencode(kw)


# ---------------------------------------------------------------- 1. accession -> PubMed IDs
def acc_to_pmids(acc):
    """GSE accession -> (gds_uid, [pmids], taxon)"""
    refresh = ARGS.refresh
    tag = "es_gds_%s" % acc
    u = eurl("esearch.fcgi", db="gds", term="%s[Accession] AND gse[ETYP]" % acc, retmode="json")
    d = fetch_json(tag + ".json", u, refresh)
    time.sleep(PAUSE)
    uids = (d.get("esearchresult") or {}).get("idlist") or []
    if not uids:                      # 回退：不加 ETYP 过滤
        u = eurl("esearch.fcgi", db="gds", term="%s[Accession]" % acc, retmode="json")
        d = fetch_json(tag + "_all.json", u, refresh)
        time.sleep(PAUSE)
        uids = (d.get("esearchresult") or {}).get("idlist") or []
        # 只保留形如 200xxxxxx 的 GSE 记录
        uids = [x for x in uids if str(x).startswith("200")]
    if not uids:
        return None, [], None
    uid = uids[0]
    su = eurl("esummary.fcgi", db="gds", id=uid, retmode="json")
    s = fetch_json("es_sum_%s.json" % uid, su, refresh)
    time.sleep(PAUSE)
    r = (s.get("result") or {}).get(str(uid)) or {}
    pmids = [str(x) for x in (r.get("pubmedids") or [])]
    return uid, pmids, r.get("taxon")


# ---------------------------------------------------------------- 2. PubMed -> MeSH
MH_RE = re.compile(r'<DescriptorName UI="[^"]*" MajorTopicYN="([YN])">([^<]+)</DescriptorName>')
PT_RE = re.compile(r"<PublicationType[^>]*>([^<]+)</PublicationType>")


def pmid_mesh(pmid):
    xml = fetch_text("pm_%s.xml" % pmid, eurl("efetch.fcgi", db="pubmed", id=pmid, retmode="xml"))
    time.sleep(PAUSE)
    mesh = [(m.group(2), m.group(1)) for m in MH_RE.finditer(xml)]
    pts = [x for x in PT_RE.findall(xml)]
    title = ""
    m = re.search(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", xml, re.S)
    if m:
        title = re.sub(r"<[^>]+>", "", m.group(1)).strip()
    year = ""
    m = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", xml, re.S)
    if m:
        year = m.group(1)
    return {"mesh": mesh, "pub_types": pts, "title": title[:300], "year": year}


UC = "colitis, ulcerative"
IBD_WORDS = ("colitis", "crohn", "inflammatory bowel")


def classify(mesh_lists):
    all_terms = [(t, f) for ml in mesh_lists for (t, f) in ml]
    low = [(t.lower(), f) for t, f in all_terms]
    uc_major = any(t == UC and f == "Y" for t, f in low)
    uc_any = any(t == UC for t, f in low)
    ibd = any(any(w in t for w in IBD_WORDS) for t, _ in low)
    if uc_major:
        return "UC确认-主要主题", uc_major, uc_any, ibd
    if uc_any:
        return "UC确认-次要", uc_major, uc_any, ibd
    if ibd:
        return "IBD但非UC", uc_major, uc_any, ibd
    return "无UC/IBD归因", uc_major, uc_any, ibd


def main():
    idx = json.load(open(os.path.join(ROOT, "data", "raw", "gemma_a2", "cohort_index.json"),
                         encoding="utf-8"))
    import pandas as pd
    coh = pd.read_csv(os.path.join(OUT, "partA2_cohorts.tsv"), sep="\t")
    arms = {int(r.gid): (r.n_uc, r.n_ctrl) for r in coh.itertuples()}

    log("队列数 =", len(idx))
    rows = []
    pmids_all = set()
    for i, (gid, e) in enumerate(sorted(idx.items(), key=lambda kv: int(kv[0]))):
        acc = str(e.get("accession") or "")
        uid, pmids, taxon = (None, [], None)
        if acc.startswith("GSE"):
            try:
                uid, pmids, taxon = acc_to_pmids(acc)
            except Exception as ex:
                log("  !! %s %s %s" % (acc, type(ex).__name__, ex))
        pmids_all.update(pmids)
        rows.append({"gid": int(gid), "acc": acc, "gds_uid": uid, "pmids": ";".join(pmids),
                     "n_pmid": len(pmids), "taxon_gds": taxon})
        if (i + 1) % 20 == 0:
            log("  [1] %d/%d  (累计 PMID %d)" % (i + 1, len(idx), len(pmids_all)))

    log("唯一 PMID = %d，开始取 MeSH（每条 %.2fs）" % (len(pmids_all), PAUSE))
    mesh_of = {}
    for j, pmid in enumerate(sorted(pmids_all)):
        try:
            mesh_of[pmid] = pmid_mesh(pmid)
        except Exception as ex:
            log("  !! PMID %s %s %s" % (pmid, type(ex).__name__, ex))
            mesh_of[pmid] = {"mesh": [], "pub_types": [], "title": "", "year": ""}
        if (j + 1) % 25 == 0:
            log("  [2] %d/%d" % (j + 1, len(pmids_all)))

    # ------------------------------------------------------------ 归类
    for r in rows:
        pms = [x for x in r["pmids"].split(";") if x]
        ml = [mesh_of.get(p, {}).get("mesh", []) for p in pms]
        verdict, uc_major, uc_any, ibd = classify(ml) if pms else ("无法核验-无PMID", False, False, False)
        r.update({"verdict": verdict, "uc_major": uc_major, "uc_any": uc_any, "ibd_any": ibd})
        pts = [x for p in pms for x in mesh_of.get(p, {}).get("pub_types", [])]
        r["is_review"] = any("review" in str(x).lower() for x in pts)
        r["years"] = ";".join(sorted({mesh_of.get(p, {}).get("year", "") for p in pms if mesh_of.get(p, {}).get("year")}))
        r["mesh_uc_major_n"] = sum(1 for t, f in [x for p in pms for x in mesh_of.get(p, {}).get("mesh", [])]
                                   if t.lower() == UC and f == "Y")
        nu, nc = arms.get(r["gid"], (None, None))
        r["n_uc"], r["n_ctrl"] = nu, nc
        r["arms_ok"] = bool(nu is not None and nc is not None and nu >= 4 and nc >= 4)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "partA2_mesh_attribution.tsv"), sep="\t",
              index=False, float_format="%.5g")

    # ------------------------------------------------------------ 汇总
    n = len(df)
    vc = df["verdict"].value_counts().to_dict()
    pub = collections.defaultdict(list)
    for r in rows:
        for p in [x for x in r["pmids"].split(";") if x]:
            pub[p].append(r["gid"])
    multi = {p: g for p, g in pub.items() if len(g) > 1}

    summary = {
        "n_keyword_hits": n,
        "verdict_counts": vc,
        "n_with_pmid": int((df.n_pmid > 0).sum()),
        "n_without_pmid": int((df.n_pmid == 0).sum()),
        "n_uc_any": int(df.uc_any.sum()),
        "n_uc_major": int(df.uc_major.sum()),
        "n_ibd_only": int((df.verdict == "IBD但非UC").sum()),
        "n_no_attribution": int((df.verdict == "无UC/IBD归因").sum()),
        "n_unique_pmids": len(pub),
        "n_pmids_multi_dataset": len(multi),
        "gids_in_multi_pmid_datasets": int(sum(len(g) for g in multi.values())),
        "n_review_datasets": int(df.is_review.sum()),
        "arms_ok_keyword": int(df.arms_ok.sum()),
        "arms_ok_and_uc_any": int((df.arms_ok & df.uc_any).sum()),
        "arms_ok_and_uc_major": int((df.arms_ok & df.uc_major).sum()),
    }
    # 供体级口径：UC 归因成立 + 两臂够 的数据集，按 PMID 去重后的研究数
    sel = df[df.arms_ok & df.uc_any]
    pubs_sel = {p for r in sel.to_dict("records") for p in r["pmids"].split(";") if p}
    summary["study_level_arms_ok_uc_any"] = len(pubs_sel) if pubs_sel else len(sel)
    sel2 = df[df.arms_ok & df.uc_major]
    pubs_sel2 = {p for r in sel2.to_dict("records") for p in r["pmids"].split(";") if p}
    summary["study_level_arms_ok_uc_major"] = len(pubs_sel2) if pubs_sel2 else len(sel2)

    json.dump(summary, open(os.path.join(OUT, "partA2_mesh_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1, default=str)

    print("\n" + "=" * 74)
    print("A2 步0：关键词归因的独立核验结果")
    print("=" * 74)
    for k, v in summary.items():
        print("  %-34s %s" % (k, v))
    print("\n  归因类别分布：")
    for k, v in sorted(vc.items(), key=lambda kv: -kv[1]):
        print("    %-18s %d (%.1f%%)" % (k, v, 100 * v / n))
    log("完成")
    return df


if __name__ == "__main__":
    main()

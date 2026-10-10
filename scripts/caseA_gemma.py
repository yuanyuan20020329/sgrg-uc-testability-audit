# -*- coding: utf-8 -*-
"""
Gemma REST v2 采集器 —— 案例 A 的数据获取层（替代 NCBI GEO，后者在本沙箱被域级阻断）。

提供：
  list_datasets(query)            -> 数据集清单（含 batchEffect / bioAssayCount）
  get_dataset(did)                -> 单数据集详情
  all_samples(did)                -> 全部样本（游标分页），含 characteristics + sample.name
  expression(did)                 -> 表达矩阵（TSV，需 gemma 的 /datasets/{id}/data 端点）
  sample_labels(did)              -> 从 characteristics 提取 (sample, 分组标签, 患者配对键)

全部响应可能是 gzip —— 统一先解压再解析。
★ Gemma 分页用游标（nextCursor），不接受 offset（会 400）。
"""
import urllib.request, urllib.parse, json, gzip, time, os, io, re, collections

BASE = "https://gemma.msl.ubc.ca/rest/v2"
UA = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
CACHE = r"C:\Users\H-P\WorkBuddy\应激颗粒\data\raw\caseA\gemma"
os.makedirs(CACHE, exist_ok=True)


def _fetch(url, binary=False, tries=5, timeout=180):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            raw = urllib.request.urlopen(req, timeout=timeout).read()
            if raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            return raw if binary else raw.decode("utf-8", "replace")
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise last


def getj(path_or_url, **params):
    url = path_or_url if path_or_url.startswith("http") else BASE + path_or_url
    if params:
        url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
    return json.loads(_fetch(url))


def list_datasets(query, limit=100, max_pages=40):
    """遍历全部命中（游标分页）。"""
    out, cursor = [], None
    for _ in range(max_pages):
        p = {"query": query, "limit": limit}
        if cursor:
            p["cursor"] = cursor
        d = getj("/datasets", **p)
        out.extend(d.get("data", []))
        cursor = d.get("nextCursor")
        if not cursor or not d.get("data"):
            break
    return out


def get_dataset(did):
    d = getj("/datasets/%s" % did)
    return (d.get("data") or [None])[0]


def all_samples(did, limit=100, max_pages=60):
    """游标分页取全部样本。"""
    out, cursor = [], None
    for _ in range(max_pages):
        p = {"limit": limit}
        if cursor:
            p["cursor"] = cursor
        d = getj("/datasets/%s/samples" % did, **p)
        batch = d.get("data", [])
        out.extend(batch)
        cursor = d.get("nextCursor")
        if not cursor or not batch:
            break
    return out


def sample_labels(did, key_categories=("tissue", "BioSource", "disease", "disease state")):
    """
    返回 list[dict(sample_id, shortName, label, patient_key, raw_chars)]。
    patient_key 从 sample.name（形如 GSE36376_Biomat_49）抽取，用于配对识别。
    """
    rows = []
    for s in all_samples(did):
        smp = s.get("sample") or {}
        chars = {}
        for c in smp.get("characteristics") or []:
            chars.setdefault(c.get("category"), []).append(c.get("value"))
        label = None
        for k in key_categories:
            if k in chars:
                label = ";".join(str(v) for v in chars[k])
                break
        nm = smp.get("name") or s.get("name") or ""
        # 患者键：去掉 GSE 前缀后取实体段，如 "GSE36376_Biomat_49" -> "49"
        m = re.search(r"_([^_]+)_(\d+)$", nm)
        pkey = m.group(2) if m else nm
        rows.append({
            "sample_id": s.get("id"),
            "shortName": s.get("shortName"),
            "name": nm,
            "label": label,
            "patient_key": pkey,
            "raw_chars": chars,
            "platform": (s.get("arrayDesign") or {}).get("shortName"),
        })
    return rows


def download_expression(did, dest=None, allow_small=False):
    """
    Gemma 数据集表达矩阵。默认走 /datasets/{id}/data?fileType=data，
    返回 TSV（行=探针/基因，列=样本）。
    """
    if dest is None:
        dest = os.path.join(CACHE, "dataset_%s_expr.tsv" % did)
    url = "%s/datasets/%s/data?fileType=data" % (BASE, did)
    try:
        txt = _fetch(url, timeout=600)
    except Exception as e:
        return None, "ERR %r" % (repr(e)[:150])
    if len(txt) < 200 and not allow_small:
        return None, "TOO-SMALL(%d bytes): %s" % (len(txt), txt[:150])
    with open(dest, "w", encoding="utf-8") as f:
        f.write(txt)
    return dest, "ok %d bytes" % len(txt)


if __name__ == "__main__":
    import sys
    did = int(sys.argv[1]) if len(sys.argv) > 1 else 7673
    d = get_dataset(did)
    print("%s | n=%s | %s" % (d.get("shortName"), d.get("bioAssayCount"), d.get("batchEffect")))
    rows = sample_labels(did)
    print("samples parsed:", len(rows))
    cnt = collections.Counter(r["label"] for r in rows)
    print("labels:", dict(cnt))

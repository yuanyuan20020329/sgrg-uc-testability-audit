# -*- coding: utf-8 -*-
"""
Part A 共享工具：Gemma 矩阵解析、列名→分组映射、SGRG 基因集
所有 Step 2 相关脚本统一从这里取数据，避免重复实现导致的口径不一致。
"""
import os, re, json
import numpy as np

BASE = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RAW  = os.path.join(BASE, "data", "raw")
OUT  = os.path.join(BASE, "results")
SGRG_CLEAN  = os.path.join(OUT, "sgrg_gene_set_clean.tsv")
MAT_38713   = os.path.join(RAW, "gemma_GSE38713.tsv")
MAT_11223   = os.path.join(RAW, "gemma_GSE11223.tsv")
MAT_22619   = os.path.join(RAW, "gemma_GSE22619.tsv")
MAT_22619_PARTIAL = os.path.join(RAW, "gemma_GSE22619.partial.tsv")
DESIGN_38713 = os.path.join(RAW, "gemma_5873_design.json")
DESIGN_11223 = os.path.join(RAW, "gemma_8385_design.json")
DESIGN_22619 = os.path.join(RAW, "gemma_6731_design.json")
SAMPLES_22619 = os.path.join(RAW, "gemma_6731_samples.json")

def log(*a): print(*a, flush=True)

# ------------------------------------------------------------------ SGRG
def load_sgrg():
    """{symbol: evidence}"""
    d = {}
    with open(SGRG_CLEAN, encoding="utf-8") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                d[p[0]] = p[1]
    return d

# ------------------------------------------------------------------ 矩阵
def load_matrix(path):
    header, rows = None, []
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts; continue
            if len(parts) < 4:
                continue
            rows.append(parts)
    return header, rows

def find_sym_col(header):
    for i, h in enumerate(header[:8]):
        if h.strip().lower() in ("genesymbol", "gene_symbol", "symbol", "gene symbols"):
            return i
    return 2

def col_keys(c):
    """一列样本名可能对应多个可用键：原始串 / '___' 前的 Biomat 名 / 'Name=' 后的样本名"""
    ks = [c.strip()]
    if "___" in c:
        ks.append(c.split("___")[0].strip())
    m = re.search(r"Name=([^=]*)$", c)
    if m:
        ks.append(m.group(1).strip())
    return ks

def resolve_label(keys, lab2g):
    """按 精确 → 尾部 '.' 还原 '-' → 子串 的顺序匹配"""
    for k in keys:
        if k in lab2g:
            return lab2g[k]
    for k in keys:
        a = re.sub(r"\.$", "-", k)
        if a in lab2g:
            return lab2g[a]
    for k in keys:
        for kk, v in lab2g.items():
            if kk and (kk == k or kk in k or k in kk):
                return v
    return "Unknown"

def collapse_genes(header, rows, sym_col, sample_cols):
    """基因级聚合（多基因共享探针同时登记到每个候选符号，取均值）"""
    gv, gp = {}, {}
    for r in rows:
        raw = (r[sym_col] or "").strip()
        if not raw or raw == "NA":
            continue
        cands = [x.strip() for x in raw.split("|") if x.strip()]
        try:
            v = np.array([float(r[c]) for c in sample_cols], float)
        except (ValueError, IndexError):
            continue
        for c in cands:
            gv.setdefault(c, []).append(v)
            gp.setdefault(c, []).append(len(cands))
    expr = {g: np.nanmean(np.vstack(vs), axis=0) for g, vs in gv.items()}
    return expr, gp

# ------------------------------------------------------------------ 分组（一律取自 design 的权威因子）
def _design_map(path, group_fn):
    d = json.load(open(path, encoding="utf-8"))["data"]
    fv = {}
    for f in d["experimentalFactors"]:
        for v in (f.get("values") or []):
            fv[v["id"]] = (f["name"], (v.get("summary") or "").strip())
    out = {}
    for a in d.get("bioMaterialAssignments", []):
        vals = {}
        for fid in (a.get("factorValueIds") or []):
            nm, sm = fv.get(fid, ("?", "?"))
            vals[nm] = sm
        key = (a.get("bioMaterialName") or "").strip()
        out[key] = group_fn(vals)
    return out

def groups_gse38713():
    """DiseaseState: reference / inactive UC / UC / UC with non-involved；另带 batch 因子"""
    def g(v):
        ds = v.get("DiseaseState", "")
        if ds.startswith("reference"):        return "Control"
        if ds.startswith("inactive"):         return "UC_remission"
        if "non-involved" in ds:              return "UC_active_noninflamed"
        if ds.startswith("ulcerative"):       return "UC_active_inflamed"
        return "Unknown"
    return _design_map(DESIGN_38713, g)

def batch_gse38713():
    """bioMaterialName -> batch 标签"""
    return _design_map(DESIGN_38713, lambda v: v.get("batch", "?"))

def groups_gse11223():
    """disease_1: ulcerative colitis vs reference；严格对照要求 disease_2 也是 reference"""
    def g(v):
        if v.get("disease_1") == "ulcerative colitis":
            return "UC"
        if v.get("disease_1") == "reference subject role" and \
           v.get("disease_2") == "reference subject role":
            return "Control"
        return "Excluded"
    return _design_map(DESIGN_11223, g)

# --------------------------------------------------- GSE22619（UC 不一致型双生子，批次平衡）
def groups_gse22619():
    """disease: ulcerative colitis（患病同胞） vs reference subject role（健康同胞）"""
    def g(v):
        ds = v.get("disease", "")
        if ds.startswith("ulcerative"):  return "UC"
        if ds.startswith("reference"):   return "Control"
        return "Unknown"
    return _design_map(DESIGN_22619, g)

def batch_gse22619():
    """bioMaterialName -> batch 标签（Batch_01/02/03）"""
    return _design_map(DESIGN_22619, lambda v: v.get("batch", "?"))

def pairs_gse22619():
    """bioMaterialName -> 双生子对编号（如 '#12'），用于配对差。
    取自 /samples 的 sample.characteristics（category='twin pair', value='#12'）。
    该队列每个批次都含完整配对，因此配对差可直接消掉批次效应。
    注：只返回对编号字符串，避免 resolve_label 把元组摊成 2D 数组。
    """
    d = json.load(open(SAMPLES_22619, encoding="utf-8"))
    out = {}
    for it in d.get("data", []):
        s = it.get("sample") or {}
        nm = (s.get("name") or "").strip()
        if not nm:
            continue
        pair = "?"
        for c in (s.get("characteristics") or []):
            if (c.get("category") or "").strip().lower() == "twin pair":
                pair = (c.get("value") or "?").strip()
                break
        out[nm] = pair
    return out

# ------------------------------------------------------------------ 矩阵完整性守卫
def decompress_gemma(gz_path, tsv_path, allow_partial=False):
    """Gemma 的 /data?download=true 返回 '先 10 字节 ' + gzip 流。通道会在中途掐断，
    因此必须显式校验 gzip 流是否正常结束（eof），否则会把截断矩阵当完整数据用。
    allow_partial=True 时：流未结束时，仍把「列数完整的行」写成矩阵（按行截断的前缀），
    并在返回值里标注 TRUNCATED_PREFIX —— 调用方必须把这一状态带进结论。
    返回 dict(status, bytes, lines, cols)。
    """
    import zlib
    if not os.path.exists(gz_path):
        return dict(status="MISSING")
    data = open(gz_path, "rb").read()
    d = zlib.decompressobj(-zlib.MAX_WBITS)
    try:
        out = d.decompress(data[10:])
        out += d.flush()
    except Exception as e:
        return dict(status="CORRUPT", bytes=len(data), error=str(e))

    def _shape(txt):
        lines = txt.split("\n")
        hdr = None
        for l in lines:
            if l.startswith("#") or not l.strip():
                continue
            hdr = l.rstrip().split("\t"); break
        ncol = len(hdr) if hdr else 0
        complete = sum(1 for l in lines[8:] if l.strip() and l.count("\t") + 1 == ncol)
        return len(lines), ncol, complete

    if d.eof:
        open(tsv_path, "wb").write(out)
        n, ncol, complete = _shape(out.decode("utf-8", errors="replace"))
        return dict(status="OK", bytes=len(data), decompressed=len(out),
                    lines=n, cols=ncol, complete_rows=complete)

    n, ncol, complete = _shape(out.decode("utf-8", errors="replace"))
    info = dict(status="TRUNCATED", bytes=len(data), decompressed=len(out),
                lines=n, cols=ncol, complete_rows=complete,
                ratio=round(len(out) / max(len(data), 1), 2))
    if allow_partial:
        # 丢掉最后一条不完整行，落盘为前缀矩阵
        txt = out.decode("utf-8", errors="replace")
        keep = [l for l in txt.split("\n")[:-1]
                if l.startswith("#") or not l.strip() or l.count("\t") + 1 == ncol]
        open(tsv_path, "w", encoding="utf-8").write("\n".join(keep) + "\n")
        info["status"] = "TRUNCATED_PREFIX"
        info["written"] = tsv_path
        info["rows_written"] = complete
    return info


def matrix_ok(path):
    """已解压矩阵的健全性：首行表头存在且数据行列数一致。"""
    if not os.path.exists(path):
        return False, "missing"
    header, rows = load_matrix(path)
    if header is None or len(rows) < 1000:
        return False, f"rows={0 if header is None else len(rows)}"
    want = len(header)
    good = sum(1 for r in rows if len(r) == want)
    if good < len(rows) * 0.5:
        return False, f"ragged ({good}/{len(rows)})"
    return True, f"rows={len(rows)} cols={want}"


# ------------------------------------------------------------------ 统计
def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    o = np.argsort(p)
    r = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    q = np.empty(n); q[o] = np.clip(r, 0, 1)
    return q
# -*- coding: utf-8 -*-
"""A2 第一步：抓取 Gemma 中全部 UC 相关队列的质控元数据 + 存算结果（不下载表达矩阵）。

产出：
  data/raw/gemma_a2/list_*.json   队列列表分页（原始）
  data/raw/gemma_a2/an_{id}.json  存算差异表达分析（含基因/探针数、DE 数、用到的因子）
  data/raw/gemma_a2/svd_{id}.json SVD（PC 方差 + 样本载荷）——用于算"主轴可分性"
  data/raw/gemma_a2/samp_{id}.json 样本特征（用于给 SVD 载荷贴疾病/批次标签）
  results/partA2_cohorts_raw.tsv  汇总宽表

用法: python partA2_harvest_meta.py [--no-svd] [--threads 6]
"""
import argparse
import gzip
import io
import json
import os
import ssl
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
CACHE = os.path.join(ROOT, "data", "raw", "gemma_a2")
os.makedirs(CACHE, exist_ok=True)
BASE = "https://gemma.msl.ubc.ca/rest/v2"
QUERIES = ["ulcerative colitis", "colitis"]

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

_lock = threading.RLock()  # ★ 必须是 RLock：work() 持锁时还会调 _log()，普通 Lock 会自锁死
_stats = {"ok": 0, "err": 0}
_tls = threading.local()


def _log(msg):
    with _lock:
        print(msg, flush=True)


def _opener():
    """每线程一个 opener（urllib 非线程安全）"""
    if not hasattr(_tls, "op"):
        _tls.op = urllib.request.build_opener(urllib.request.HTTPSHandler(context=CTX))
    return _tls.op


def fetch(url, tag, retries=3):
    """带缓存的 urllib 抓取（显式 gunzip）。返回 dict 或 None。

    注意：曾经用 subprocess+curl，在后台 shell 里会死锁（Windows 管道继承），
    本沙箱内一律改用进程内 urllib。
    """
    p = os.path.join(CACHE, tag + ".json")
    if os.path.exists(p) and os.path.getsize(p) > 2:
        try:
            return json.load(open(p, encoding="utf-8"))
        except Exception:
            try:
                os.remove(p)
            except Exception:
                pass
    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0 (research; QC audit)",
        "Accept": "application/json",
        "Accept-Encoding": "gzip",
    })
    for k in range(retries):
        try:
            with _opener().open(req, timeout=90) as r:
                raw = r.read()
                if r.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.decompress(raw)
            d = json.loads(raw.decode("utf-8"))
            if isinstance(d, dict) and d.get("error"):
                raise ValueError("api error: %s" % str(d["error"])[:120])
            with open(p, "w", encoding="utf-8") as f:
                json.dump(d, f, ensure_ascii=False)
            with _lock:
                _stats["ok"] += 1
            return d
        except Exception as e:
            if k == retries - 1:
                with _lock:
                    _stats["err"] += 1
                _log("  !! 失败 %s : %s %s" % (tag, type(e).__name__, str(e)[:80]))
                return None
            time.sleep(1.2 * (k + 1))
    return None


# ---------------------------------------------------------------- 1. 队列列表
def harvest_lists():
    ids = {}
    pages = 0
    for q in QUERIES:
        off = 0
        while True:
            tag = "list_%s_%d" % (q.replace(" ", "_"), off)
            d = fetch("%s/datasets?query=%s&limit=100&offset=%d" % (BASE, urllib.parse.quote(q), off), tag)
            if d is None:
                break
            data = d.get("data") or []
            for e in data:
                ids.setdefault(e["id"], e)
            pages += 1
            tot = d.get("totalElements") or 0
            _log("  [list] q=%s offset=%d 本页 %d 累计 %d / %d" % (q, off, len(data), len(ids), tot))
            off += 100
            if len(data) < 100 or off >= tot:
                break
    _log("[1] 队列列表：%d 个唯一数据集（%d 页）" % (len(ids), pages))
    return ids


# ---------------------------------------------------------------- 2. 逐队列
def harvest_one(did):
    fetch("%s/datasets/%s/analyses/differential" % (BASE, did), "an_%s" % did)
    if not ARGS.no_svd:
        fetch("%s/datasets/%s/svd" % (BASE, did), "svd_%s" % did)
    # 样本特征（cursor 分页）
    if not ARGS.no_svd:
        tag0 = "samp_%s_0" % did
        d = fetch("%s/datasets/%s/samples?limit=100" % (BASE, did), tag0)
        cur = (d or {}).get("nextCursor")
        n = 1
        while cur and n < 8:
            d2 = fetch("%s/datasets/%s/samples?limit=100&cursor=%s" % (BASE, did, cur), "samp_%s_%d" % (did, n))
            cur = (d2 or {}).get("nextCursor")
            n += 1


def main():
    ids = harvest_lists()
    todo = [i for i in sorted(ids)]
    if ARGS.limit:
        todo = todo[:ARGS.limit]
    _log("[2] 开始抓取 %d 个队列的存算结果（线程 %d）" % (len(todo), ARGS.threads))
    from concurrent.futures import ThreadPoolExecutor
    t0 = time.time()
    done = [0]
    def work(did):
        t1 = time.time()
        try:
            harvest_one(did)
        except Exception as e:
            _log("  !! %s %s %s" % (did, type(e).__name__, e))
        with _lock:
            done[0] += 1
            _log("  [2] %d/%d gid=%s %.1fs (ok=%d err=%d)" % (done[0], len(todo), did, time.time() - t1,
                                                             _stats["ok"], _stats["err"]))
    if ARGS.threads <= 1:
        for d in todo:
            work(d)
    else:
        with ThreadPoolExecutor(max_workers=ARGS.threads) as ex:
            list(ex.map(work, todo))
    _log("[2] 抓取完成 %.0fs  ok=%d err=%d" % (time.time() - t0, _stats["ok"], _stats["err"]))
    json.dump({str(k): v for k, v in ids.items()}, open(os.path.join(CACHE, "cohort_index.json"), "w", encoding="utf-8"),
              ensure_ascii=False)
    _log("[3] 索引写出 data/raw/gemma_a2/cohort_index.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-svd", action="store_true")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ARGS = ap.parse_args()
    main()

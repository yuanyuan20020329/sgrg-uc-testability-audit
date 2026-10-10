# -*- coding: utf-8 -*-
"""GSE133060 检验的数据准备（ resilient 下载，网络限流下可跑）

目标
----
GSE133060 = 全库唯一「UC 且 ≥6 对同供体炎症/非炎症配对且测 mRNA」的公共资源
（8 对，Agilent GPL4133 全基因组表达芯片）。本脚本把它需要的两个文件 resilient 地拉下来：
  1) GSE133060_RAW.tar   42 MB   16 个 Agilent Feature Extraction 文本
  2) GPL4133_family.soft.gz      探针 → 基因符号注释

为什么不用 dl_parallel.py 多线程：NCBI 今天对手握手严重限流，并发只会更糟。
这里用单流 + 指数退避重试 + 续传（Range）。
"""
import os
import ssl
import socket
import sys
import time
import urllib.request

socket.setdefaulttimeout(180)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}

URLS = [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133060/suppl/GSE133060_RAW.tar",
     "data/raw/gse133060/GSE133060_RAW.tar"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL413nnn/GPL4133/soft/GPL4133_family.soft.gz",
     "data/raw/gse133060/GPL4133_family.soft.gz"),
]


def dl(url, out, tries=12):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    have = os.path.getsize(out) if os.path.exists(out) else 0
    # 先 HEAD 拿总大小
    total = None
    for k in range(tries):
        try:
            r = urllib.request.Request(url, method="HEAD", headers=UA)
            with urllib.request.urlopen(r, context=CTX, timeout=120) as x:
                total = int(x.headers.get("Content-Length", 0))
            break
        except Exception as e:
            wait = min(120, 10 * (2 ** min(k, 4)))
            print("HEAD 失败(%s) 第%d次，%ds 后重试" % (type(e).__name__, k + 1, wait), flush=True)
            time.sleep(wait)
    if not total:
        print("HEAD 全失败，放弃", url)
        return False
    if have >= total:
        print("[已存在] %s (%d bytes)" % (out, have))
        return True
    print("目标 %d bytes，已有 %d" % (total, have), flush=True)
    for k in range(tries):
        try:
            start = os.path.getsize(out) if os.path.exists(out) else 0
            if start >= total:
                print("[完成] %s" % out)
                return True
            req = urllib.request.Request(url, headers={**UA, "Range": "bytes=%d-" % start})
            with urllib.request.urlopen(req, context=CTX, timeout=180) as x, \
                 open(out, "ab") as f:
                while True:
                    chunk = x.read(1 << 17)
                    if not chunk:
                        break
                    f.write(chunk)
            size = os.path.getsize(out)
            print("  本轮得 %d / %d" % (size, total), flush=True)
            if size >= total:
                print("[完成] %s (%d bytes, %d 轮)" % (out, size, k + 1))
                return True
        except Exception as e:
            got = os.path.getsize(out) if os.path.exists(out) else 0
            wait = min(180, 20 * (2 ** min(k, 4)))
            print("  中断于 %d/%d (%s)，%ds 后续传" % (got, total, type(e).__name__, wait), flush=True)
            time.sleep(wait)
    return False


if __name__ == "__main__":
    ok = all(dl(u, o) for u, o in URLS)
    print("\nDONE", ok)
    sys.exit(0 if ok else 1)

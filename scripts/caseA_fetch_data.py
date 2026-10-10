#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
案例 A 数据获取：GSE65372（发现）/ GSE84402（验证）
=================================================
带 HTTP Range 续传 + 指数退避。GEO FTP 在本沙箱单流偶发超时，
故按 Range 分块续传（Part D 的 dl_parallel 同思路，但单文件串行加续传）。

输出：
  data/raw/caseA/GSE65372_series_matrix.txt.gz
  data/raw/caseA/GSE84402_series_matrix.txt.gz
"""
import os, sys, time, urllib.request

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(ROOT, "data", "raw", "caseA")
os.makedirs(DST, exist_ok=True)


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", 0))


def fetch(url, out, chunk=8 << 20, tries=12):
    total = None
    for i in range(tries):
        try:
            total = remote_size(url)
            break
        except Exception as e:
            print("  HEAD retry %d: %r" % (i, e))
            time.sleep(2 + i * 2)
    if not total:
        print("  [ERR] 无法取到远端大小")
        return 0

    have = os.path.getsize(out) if os.path.exists(out) else 0
    print("  远端 %d bytes，本地已有 %d" % (total, have))
    if have >= total:
        print("  已完成")
        return have

    mode = "ab" if have else "wb"
    f = open(out, mode)
    try:
        while have < total:
            end = min(have + chunk - 1, total - 1)
            ok = False
            for i in range(tries):
                try:
                    req = urllib.request.Request(url, headers={
                        "User-Agent": "Mozilla/5.0",
                        "Range": "bytes=%d-%d" % (have, end),
                    })
                    with urllib.request.urlopen(req, timeout=180) as r:
                        data = r.read()
                    f.write(data)
                    f.flush()
                    have += len(data)
                    ok = True
                    print("    %d / %d  (%.1f%%)" % (have, total, 100.0 * have / total))
                    break
                except Exception as e:
                    print("    retry %d @%d: %r" % (i, have, e))
                    time.sleep(2 + i * 3)
            if not ok:
                print("  [ERR] 分块失败，中止")
                break
    finally:
        f.close()
    return have


def main():
    for gse in ["GSE65372", "GSE84402"]:
        url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/%s_series_matrix.txt.gz"
               % (gse[:-3], gse, gse))
        out = os.path.join(DST, "%s_series_matrix.txt.gz" % gse)
        print("=== %s" % gse)
        n = fetch(url, out)
        print("  => %s : %d bytes" % (out, n))


if __name__ == "__main__":
    main()

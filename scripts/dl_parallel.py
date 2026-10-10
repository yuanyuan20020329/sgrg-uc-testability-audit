# -*- coding: utf-8 -*-
"""多线程分片下载大文件：每片独立 .partN 文件（可单独续传），全部完成后合并。
避免「两进程同时追加同一文件」造成的静默损坏。
用法: python dl_parallel.py <url> <out> [nthreads]
"""
import os
import sys
import ssl
import time
import threading
import urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}

url, out = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 8

os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

req0 = urllib.request.Request(url, method="HEAD", headers=UA)
with urllib.request.urlopen(req0, context=CTX, timeout=90) as r:
    total = int(r.headers.get("Content-Length", 0))
    ranges = r.headers.get("Accept-Ranges", "")
print("target = %d bytes | Accept-Ranges = %s" % (total, ranges), flush=True)
if total <= 0:
    print("无法获取文件大小"); sys.exit(1)

base = total // N
bounds = []
for i in range(N):
    s = i * base
    e = total - 1 if i == N - 1 else (i + 1) * base - 1
    bounds.append((s, e))

T0 = time.time()
lock = threading.Lock()


def worker(i, s, e):
    part = "%s.part%d" % (out, i)
    need = e - s + 1
    attempt = 0
    while True:
        done = os.path.getsize(part) if os.path.exists(part) else 0
        if done >= need:
            with lock:
                print("  part%d DONE %.1f MB" % (i, done / 1e6), flush=True)
            return
        attempt += 1
        if attempt > 60:
            with lock:
                print("  part%d ABORT %d/%d" % (i, done, need), flush=True)
            return
        try:
            h = dict(UA)
            h["Range"] = "bytes=%d-%d" % (s + done, e)
            rq = urllib.request.Request(url, headers=h)
            with urllib.request.urlopen(rq, context=CTX, timeout=300) as r:
                with open(part, "ab" if done else "wb") as f:
                    last = time.time()
                    while True:
                        c = r.read(1 << 20)
                        if not c:
                            break
                        f.write(c)
                        done += len(c)
                        if time.time() - last > 20:
                            last = time.time()
                            with lock:
                                got = sum(os.path.getsize("%s.part%d" % (out, k))
                                          for k in range(N) if os.path.exists("%s.part%d" % (out, k)))
                                print("  total %.1f / %.1f MB (%.2f MB/s)" % (
                                    got / 1e6, total / 1e6, got / 1e6 / max(time.time() - T0, 1e-6)), flush=True)
        except Exception as ex:
            with lock:
                print("  part%d 重试#%d (%s) %.1f MB" % (i, attempt, type(ex).__name__, done / 1e6), flush=True)
            time.sleep(2)


ths = [threading.Thread(target=worker, args=(i, s, e), daemon=True) for i, (s, e) in enumerate(bounds)]
for t in ths:
    t.start()
for t in ths:
    t.join()

sizes = [os.path.getsize("%s.part%d" % (out, i)) for i in range(N)]
print("分片完成:", [round(x / 1e6, 1) for x in sizes], "合计 %.1f MB" % (sum(sizes) / 1e6), flush=True)
if sum(sizes) != total:
    print("!! 分片大小不足，未合并。缺 %d 字节" % (total - sum(sizes))); sys.exit(2)

with open(out, "wb") as fo:
    for i in range(N):
        part = "%s.part%d" % (out, i)
        with open(part, "rb") as fi:
            while True:
                b = fi.read(1 << 22)
                if not b:
                    break
                fo.write(b)
if os.path.getsize(out) != total:
    print("!! 合并后大小 %d != %d" % (os.path.getsize(out), total)); sys.exit(3)
for i in range(N):
    os.remove("%s.part%d" % (out, i))
print("DONE %s  %.1f MB  %.1fs" % (out, total / 1e6, time.time() - T0), flush=True)

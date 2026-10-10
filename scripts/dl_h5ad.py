# -*- coding: utf-8 -*-
"""下载 CELLxGENE h5ad：断点续传 + 超时自动重试"""
import os
import sys
import time
import ssl
import urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

url, out = sys.argv[1], sys.argv[2]
os.makedirs(os.path.dirname(out), exist_ok=True)

# 目标大小
req0 = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
try:
    with urllib.request.urlopen(req0, context=CTX, timeout=60) as r:
        target = int(r.headers.get("Content-Length", 0))
except Exception:
    target = 0
print("target = %d bytes" % target, flush=True)

t0 = time.time()
attempt = 0
while True:
    done = os.path.getsize(out) if os.path.exists(out) else 0
    if target and done >= target:
        print("DONE %s  %.1f MB  %.1fs" % (out, done / 1e6, time.time() - t0), flush=True)
        break
    attempt += 1
    if attempt > 40:
        print("ABORT 重试过多，停在 %d / %d" % (done, target), flush=True)
        sys.exit(2)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        if done:
            req.add_header("Range", "bytes=%d-" % done)
        with urllib.request.urlopen(req, context=CTX, timeout=300) as r:
            with open(out, "ab" if done else "wb") as f:
                last = time.time()
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    done += len(chunk)
                    if time.time() - last > 15:
                        last = time.time()
                        print("  %.1f / %.1f MB  (%.2f MB/s)" % (
                            done / 1e6, target / 1e6, done / 1e6 / max(time.time() - t0, 1e-6)), flush=True)
    except Exception as e:
        print("  重试 #%d（%s: %s）已下载 %.1f MB" % (
            attempt, type(e).__name__, str(e)[:80], done / 1e6), flush=True)
        time.sleep(3)
        continue
done = os.path.getsize(out)
if target and done < target:
    print("未完成：%d / %d" % (done, target), flush=True)
    sys.exit(3)
print("DONE %s  %.1f MB  %.1fs" % (out, done / 1e6, time.time() - t0), flush=True)

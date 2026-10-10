#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parallel, resumable chunked downloader for the UC GWAS (de Lange 2017, GCST004133).

Why: this sandbox throttles each TCP connection to roughly 9-12 KB/s and does not
let a single connection finish a multi-hundred-MB transfer. ftp.ebi.ac.uk DOES
advertise `Accept-Ranges: bytes`, and throughput scales ~linearly with the number
of concurrent connections. So we split the object into fixed-size chunks, fetch
them in parallel with per-chunk resume, then concatenate.

Idempotent: re-running picks up from where it left off (complete chunk files are
skipped). Safe to interrupt.
"""
import os
import sys
import time
import json
import urllib.request
import urllib.error
import threading

URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST004001-GCST005000/GCST004133/harmonised/"
       "28067908-GCST004133-EFO_0000729.h.tsv.gz")

ROOT = r"C:/Users/H-P/WorkBuddy/应激颗粒"
RAW = os.path.join(ROOT, "data", "raw", "gwas")
OUT = os.path.join(RAW, "GCST004133_uc_harmonised.h.tsv.gz")
CHUNKDIR = os.path.join(RAW, "gcst004133_chunks")
STATE = os.path.join(RAW, "gcst004133_state.json")
LOG = os.path.join(ROOT, "logs", "download_gcst004133.log")

CHUNK = 2 * 1024 * 1024          # 2 MB
NWORKERS = 16
PER_CONN_TIMEOUT = 120           # seconds per HTTP read attempt
MAX_TRIES_PER_CHUNK = 40
UA = "Mozilla/5.0 (compatible; research-script/1.0)"


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    try:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except Exception:
        pass


def head_total():
    req = urllib.request.Request(URL, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def fetch_range(start, end, dest, expected):
    """Download [start,end] inclusive into dest, resuming a partial dest file."""
    tries = 0
    while True:
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if have >= expected:
            return True
        tries += 1
        if tries > MAX_TRIES_PER_CHUNK:
            return False
        s = start + have
        req = urllib.request.Request(
            URL, headers={"User-Agent": UA, "Range": "bytes=%d-%d" % (s, end)})
        try:
            with urllib.request.urlopen(req, timeout=PER_CONN_TIMEOUT) as r:
                mode = "ab" if have else "wb"
                with open(dest, mode) as f:
                    while True:
                        try:
                            b = r.read(32768)
                        except Exception:
                            break
                        if not b:
                            break
                        f.write(b)
        except Exception as e:
            log("  chunk retry (%d bytes so far): %s" % (have, type(e).__name__))
            time.sleep(1.5)
            continue
        if os.path.getsize(dest) >= expected:
            return True
        time.sleep(0.5)
    return False


def main():
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(CHUNKDIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)

    total = head_total()
    nchunks = (total + CHUNK - 1) // CHUNK
    log("URL OK; total=%d bytes; chunks=%d (%.0f MB each)" % (total, nchunks, CHUNK / 1e6))

    # idempotency / consistency guard
    if os.path.exists(STATE):
        st = json.load(open(STATE))
        if st.get("total") != total:
            log("WARNING: size changed since last run; resetting chunk dir")
            for f in os.listdir(CHUNKDIR):
                os.remove(os.path.join(CHUNKDIR, f))
        else:
            log("resuming previous run")
    json.dump({"url": URL, "total": total, "chunk": CHUNK, "n": nchunks},
              open(STATE, "w"), indent=1)

    done = [i for i in range(nchunks)
            if os.path.exists(os.path.join(CHUNKDIR, "chunk_%05d" % i))
            and os.path.getsize(os.path.join(CHUNKDIR, "chunk_%05d" % i)) ==
                min(CHUNK, total - i * CHUNK)]
    log("already complete: %d / %d" % (len(done), nchunks))

    todo = [i for i in range(nchunks) if i not in set(done)]
    lock = threading.Lock()
    counter = {"done": len(done), "bytes": sum(
        os.path.getsize(os.path.join(CHUNKDIR, "chunk_%05d" % i)) for i in done)}
    t0 = time.time()

    def worker():
        while True:
            with lock:
                if not todo:
                    return
                ci = todo.pop(0)
            start = ci * CHUNK
            end = min(start + CHUNK, total) - 1
            expected = end - start + 1
            dest = os.path.join(CHUNKDIR, "chunk_%05d" % ci)
            ok = fetch_range(start, end, dest, expected)
            with lock:
                if ok:
                    counter["done"] += 1
                    counter["bytes"] += expected
                    el = time.time() - t0
                    spd = counter["bytes"] / max(el, 1) / 1024.0
                    log("chunk %d/%d done | %d/%d | %.0f KB/s | eta %.1f min"
                        % (ci, nchunks, counter["done"], nchunks, spd,
                           (total - counter["bytes"]) / max(spd * 1024, 1) / 60))
                else:
                    todo.append(ci)
                    log("chunk %d FAILED after retries, requeued" % ci)

    threads = [threading.Thread(target=worker) for _ in range(NWORKERS)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    missing = [i for i in range(nchunks)
               if not os.path.exists(os.path.join(CHUNKDIR, "chunk_%05d" % i))
               or os.path.getsize(os.path.join(CHUNKDIR, "chunk_%05d" % i)) !=
               min(CHUNK, total - i * CHUNK)]
    if missing:
        log("INCOMPLETE: %d chunks missing: %s" % (len(missing), missing[:20]))
        sys.exit(2)

    log("concatenating ...")
    with open(OUT, "wb") as out:
        for i in range(nchunks):
            with open(os.path.join(CHUNKDIR, "chunk_%05d" % i), "rb") as f:
                while True:
                    b = f.read(1 << 20)
                    if not b:
                        break
                    out.write(b)
    sz = os.path.getsize(OUT)
    log("ASSEMBLED %s  (%d bytes; expected %d)" % (OUT, sz, total))
    if sz != total:
        log("SIZE MISMATCH")
        sys.exit(3)
    log("DONE in %.1f min" % ((time.time() - t0) / 60))


if __name__ == "__main__":
    main()

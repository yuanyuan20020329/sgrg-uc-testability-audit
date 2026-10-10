# -*- coding: utf-8 -*-
"""
coloc.susie 的 LD 参考构建 —— 从 1000G phase3 (GRCh37) 按区域抽取 EUR 基因型并算 LD。

为什么要自己写：本机无 plink / bcftools / tabix，单个染色体 VCF 近 1 GB，
而候选区域只有 ~1.1 Mb。做法是解析 .tbi 索引 + HTTP Range 只取覆盖该区域的 BGZF 块。

用法：
    python fetch_regions.py --budget        # 只打印各区域的下载预算
    python fetch_regions.py --run           # 实际下载并抽取
"""
import argparse, gzip, json, os, subprocess, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ld_bgzf_tabix import read_tbi, region_file_range, decompress_range

BASE_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
            "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
PANEL_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
             "integrated_call_samples_v3.20130502.ALL.panel")

WORK = r"D:\ld_ref_1000g"
OUT = r"D:\ld_ref_1000g\regions"

# 与 09_coloc_fixed 的区域定义完全一致（GRCh37，1-based 闭区间）
REGIONS = [
    ("chr16_ogfod1", 16, 55000000, 58000000),
]
EUR_POPS = {"CEU", "TSI", "FIN", "GBR", "IBS"}


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, **kw)


def curl_range(url, start, stop, out_path, timeout=1800):
    """下载 [start, stop] 字节区间。"""
    tmp = out_path + ".part"
    r = sh("curl -s --fail --max-time %d -r %d-%d -o \"%s\" \"%s\""
           % (timeout, start, stop, tmp, url))
    if r.returncode != 0 or not os.path.exists(tmp):
        raise RuntimeError("range 下载失败 %s [%d,%d]: %s" % (url, start, stop, r.stderr[:300]))
    os.replace(tmp, out_path)
    return os.path.getsize(out_path)


def ensure_tbi(chrom):
    p = os.path.join(WORK, "chr%d.tbi" % chrom)
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return p
    url = BASE_URL.format(chr=chrom) + ".tbi"
    r = sh("curl -s --fail --max-time 300 -o \"%s\" \"%s\"" % (p, url))
    if r.returncode != 0:
        raise RuntimeError("tbi 下载失败: " + url)
    return p


def ensure_panel():
    """EUR 样本列表。★ panel 文件列序为 sample / pop / super_pop / gender，
    故人群列是索引 1（曾误用索引 2 的 super_pop，导致 EUR 列表恒为空）。"""
    p = os.path.join(WORK, "integrated_call_samples_v3.20130502.ALL.panel")
    if not os.path.exists(p):
        sh("curl -s --fail --max-time 300 -o \"%s\" \"%s\"" % (p, PANEL_URL))
    eur = []
    with open(p, encoding="utf-8") as fh:
        fh.readline()                               # 跳过表头
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1] in EUR_POPS:
                eur.append(f[0])
    if not eur:
        raise RuntimeError("EUR 样本列表为空 —— 检查 panel 文件列序：%s" % p)
    return eur


def ensure_samples(chrom):
    """样本名列表。区域 chunk 是文件中段，不含 #CHROM 行，故须单独取文件头。

    1000G 各染色体的样本集合完全一致，故取一次后缓存复用。
    """
    p = os.path.join(WORK, "vcf_samples.txt")
    if os.path.exists(p):
        return open(p, encoding="utf-8").read().rstrip("\n").split("\t")
    url = BASE_URL.format(chr=chrom)
    tmp = os.path.join(WORK, "vcf_head_probe.gz")
    curl_range(url, 0, 400000, tmp)
    txt = decompress_range(open(tmp, "rb").read(), 0).decode("utf-8", "replace")
    os.remove(tmp)
    for line in txt.split("\n"):
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            with open(p, "w", encoding="utf-8") as fh:
                fh.write("\t".join(samples))
            return samples
    raise RuntimeError("未能从 VCF 头部取到 #CHROM 行")


def fetch_region(chrom, start, end):
    """返回 (samples, list_of_records, bytes_downloaded)；records = (pos, ref, alt, gt_list)。"""
    tbi = read_tbi(ensure_tbi(chrom))
    ref_name = str(chrom) if str(chrom) in tbi["bins"] else tbi["names"][0]
    rng = region_file_range(tbi, ref_name, start - 1, end)
    if rng is None:
        raise RuntimeError("无法定位区域 %d:%d-%d" % (chrom, start, end))
    raw_path = os.path.join(WORK, "chunk_chr%d_%d_%d.gz" % (chrom, start, end))
    n = curl_range(BASE_URL.format(chr=chrom), rng[0], rng[1] - 1, raw_path)
    data = open(raw_path, "rb").read()
    txt = decompress_range(data, 0).decode("utf-8", "replace")

    # ★ 区域 chunk 取自文件中段，不含 #CHROM 行 → 样本名单独取，否则会静默拿到空列表
    samples = ensure_samples(chrom)
    recs = []
    for line in txt.split("\n"):
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or len(f) != len(samples) + 9:
            continue
        try:
            pos = int(f[1])
        except ValueError:
            continue
        if pos < start or pos > end:
            continue
        alt = f[4]
        if "," in alt:          # 只保留双等位
            continue
        recs.append((pos, f[3].upper(), alt.upper(), f[9:], f[2]))
    try:
        os.remove(raw_path)
    except OSError:
        pass
    if not recs:
        raise RuntimeError("区域 %d:%d-%d 抽取到 0 条记录（range=%s）" % (chrom, start, end, rng))
    return samples, recs, n


def gt_dosage(gt):
    """'0|1' -> 1（ALT 拷贝数）；缺失返回 None。"""
    if not gt or gt[0] in ".|":
        return None
    a = gt.split(":")[0].replace("|", "/")
    parts = a.split("/")
    if len(parts) != 2 or "." in parts:
        return None
    try:
        return int(parts[0]) + int(parts[1])
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", action="store_true")
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    if args.budget or not args.run:
        print("%-16s %-6s %-12s %-12s %s" % ("region", "chr", "start", "end", "预算下载量"))
        tot = 0
        for rid, c, s, e in REGIONS:
            tbi = read_tbi(ensure_tbi(c))
            rng = region_file_range(tbi, str(c), s - 1, e)
            mb = (rng[1] - rng[0]) / 1e6
            tot += mb
            print("%-16s %-6d %-12d %-12d %.2f MB" % (rid, c, s, e, mb))
        print("合计 %.2f MB" % tot)
        if not args.run:
            return

    eur = ensure_panel()
    print("EUR 样本数 =", len(eur), flush=True)

    manifest = []
    for rid, c, s, e in REGIONS:
        t0 = time.time()
        samples, recs, nbytes = fetch_region(c, s, e)
        idx = [samples.index(x) for x in eur if x in samples]
        if len(idx) < 0.5 * len(eur):
            raise RuntimeError("%s: EUR 样本仅匹配到 %d/%d —— 样本名对不上"
                               % (rid, len(idx), len(eur)))
        eur_used = [samples[i] for i in idx]
        out = os.path.join(OUT, "%s.tsv" % rid)
        nsnv = 0
        with open(out, "w", newline="", encoding="utf-8") as fh:
            fh.write("\t".join(["pos", "ref", "alt", "rsid"] + eur_used) + "\n")
            for pos, r, a, gts, rsid in recs:
                ds = [gt_dosage(gts[i]) for i in idx]
                ok = sum(1 for d in ds if d is not None)
                if ok < 0.8 * len(ds):     # 缺失率 >20% 丢弃
                    continue
                ac = sum(d for d in ds if d is not None)
                af = ac / (2 * ok)
                if af < 0.01 or af > 0.99:  # 单态/近单态对 LD 无贡献
                    continue
                fh.write("\t".join([str(pos), r, a, rsid] + ["." if d is None else str(d) for d in ds]) + "\n")
                nsnv += 1
        manifest.append(dict(region=rid, chr=c, start=s, end=e,
                             downloaded_MB=round(nbytes / 1e6, 2),
                             n_eur=len(eur_used), n_snv=nsnv, seconds=round(time.time() - t0, 1)))
        print("  %-16s 下载 %.2f MB | EUR=%d | 保留 SNP=%d | %.1fs"
              % (rid, nbytes / 1e6, len(eur_used), nsnv, time.time() - t0), flush=True)

    json.dump(manifest, open(os.path.join(OUT, "fetch_manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n完成 →", OUT)


if __name__ == "__main__":
    main()

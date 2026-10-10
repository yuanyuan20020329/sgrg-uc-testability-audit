# -*- coding: utf-8 -*-
"""
BGZF + Tabix(.tbi) 最小实现 —— 用于对 1000G VCF 做 HTTP range 按区域抽取。

设计动机：本机无 plink / bcftools / tabix，且 1000G 单条染色体 VCF 近 1 GB，
全量下载不可接受。改为
  1) 下载体积很小的 .tbi 索引（~170 KB）
  2) 解析 binning index，把目标区域映射到少量 BGZF 压缩块的虚拟偏移
  3) 用 HTTP Range 只下载这些块，就地流式解压

依赖：仅标准库 + 已装 curl（用于 range 请求，由调用方负责）。
"""
import struct, zlib, io, os

BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


# ----------------------------------------------------------------------
# BGZF
# ----------------------------------------------------------------------
def bgzf_blocks(data, start=0):
    """切分 BGZF 流为块列表 [(file_offset, compressed_bytes)]。

    遇到不完整块（HTTP range 末尾截断）或 EOF 块即停止。
    """
    blocks = []
    pos = start
    n = len(data)
    while pos + 12 <= n:
        if data[pos] != 0x1F or data[pos + 1] != 0x8B:
            break
        xlen = struct.unpack("<H", data[pos + 10:pos + 12])[0]
        if pos + 12 + xlen > n:
            break
        extra = data[pos + 12:pos + 12 + xlen]
        bsize = None
        k = 0
        while k + 4 <= len(extra):
            si1, si2, slen = extra[k], extra[k + 1], struct.unpack("<H", extra[k + 2:k + 4])[0]
            if si1 == 66 and si2 == 67 and slen == 2:
                bsize = struct.unpack("<H", extra[k + 4:k + 6])[0]
                break
            k += 4 + slen
        if bsize is None:
            break
        total = bsize + 1
        if pos + total > n:
            break  # 截断
        blocks.append((pos, data[pos:pos + total]))
        pos += total
    return blocks


def bgzf_decompress(block):
    """解压单个 BGZF 块，返回 bytes（失败返回 None）。

    ★ 两个坑（都真实踩过）：
      1) 不能把整个块丢给 zlib —— 块首是 18 字节 gzip header，zlib 会把首字节
         0x1f 当成 deflate 块头（BTYPE=3，非法块类型）而报 'invalid block type'。
         必须从 12+XLEN 处开始，到末尾 8 字节 CRC32+ISIZE 之前结束。
      2) 用 decompressobj 而非 decompress，前者不会因尾部残留而报错。
    """
    try:
        xlen = struct.unpack("<H", block[10:12])[0]
        body = block[12 + xlen:-8]
        return zlib.decompressobj(-15).decompress(body)
    except (zlib.error, struct.error):
        return None


def decompress_range(data, start_offset):
    """把一段从 start_offset 起的 BGZF 字节解压成文本行列表。"""
    out = []
    for off, blk in bgzf_blocks(data):
        if blk == BGZF_EOF:
            break
        dec = bgzf_decompress(blk)
        if dec is None:
            continue
        out.append(dec)
    return b"".join(out)


# ----------------------------------------------------------------------
# Tabix
# ----------------------------------------------------------------------
def read_tbi(path):
    """解析 .tbi 索引，返回 dict：names, bins[refname]->{bin: [(voff_beg, voff_end)]},
    linear[refname]->[voff,...], n_ref, col_seq/col_beg/col_end 等 meta。"""
    raw = open(path, "rb").read()
    # 1000G 分发的 .tbi 自身被 gzip 压缩（部分是 BGZF，部分是普通 gzip）→ 先解开
    if raw[:2] == b"\x1f\x8b":
        try:
            import gzip as _gz
            raw = _gz.decompress(raw)
        except Exception:
            raw = decompress_range(raw, 0)
    if raw[:4] != b"TBI\x01":
        raise ValueError("不是 TBI 索引（magic 不符）: %r" % raw[:4])
    p = 4
    (n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm) = struct.unpack("<8i", raw[p:p + 32])
    p += 32
    names = raw[p:p + l_nm].split(b"\x00")
    p += l_nm
    names = [x.decode() for x in names if x]

    bins_all, linear_all = {}, {}
    for ref in range(n_ref):
        n_bin = struct.unpack("<i", raw[p:p + 4])[0]; p += 4
        bins = {}
        for _ in range(n_bin):
            bin_id = struct.unpack("<I", raw[p:p + 4])[0]; p += 4
            n_chunk = struct.unpack("<i", raw[p:p + 4])[0]; p += 4
            chunks = []
            for _ in range(n_chunk):
                voff_beg, voff_end = struct.unpack("<QQ", raw[p:p + 16]); p += 16
                chunks.append((voff_beg, voff_end))
            bins[bin_id] = chunks
        n_intv = struct.unpack("<i", raw[p:p + 4])[0]; p += 4
        lin = list(struct.unpack("<%dQ" % n_intv, raw[p:p + 8 * n_intv])); p += 8 * n_intv
        rname = names[ref] if ref < len(names) else str(ref)
        bins_all[rname] = bins
        linear_all[rname] = lin
    return dict(names=names, bins=bins_all, linear=linear_all,
                col_seq=col_seq, col_beg=col_beg, col_end=col_end, meta=meta, skip=skip)


def reg2bins(beg, end):
    """Tabix 的 region→bins 映射（beg/end 为 0-based 半开区间）。"""
    out = [0]
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1):
        out.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1):
        out.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1):
        out.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1):
        out.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1):
        out.append(k)
    return out


def chunks_for_region(tbi, ref, beg, end):
    """返回覆盖 [beg,end) 的最小/最大压缩文件偏移（用于 HTTP Range）。"""
    bins = tbi["bins"].get(ref)
    if bins is None:
        return None
    want = set(reg2bins(beg, end))
    (vmin, vmax) = (1 << 48, 0)
    for b in want:
        for (vb, ve) in bins.get(b, []):
            # 只保留与线性索引区间有交集的 chunk
            vmin = min(vmin, vb)
            vmax = max(vmax, ve)
    if vmax == 0:
        return None
    # 虚拟偏移 = (coff<<16) | uoff
    return (vmin >> 16, ((vmax >> 16) + 65536))


def linear_min_offset(tbi, ref, beg):
    """用线性索引给出该位置之前一定落在目标区域的起始偏移（粗切）。"""
    lin = tbi["linear"].get(ref)
    if not lin:
        return None
    idx = beg >> 14
    for i in range(idx, min(idx + 64, len(lin))):
        if lin[i]:
            return lin[i] >> 16
    return None


def region_file_range(tbi, ref, beg, end, pad=131072):
    """给出覆盖 [beg,end) 所需的最小文件字节区间 (start, stop)。

    以 linear index 为主（每个 16 kb 窗口首记录的虚拟偏移，粒度远细于 binning index）：
      下界 = 起点所在窗口起第一个非零 linear（其所在块一定不晚于区域首记录）
      上界 = 终点窗口之后第一个非零 linear（其后已无区域记录）
    beg/end 为 0-based 半开区间。返回的 stop 再放宽一个块长以免截断。
    """
    lin = tbi["linear"].get(ref)
    bins = tbi["bins"].get(ref)
    if not lin:
        c = chunks_for_region(tbi, ref, beg, end)
        return c
    i0, i1 = beg >> 14, end >> 14
    n = len(lin)

    start = None
    for i in range(i0, min(i1 + 1, n)):
        if lin[i]:
            start = lin[i] >> 16
            break
    if start is None:  # 区间内全空 → 向前找最近的锚点
        for i in range(i0 - 1, max(-1, i0 - 4096), -1):
            if i >= 0 and lin[i]:
                start = lin[i] >> 16
                break
    if start is None and bins:
        vmins = [vb for b in reg2bins(beg, end) for (vb, ve) in bins.get(b, [])]
        start = (min(vmins) >> 16) if vmins else None
    if start is None:
        return None

    stop = None
    for i in range(i1 + 1, min(i1 + 4096, n)):
        if lin[i]:
            stop = lin[i] >> 16
            break
    if stop is None and bins:
        vmaxs = [ve for b in reg2bins(beg, end) for (vb, ve) in bins.get(b, [])]
        if vmaxs:
            stop = (max(vmaxs) >> 16) + 65536
    if stop is None:
        stop = start + 8 * 1024 * 1024
    return (start, max(stop, start + pad))

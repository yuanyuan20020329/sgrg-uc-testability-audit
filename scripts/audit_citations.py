# -*- coding: utf-8 -*-
"""Citation audit: every reference entry must be cited, every citation must resolve."""
import re

P = "C:/Users/H-P/WorkBuddy/应激颗粒/results/主论文初稿_v1_EN.md"
txt = open(P, encoding="utf-8").read()

# split body vs reference list
marker = "## 7. References"
body, refs = txt.split(marker, 1)

# reference entries:  ^N. ...
entries = {}
for m in re.finditer(r"^(\d+)\.\s+(.*)$", refs, flags=re.M):
    entries[int(m.group(1))] = m.group(2)[:110]

# citations in body:  [1] [2,3] [4,5,6] [1-3] [12–17]  (en/em dash ranges included)
cites = set()
for m in re.finditer(r"\[([\d,\s\u2010-\u2015]+)\]", body):
    for part in re.split(r"[,]", m.group(1)):
        part = part.strip()
        if re.search(r"[\u2010-\u2015-]", part):
            a, b = re.split(r"[\u2010-\u2015-]", part, maxsplit=1)
            try:
                for n in range(int(a), int(b) + 1):
                    cites.add(n)
            except ValueError:
                pass
        elif part.isdigit():
            cites.add(int(part))

print("entries defined : %d  (%s)" % (len(entries), sorted(entries)))
print("citations used  : %d  (%s)" % (len(cites), sorted(cites)))
print()
missing_entry = sorted(c for c in cites if c not in entries)
orphan_entry = sorted(e for e in entries if e not in cites)
print("cited but NOT defined  :", missing_entry)
print("defined but NOT cited  :", orphan_entry)
print()
if orphan_entry:
    for e in orphan_entry:
        print("  [%d] %s" % (e, entries[e]))
print()
print("--- entries still containing placeholders/brackets ---")
for n, t in sorted(entries.items()):
    if "〔" in t or "待补" in t or "⚠" in t or "重复" in t or "确认" in t:
        print("  [%d] %s" % (n, t))

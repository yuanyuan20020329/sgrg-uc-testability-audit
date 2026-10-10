# -*- coding: utf-8 -*-
"""检查 CELLxGENE h5ad 结构"""
import sys
import numpy as np
import h5py

path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/sc/CxG_IBD_HC_10x.h5ad"
print("FILE:", path)
f = h5py.File(path, "r")


def walk(name, obj):
    if isinstance(obj, h5py.Dataset):
        print("  D %-55s %-14s %s" % (name, obj.shape, obj.dtype))
    else:
        print("  G %-55s attrs=%s" % (name, dict(obj.attrs)))


print("--- 结构 ---")
f.visititems(walk)

print("\n--- 关键 ---")
print("X attrs:", dict(f["X"].attrs) if "X" in f else None)
if isinstance(f.get("X"), h5py.Group):
    ks = list(f["X"].keys())
    print("X group:", ks, "dtype=", getattr(f["X"][ks[0]] if ks else None, "dtype", None))
    for k in ks:
        print("   ", k, f["X"][k].shape, f["X"][k].dtype)
print("obs cols:", [k for k in f["obs"].keys()])
print("var cols:", [k for k in f["var"].keys()])
print("var index name attrs:", dict(f["var"].attrs), dict(f["var"]["_index"].attrs) if "_index" in f["var"] else "")
print("obs index name attrs:", dict(f["obs"].attrs))

# 主键取值
def cats(g, name, n=25):
    if name not in g:
        return
    d = g[name]
    if isinstance(d, h5py.Group):
        c = d["categories"][:]
        print("  %s: %d 类 ->" % (name, len(c)), [x.decode() if isinstance(x, bytes) else x for x in c[:n]])
    else:
        print("  %s: raw" % name)

for k in ["cell_type", "cell_type_ontology_term_id", "disease", "disease_ontology_term_id",
          "tissue", "tissue_ontology_term_id", "assay", "donor_id", "sample_id",
          "suspension_type", "self_reported_ethnicity", "development_stage", "sex", "is_primary_data"]:
    cats(f["obs"], k)

if "raw" in f:
    print("\n--- raw ---")
    f["raw"].visititems(lambda n, o: print("  %-50s %s" % (n, getattr(o, "shape", ""))))
f.close()

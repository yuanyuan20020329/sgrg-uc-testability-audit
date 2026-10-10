# -*- coding: utf-8 -*-
"""Part B 辅助 —— 统计各 MDD 队列/子集的样本构成（疾病 × 脑区）"""
import json, gzip, io, time, urllib.request, urllib.parse

def get(url, tries=5):
    last=None
    for i in range(tries):
        try:
            req=urllib.request.Request(url, headers={"Accept-Encoding":"gzip","User-Agent":"curl/8.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                raw=r.read()
                if r.headers.get("Content-Encoding")=="gzip":
                    raw=gzip.GzipFile(fileobj=io.BytesIO(raw)).read()
                return json.loads(raw)
        except Exception as e:
            last=e; time.sleep(3+4*i)
    raise RuntimeError(f"{url} :: {last}")

def counts(dsid, cats_disease, cats_region):
    cursor=""; res={}
    while True:
        url=f"https://gemma.msl.ubc.ca/rest/v2/datasets/{dsid}/samples?limit=50"
        if cursor: url+=f"&cursor={urllib.parse.quote(cursor)}"
        d=get(url)
        for s in d.get("data",[]):
            ch={c.get("category"):c.get("value") for c in (s.get("sample",{}).get("characteristics") or [])}
            dis=next((ch[c] for c in cats_disease if c in ch), None)
            reg=next((ch[c] for c in cats_region if c in ch), None)
            key=(dis, reg)
            res[key]=res.get(key,0)+1
        nc=d.get("nextCursor")
        if not nc: break
        cursor=nc
    return res

out={}
def sk(d):
    return {f"{k[0]} || {k[1]}": v for k, v in d.items()}
out["GSE19738"]=sk(counts(3614, ["disease"], []))
out["GSE80655"]=sk(counts(13014, ["clinical diagnosis"], ["brain region","BioSource","organism part"]))
out["GSE53987"]=sk(counts(8359, ["disease"], ["organism part"]))
print(json.dumps(out, ensure_ascii=False, indent=2))
with open(r"C:/Users/H-P/WorkBuddy/应激颗粒/results/partB_cohort_counts.json","w",encoding="utf-8") as f:
    json.dump(out,f,ensure_ascii=False,indent=2)

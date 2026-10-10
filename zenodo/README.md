# Zenodo deposit

This directory holds the metadata and checksums for the archived release. The
release itself is the whole repository snapshot; Zenodo mints the DOI.

## Files

| File | Purpose |
|---|---|
| `manifest.json` | SHA-256 for every shipped result table, plus the payload file list |
| `README.md` | this file |
| `deposit_metadata.md` | text ready to paste into the Zenodo web form |

## How to mint the DOI

1. In the GitHub repository, create a release, e.g. tag `v1.0.0`.
2. Log in to Zenodo → *Settings* → *GitHub* → toggle the repository **On**.
3. Create the GitHub release; Zenodo archives the tag automatically and mints a
   DOI. (Alternatively upload a ZIP of this directory manually.)
4. Paste the minted DOI into
   * `CITATION.cff` → `identifiers:`
   * the manuscript's *Data and code availability* section, replacing the
     repository-URL placeholder.

## Scope of the deposit

* **Included:** all analysis code, all key result tables, the manuscript
  sources, and the pinned environments.
* **Not included:** any primary public data (`data/` is excluded). Every primary
  dataset is identified by accession in `DATA_ACCESS.md`, together with the
  script that retrieves it.

## Integrity

`manifest.json` lists a SHA-256 for each file under `results/`. To verify a
downloaded copy:

```bash
python - <<'PY'
import hashlib, json, pathlib
man = json.load(open("zenodo/manifest.json", encoding="utf-8"))
bad = []
for name, want in man["sha256"].items():
    p = pathlib.Path("results") / name
    got = hashlib.sha256(p.read_bytes()).hexdigest()
    if got != want:
        bad.append(name)
print("checked", len(man["sha256"]), "files;", "mismatches:", bad or "none")
PY
```

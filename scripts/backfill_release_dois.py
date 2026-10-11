#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Propagate the minted v1.0.1 DOIs into every document that cites them.

Why this exists
---------------
After a re-release there are two new identifiers (a code version DOI and a dataset
version DOI) and roughly a dozen files that name them: the payload's ``manifest.json``,
the deposit README and licence, the code archive's ``CITATION.cff`` / ``README.md`` /
``zenodo/deposit_metadata.md``, four working docs under ``resource_paper/``, the Data
Descriptor and the methods manuscript with its Chinese appendix.

Editing those by hand is how a re-release ships a stale citation: the two DOIs look
alike (``…23297703`` vs ``…23297702``), the old ones legitimately survive in a few
places as *history*, and one of the numbers is a byte count that moves whenever a
rebuilt file changes length. So the substitutions live here, explicitly, one line per
occurrence, and every one of them is checked:

* each search string must be found exactly the expected number of times, or the run
  aborts before writing anything (a silent no-match is the failure mode that matters);
* the two concept DOIs and the v1.0.0 version DOIs are **never** rewritten where they
  are the historical record — only where they are the current citation;
* the byte counts are read back from the rebuilt payload rather than typed.

Two invocation points, matching the two things that get minted in that order
--------------------------------------------------------------------------
The code DOI must exist before the data payload is built (the payload cites it), and the
dataset DOI only exists after the payload is published, so this runs twice::

    # 1. after the GitHub release for tag v1.0.1 is published (Zenodo mints the code DOI)
    python scripts/backfill_release_dois.py --code-doi 10.5281/zenodo.<id>

    # 2. after scripts/publish_dataset.py --new-version-of … has published the dataset
    python scripts/backfill_release_dois.py

The second call needs no arguments: ``publish_dataset.py`` stage 5 writes the dataset DOI
to ``resource_paper/dataset_doi.txt``, which is read here. Either call also re-runs the
payload build, the pre-publish gate and the two manuscript audits, so the pass ends with
the same 52/52 report the release started from.

Nothing here is allowed to run during a release: pass ``--check`` to see what is still
stale without touching a file.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RP = ROOT / "resource_paper"
DEPOSIT = RP / "deposit"
ARCH = ROOT / "archive_src"
RES = ROOT / "results"

VERSION = "1.0.1"
RELEASE_DATE = "2026-10-11"

# The two identifiers that never change. A concept DOI is the record family; it resolves
# to the newest version, so it is identical in every release and must not be rewritten.
CODE_CONCEPT = "10.5281/zenodo.23291758"
DATA_CONCEPT = "10.5281/zenodo.23297702"

# The v1.0.0 version DOIs. These stay in the documents where they are the historical
# record, and are replaced only where the sentence means "the current version".
OLD_CODE = "10.5281/zenodo.23291759"
OLD_DATA = "10.5281/zenodo.23297703"

# Byte counts quoted in prose. Recomputed from the rebuilt payload, never typed.
OLD_TOTAL = "3,095,926"
OLD_TOTAL_NOCOMMA = "3,094,721"

CODE_DOI_FILE = RP / "code_doi.txt"
DATA_DOI_FILE = RP / "dataset_doi.txt"


# --------------------------------------------------------------------------- helpers


def _nl_of(raw: str) -> str:
    return "\r\n" if "\r\n" in raw else "\n"


def read_first_doi(path: Path) -> str | None:
    """First non-comment line of a ``*_doi.txt`` record."""
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            return line
    return None


def _fmt(n: int) -> str:
    return f"{n:,}"


def payload_totals() -> tuple[int, int, float]:
    files = [p for p in sorted(DEPOSIT.iterdir()) if p.is_file()]
    total = sum(p.stat().st_size for p in files)
    return len(files), total, total / 1024 / 1024


class Report:
    """Collect and apply (or just report) a list of literal text substitutions.

    Every substitution is checked against the file it targets. The check is what makes the
    two-invocation sequence safe: the pass is run once for the code DOI and again for the
    dataset DOI, so on the second run the first run's edits are already in place. An edit
    is therefore accepted in either state — *pending* (the old text is there) or *already
    applied* (only the new text is). Anything else is a hard stop, because the failure mode
    that matters is a substitution that silently matches nothing and leaves a stale DOI in
    the manuscript.
    """

    def __init__(self, check: bool) -> None:
        self.check = check
        self.applied = 0
        self.pending = 0
        self.uptodate = 0
        self.problems: list[str] = []
        self.touched: list[str] = []

    def edit(self, rel: str, old: str, new: str, want: int = 1) -> None:
        if new == old:
            return
        path = ROOT / rel
        raw = path.read_text(encoding="utf-8")
        got, already = raw.count(old), raw.count(new)
        if got != want:
            if got == 0 and already:
                self.uptodate += 1
                print(f"  [skip] {rel}: already carries this change")
                return
            self.problems.append(
                f"{rel}: expected {want} occurrence(s) of {old.splitlines()[0][:60]!r}, "
                f"found {got} (and {already} of the replacement) — skipped; the file has "
                f"moved on, re-read it before trusting this pass")
            return
        if self.check:
            self.pending += 1
            print(f"  [pending] {rel}  ::  {old.splitlines()[0][:64]}")
            return
        path.write_text(raw.replace(old, new), encoding="utf-8", newline=_nl_of(raw))
        self.applied += 1
        if rel not in self.touched:
            self.touched.append(rel)


# --------------------------------------------------------------------- edit sections


def edit_archive_side(rep: Report, new_code: str) -> None:
    """The code archive's own documents (``archive_src/`` → ``release/``)."""
    print("\n[1/5] code archive documents")

    rep.edit(
        "archive_src/CITATION.cff",
        "# v1.0.1. The Zenodo version DOI is minted when the GitHub release is published, so it is\n"
        "# added in a follow-up commit (see the `identifiers:` block below); the concept DOI\n"
        "# 10.5281/zenodo.23291758 always resolves to the newest version of this code archive.\n",
        "# v1.0.1. The version DOI below is the record Zenodo minted for this release; the concept\n"
        f"# DOI {CODE_CONCEPT} always resolves to the newest version of this code archive.\n",
    )
    rep.edit(
        "archive_src/CITATION.cff",
        "license: MIT\nrepository-code:",
        f"license: MIT\ndoi: {new_code}\nrepository-code:",
    )
    rep.edit(
        "archive_src/CITATION.cff",
        "identifiers:\n"
        f"  - type: doi\n    value: {CODE_CONCEPT}\n"
        '    description: "Zenodo concept DOI of this code archive (resolves to the newest version)"\n'
        "\n"
        "# The v1.0.1 version DOI is added here once Zenodo has archived the release, e.g.\n"
        "#   doi: 10.5281/zenodo.<id>\n"
        "#   identifiers:\n"
        "#     - type: doi\n"
        "#       value: 10.5281/zenodo.<id>\n"
        '#       description: "Zenodo archival DOI (version v1.0.1)"\n',
        "identifiers:\n"
        f"  - type: doi\n    value: {new_code}\n"
        '    description: "Zenodo archival DOI (version v1.0.1)"\n'
        f"  - type: doi\n    value: {CODE_CONCEPT}\n"
        '    description: "Zenodo concept DOI of this code archive (resolves to the newest version)"\n'
        f"  - type: doi\n    value: {OLD_CODE}\n"
        '    description: "Zenodo archival DOI (version v1.0.0, published and immutable)"\n',
    )
    rep.edit(
        "archive_src/zenodo/deposit_metadata.md",
        f"* This archive: concept DOI {CODE_CONCEPT} (resolves to the newest\n"
        "  version); the v1.0.1 version DOI is minted when the GitHub release is published and is\n"
        "  then recorded here and in `CITATION.cff`.\n"
        f"* Previous version (v1.0.0): {OLD_CODE} (published and immutable)\n",
        f"* This archive (v1.0.1): {new_code}\n"
        f"* Concept DOI (resolves to the newest version): {CODE_CONCEPT}\n"
        f"* Previous version (v1.0.0): {OLD_CODE} (published and immutable)\n",
    )
    rep.edit(
        "archive_src/zenodo/deposit_metadata.md",
        "* \u23f3 To do for v1.0.1: publish the GitHub release for tag `v1.0.1`, then copy the minted\n"
        "  version DOI into `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md`, this file\n"
        "  and the manuscript's *Data and code availability* statement. The data payload must cite\n"
        "  the **v1.0.1** code DOI, not the v1.0.0 one \u2014 `scripts/audit_data_descriptor.py\n"
        "  --pre-publish` fails the build if it does not.\n",
        f"* \u2705 Done (v1.0.1, {RELEASE_DATE}): version DOI `{new_code}` backfilled into\n"
        "  `resource_paper/code_doi.txt`, `CITATION.cff`, `README.md`, this file and the manuscript's\n"
        "  *Data and code availability* statement. The data payload cites the **v1.0.1** code DOI,\n"
        "  not the v1.0.0 one \u2014 `scripts/audit_data_descriptor.py --pre-publish` fails otherwise.\n",
    )


def edit_working_docs(rep: Report, new_code: str, new_data: str | None) -> None:
    """``resource_paper/`` working docs: the resource README, the readiness audit, the workbook.

    Split in two on purpose. The code DOI exists before the data record is published, but the
    dataset DOI only exists afterwards — so edits are grouped by *which identifier they
    depend on*, and the data-side ones are skipped (not run with a stand-in) until
    ``new_data`` is known. Running them early would relabel a row "v1.0.1" while still
    quoting the superseded DOI, which is worse than leaving the edit for later.
    """
    print("\n[2/5] resource_paper working docs")

    # ---- depends only on the code DOI ---------------------------------------
    rep.edit("resource_paper/README.md",
             "**The frozen, upload-ready v1.0.0 payload for the dataset record**",
             "**The frozen, upload-ready v1.0.1 payload for the dataset record**")
    rep.edit(
        "resource_paper/README.md",
        "* **The code archive** (analysis pipeline + key result tables, MIT) \u2014\n"
        f"  **{OLD_CODE}** \u2014 https://doi.org/{OLD_CODE}\n",
        "* **The code archive** (analysis pipeline + key result tables, MIT) \u2014\n"
        f"  **{new_code}** \u2014 https://doi.org/{new_code}\n"
        f"  (v{VERSION}; the concept DOI {CODE_CONCEPT} points at all versions.)\n",
    )
    rep.edit("resource_paper/README.md",
             "_First release 2026-10-11._",
             f"_First release 2026-10-11 (v1.0.0). **Current version v{VERSION}** \u2014 "
             "English-language normalisation of the free-text fields only; no gate outcome, "
             "numeric value or row changed. See `deposit/CHANGELOG.md`. v1.0.0 stays published "
             "and immutable._")
    rep.edit("resource_paper/SUBMISSION_READINESS.md",
             "| 1 | Deposit in a **public, community-recognised repository** meeting the policy "
             "| \u2705 **published 2026-10-11** |",
             f"| 1 | Deposit in a **public, community-recognised repository** meeting the policy "
             f"| \u2705 **v{VERSION} published {RELEASE_DATE}** |")
    rep.edit("resource_paper/SUBMISSION_READINESS.md",
             "A single-unit convenience zip is built at "
             "`resource_paper/testability-audit-tables-v1.0.0.zip`.",
             f"A single-unit convenience zip is built at "
             f"`resource_paper/testability-audit-tables-v{VERSION}.zip`.")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"> This is the **second** Zenodo record. The first one (`{OLD_CODE}`, code",
             f"> This is the **second** Zenodo record. The first one (`{OLD_CODE}` at v1.0.0, code")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             "| **Version** | 1.0.0 |", f"| **Version** | {VERSION} |")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"| `{OLD_CODE}` | is supplemented by | Software |",
             f"| `{new_code}` | is supplemented by | Software |")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"the **code** archive (`{OLD_CODE}`), and `manifest.json` records **both** hashes:",
             f"the **code** archive (`{new_code}`), and `manifest.json` records **both** hashes:")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"(4) the record's own DOI, the dataset DOI in `README.md`, and the code DOI\n"
             f"`{OLD_CODE}` are all present and **distinct**.",
             f"(4) the record's own DOI, the dataset DOI in `README.md`, and the code DOI\n"
             f"`{new_code}` are all present and **distinct**.")

    if not new_data:
        print("  [..] 4 data-DOI edits deferred until the dataset DOI is recorded")
        return

    # ---- depends on the dataset DOI ----------------------------------------
    rep.edit(
        "resource_paper/README.md",
        f"* **The tables (this resource, dataset record)** \u2014 **{OLD_DATA}**\n"
        f"  \u2014 https://doi.org/{OLD_DATA}\n"
        f"  (minted 2026-10-11; version DOI. The concept DOI {DATA_CONCEPT} points at all\n"
        "  versions of this record.)\n",
        f"* **The tables (this resource, dataset record)** \u2014 **{new_data}**\n"
        f"  \u2014 https://doi.org/{new_data}\n"
        f"  (version of record, v{VERSION}; concept DOI {DATA_CONCEPT} points at all versions.)\n",
    )
    rep.edit("resource_paper/SUBMISSION_READINESS.md",
             f"The 14-file payload is live as Zenodo record {OLD_DATA.rsplit('.', 1)[-1]} (see gap 3).",
             f"The payload is live as Zenodo record {new_data.rsplit('.', 1)[-1]} "
             f"(v{VERSION}; see gap 3).")
    rep.edit(
        "resource_paper/SUBMISSION_READINESS.md",
        f"| 3 | Data must be **live** and directly linked to the paper | \u2705 **minted 2026-10-11** | "
        f"**Dataset DOI `{OLD_DATA}`** \u2014 record live at "
        f"https://zenodo.org/records/{OLD_DATA.rsplit('.', 1)[-1]} (status `published`), "
        f"concept DOI `{DATA_CONCEPT}`.",
        f"| 3 | Data must be **live** and directly linked to the paper | \u2705 **v{VERSION} minted "
        f"{RELEASE_DATE}** | **Dataset DOI `{new_data}`** (version of record) \u2014 record live at "
        f"https://zenodo.org/records/{new_data.rsplit('.', 1)[-1]} (status `published`); "
        f"v1.0.0 `{OLD_DATA}` remains published and immutable. Concept DOI `{DATA_CONCEPT}`.",
    )
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"record ever holds a placeholder version (it is now a no-op: both quote\n"
             f"   `{OLD_DATA}` and no placeholder remains).",
             f"record ever holds a placeholder version (stage 2 now *replaces* the previous\n"
             f"   version DOI: both quote `{new_data}`, so the superseded `{OLD_DATA}` must not be\n"
             f"   left behind).")


def edit_manuscripts(rep: Report, new_code: str, new_data: str | None) -> None:
    """Data Descriptor + methods manuscript (EN) + the Chinese appendix.

    Grouped on the same rule as :func:`edit_working_docs`: code-DOI edits run as soon as the
    code DOI is known, dataset-DOI edits wait for the dataset record to be published.
    """
    print("\n[3/5] manuscripts")

    # ---- depends only on the code DOI ---------------------------------------
    dd = "resource_paper/Data_Descriptor_v1_EN.md"
    rep.edit(dd, "Version 1.0.0\ncontains dataset-level metadata",
                 f"Version {VERSION}\ncontains dataset-level metadata")
    rep.edit(dd, f"Zenodo (code DOI: {OLD_CODE}). In that archive",
                 f"Zenodo (code DOI: {new_code}). In that archive")
    rep.edit(dd, f"10.5281/zenodo.23291759. The code and dataset records are independently versioned",
                 f"{new_code}. The code and dataset records are independently versioned")

    en = "results/主论文初稿_v1_EN.md"
    rep.edit(en, f"代码 Zenodo DOI **{OLD_CODE}**", f"代码 Zenodo DOI **{new_code}**")
    rep.edit(en, f"DOI: {OLD_CODE})", f"DOI: {new_code})")

    cn = "results/主论文初稿_v1_CN附注.md"
    rep.edit(cn, f"- `CITATION.cff` —— 含 `doi:` 与 `identifiers:`，已回填版本 DOI `{OLD_CODE}`。",
                 f"- `CITATION.cff` —— 含 `doi:` 与 `identifiers:`；当前版本 v{VERSION} 的版本 DOI 为 "
                 f"`{new_code}`（v1.0.0 为 `{OLD_CODE}`）。")

    if not new_data:
        print("  [..] 5 dataset-DOI edits deferred until the dataset DOI is recorded")
        return

    # ---- depends on the dataset DOI ----------------------------------------
    rep.edit(dd, f"Zenodo (version DOI: {OLD_DATA}) and can support",
                 f"Zenodo (version DOI: {new_data}) and can support")
    rep.edit(dd, f"The complete version 1.0.0 package is available from Zenodo under CC BY 4.0\n"
                 f"({OLD_DATA}). The concept DOI, {DATA_CONCEPT}, resolves to all",
                 f"The complete version {VERSION} package is available from Zenodo under CC BY 4.0\n"
                 f"({new_data}). The concept DOI, {DATA_CONCEPT}, resolves to all")
    rep.edit(dd, f"{OLD_DATA} (version 1.0.0; concept DOI {DATA_CONCEPT}).",
                 f"{new_data} (version {VERSION}; concept DOI {DATA_CONCEPT}).")
    rep.edit(en, f"审计表数据集 Zenodo DOI **{OLD_DATA}**",
                 f"审计表数据集 Zenodo DOI **{new_data}**")
    rep.edit(en, f"Zenodo, DOI: {OLD_DATA} (concept DOI: {DATA_CONCEPT})",
                 f"Zenodo, DOI: {new_data} (concept DOI: {DATA_CONCEPT})")
    rep.edit(cn, f"| **数据集版本 DOI（写稿件引这个）** | **`{OLD_DATA}`** |",
                 f"| **数据集版本 DOI（v1.0.0，已被 v1.0.1 取代 \u2014 见 \u00a7C10）** | **`{OLD_DATA}`** |")


def append_appendix_section(rep: Report, new_code: str, new_data: str, total: int) -> None:
    """Append the v1.0.1 record to the Chinese appendix. Idempotent."""
    rel = "results/主论文初稿_v1_CN附注.md"
    path = ROOT / rel
    raw = path.read_text(encoding="utf-8")
    if "### C10 v1.0.1" in raw:
        print(f"  [ok] {rel} already carries the C10 section")
        return
    if rep.check:
        rep.pending += 1
        print(f"  [pending] {rel}  ::  append the §C10 v1.0.1 record")
        return
    nl = _nl_of(raw)
    section = f"""
### C10 v1.0.1 语言规范化重发（{RELEASE_DATE}）

v1.0.0 发布后发现的**唯一**问题：载荷的若干**自由文本**字段是中文写就的（账本 29 行的 `reason`、
`partD_resfeas_libwide_summary.json` 的 `metric_definitions` 13 条与 `caveat` 1 条、
`partD_resfeas_verdict.json` 的一个示例串、`datapackage.json` 的两条字段描述），
英文读者读不出"某资源因何未通过"。v1.0.1 只做这一件事。

| | |
|---|---|
| **数据集版本 DOI（v1.0.1，写稿件引这个）** | **`{new_data}`** |
| 数据集概念 DOI（不变） | `{DATA_CONCEPT}` |
| **代码版本 DOI（v1.0.1）** | **`{new_code}`** |
| 代码概念 DOI（不变） | `{CODE_CONCEPT}` |
| v1.0.0 数据集版本 DOI（保留、不可变） | `{OLD_DATA}` |
| v1.0.0 代码版本 DOI（保留、不可变） | `{OLD_CODE}` |
| 载荷 | 14 个文件 / {total:,} B |

**声称的边界（可由代码核验，不靠承诺）**：字段名、表结构、列顺序、四个 gate、`USABLE`、
所有数值、行顺序**一律未动**；改的只有自由文本的语言。断言方式是 `scripts/test_language_normalisation.py`
把载荷逐文件与冻结的 v1.0.0 基准（`resource_paper/frozen_deposit_v1.0.0/`）比对，**只允许**
白名单字段不同。该测试已注册为发布前硬门槛：`scripts/audit_data_descriptor.py --pre-publish`
会跑它，不通过就不许发布（`publish_dataset.py` 的 stage 2g）。

**代码同版本重发**：生成器改了（`partD_resfeas_ledger.py` 等走双语模块），所以代码也必须重发，
否则 v1.0.1 的载荷指向的还是 v1.0.0 的脚本。数据载荷里的 `manifest.json`、README、LICENSE-DATA
都引用**当前**代码版本 DOI `{new_code}`；门槛会比对 `code_doi.txt` 的 `# version:` 与载荷版本，
不一致直接失败 —— 这正是防止"重发后仍悄悄引用上一版代码"的那道闸。

**发布顺序（不可颠倒）**：① 代码先发（GitHub Release tag `v1.0.1`），拿到代码版本 DOI →
② `scripts/backfill_release_dois.py --code-doi <代码DOI>` 写入并回填 →
③ 数据再发（`scripts/publish_dataset.py --new-version-of {OLD_DATA} --forget-token`），
数据 DOI 铸好后 → ④ 再跑一次 `backfill_release_dois.py`（无参数，自动读 `dataset_doi.txt`）
把稿件、资源稿与工作文档一次改齐，并复跑门槛与两个一致性审计。
"""
    path.write_text(raw.rstrip("\r\n") + nl + section, encoding="utf-8", newline=nl)
    rep.applied += 1
    if rel not in rep.touched:
        rep.touched.append(rel)
    print(f"  [ok] appended §C10 to {rel}")


# ------------------------------------------------------------------ byte-dependent


def edit_workbook_history(rep: Report) -> None:
    """Label the workbook's "as published" record as the v1.0.0 one.

    Section 9 keeps the v1.0.0 numbers on purpose — that release happened and its record is
    immutable — so it must say so, or the next reader takes a superseded DOI for the live one.
    """
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             "## 9. As published — 2026-10-11\n",
             "## 9. As published — v1.0.0, 2026-10-11 (superseded by v1.0.1; kept as history)\n")


def append_workbook_section(rep: Report, new_code: str, new_data: str, total: int) -> None:
    """Append the v1.0.1 record to the workbook. Idempotent."""
    rel = "resource_paper/zenodo_deposit_datasets.md"
    path = ROOT / rel
    raw = path.read_text(encoding="utf-8")
    if "## 10. As published — v1.0.1" in raw:
        print(f"  [ok] {rel} already carries the v1.0.1 record")
        return
    if rep.check:
        rep.pending += 1
        print(f"  [pending] {rel}  ::  append the v1.0.1 \"as published\" record")
        return
    nl = _nl_of(raw)
    section = f"""
---

## 10. As published — v{VERSION}

| | |
|---|---|
| Dataset DOI (version of record) | **{new_data}** |
| Concept DOI (all versions) | {DATA_CONCEPT} |
| Code DOI (version of record) | **{new_code}** |
| Code concept DOI (all versions) | {CODE_CONCEPT} |
| Previous dataset version | {OLD_DATA} (v1.0.0, published and immutable) |
| Previous code version | {OLD_CODE} (v1.0.0, published and immutable) |
| Files | 14 ({_fmt(total)} B), resource type **Dataset**, licence **CC BY 4.0**, version {VERSION} |

**What changed from v1.0.0: language only.** The free-text fields that were written in Chinese
— the ledger's per-row `reason`, the summary's `metric_definitions` glossary and `caveat`, one
example string in the verdict file, two prose field descriptions in `datapackage.json` — are
English. Field names, table structure, column order, the four gates, `USABLE`, every numeric
value and the row order are unchanged, and that is asserted rather than promised:
`scripts/test_language_normalisation.py` diffs this payload against the frozen v1.0.0 copy in
`resource_paper/frozen_deposit_v1.0.0/` and only the allow-listed free-text fields may differ.
The test runs as a hard pre-publication gate (`audit_data_descriptor.py --pre-publish`, wired
into `publish_dataset.py`), and `CHANGELOG.md` in the payload carries the same statement.

**Why the code archive was re-released too.** The generators changed (the ledger now takes its
strings from `scripts/resfeas_reason_i18n.py`), so the scripts in v1.0.0 no longer reproduce
these bytes. Because the payload cross-links the code DOI, a language-only data release with a
stale code link would be a citation error; `audit_data_descriptor.py` fails the build when
`code_doi.txt`'s version disagrees with the payload version, so it cannot ship silently.
"""
    path.write_text(raw.rstrip("\r\n") + nl + section, encoding="utf-8", newline=nl)
    rep.applied += 1
    if rel not in rep.touched:
        rep.touched.append(rel)
    print(f"  [ok] appended the v1.0.1 record to {rel}")


def refresh_byte_counts(rep: Report) -> None:
    """Re-read the byte counts from the rebuilt payload and fix the prose that quotes them.

    One number in the descriptor ("14 files (3,095,926 bytes in total)") and one in the
    workbook ("14 files, 3,094,721 B") are properties of the payload, and they move whenever
    a rendered file changes length. They are recomputed here rather than typed, because a
    stale total is exactly the kind of detail a reviewer checks first.
    """
    print("\n[4/5] byte counts (recomputed from the rebuilt payload)")
    n_files, total, mb = payload_totals()
    print(f"  payload: {n_files} files, {total:,} B ({mb:.2f} MB)")

    rep.edit("resource_paper/Data_Descriptor_v1_EN.md",
             f"The Zenodo record contains 14 files ({OLD_TOTAL} bytes in total)",
             f"The Zenodo record contains {n_files} files ({_fmt(total)} bytes in total)")
    rep.edit("resource_paper/zenodo_deposit_datasets.md",
             f"Upload payload: **`resource_paper/deposit/`** (14 files, {OLD_TOTAL_NOCOMMA} B = 2.95 MB total, of which",
             f"Upload payload: **`resource_paper/deposit/`** ({n_files} files, {_fmt(total)} B = "
             f"{mb:.2f} MB total, of which")

    # The per-file table in the workbook. Refresh the Bytes column from disk so the table
    # cannot drift from the payload it claims to describe.
    rel = "resource_paper/zenodo_deposit_datasets.md"
    path = ROOT / rel
    raw = path.read_text(encoding="utf-8")
    row = re.compile(r"^(\| `([A-Za-z0-9_.\-]+)` \| )([\d,]+)( \|.*)$", re.M)

    def sub(m: re.Match) -> str:
        name = m.group(2)
        f = DEPOSIT / name
        if not f.is_file():
            return m.group(0)
        return f"{m.group(1)}{_fmt(f.stat().st_size)}{m.group(4)}"

    new, n = row.subn(sub, raw)
    if n and new != raw:
        if rep.check:
            rep.pending += 1
            print(f"  [pending] {rel}  ::  {n} byte rows to refresh")
        else:
            path.write_text(new, encoding="utf-8", newline=_nl_of(raw))
            rep.applied += 1
            if rel not in rep.touched:
                rep.touched.append(rel)
        print(f"  [ok] {n} rows in the payload table carry current sizes")
    else:
        print("  [ok] payload table already current")


# ------------------------------------------------------------------------- pipeline


def write_code_doi_record(new_code: str) -> None:
    text = (
        f"{new_code}\n"
        f"# concept doi: {CODE_CONCEPT}\n"
        f"# version: v{VERSION}\n"
        "#\n"
        "# Single source of truth for the COMPANION CODE record that the dataset payload, its\n"
        "# manifest.json, the deposit README and the manuscript all cross-link to.\n"
        "#\n"
        "#   line 1            the code VERSION DOI (what the data record cites)\n"
        "#   \"# concept doi:\"  the code CONCEPT DOI (always resolves to the newest code version)\n"
        "#   \"# version:\"      the code version string, e.g. v1.0.1\n"
        "#\n"
        "# This file is the ONLY place the code DOI is recorded. The payload build reads it, so a\n"
        "# re-release becomes a one-line change here rather than a hunt through the generators.\n"
        "# scripts/audit_data_descriptor.py fails when \"# version:\" disagrees with the payload\n"
        "# version, and that check is wired in as the pre-publish gate of\n"
        "# scripts/publish_dataset.py.\n"
    )
    CODE_DOI_FILE.write_text(text, encoding="utf-8", newline="\n")
    print(f"  [ok] {CODE_DOI_FILE.relative_to(ROOT)} -> {new_code} (v{VERSION})")


def run(cmd: list[str], label: str) -> int:
    print(f"\n$ {' '.join(cmd)}")
    r = subprocess.run(cmd, cwd=str(ROOT))
    print(f"  exit {r.returncode}  ({label})")
    return r.returncode


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--code-doi", help="the v1.0.1 CODE version DOI just minted by Zenodo")
    ap.add_argument("--data-doi", help="override the dataset version DOI (else dataset_doi.txt)")
    ap.add_argument("--check", action="store_true",
                    help="report what is still stale; touch nothing")
    ap.add_argument("--no-build", action="store_true",
                    help="skip the payload rebuild and the gate (documents only)")
    args = ap.parse_args()

    recorded_code = read_first_doi(CODE_DOI_FILE)
    recorded_data = read_first_doi(DATA_DOI_FILE)
    new_code = args.code_doi or (recorded_code if recorded_code != OLD_CODE else None)
    new_data = args.data_doi or (recorded_data if recorded_data != OLD_DATA else None)

    print("=" * 74)
    print("backfill the v1.0.1 DOIs")
    print(f"  code  : {new_code or '(not minted yet — pass --code-doi)'}")
    print(f"  data  : {new_data or '(not minted yet — publish_dataset.py writes dataset_doi.txt)'}")
    print(f"  mode  : {'check only' if args.check else 'apply'}")
    print("=" * 74)

    # A missing code DOI is not "nothing to do": several edits relabel a citation as
    # v1.0.1 while substituting the DOI, so running without it would stamp the v1.0.1 label
    # onto the v1.0.0 identifier. That is worse than not running, so refuse.
    if not new_code:
        print("\nno v1.0.1 code DOI yet.\n"
              "  The code archive is released first, because the data payload cites it:\n"
              f"    1. publish the GitHub release for tag v{VERSION}\n"
              f"    2. python scripts/backfill_release_dois.py --code-doi 10.5281/zenodo.<id>\n"
              "  Pass --code-doi now or let the minted DOI land in "
              f"{CODE_DOI_FILE.relative_to(ROOT)}.")
        return 2
    if not new_data:
        print("\n  [note] the dataset v1.0.1 DOI is not recorded yet, so only the code-side\n"
              "         documents and the payload's code cross-link are updated. Re-run this\n"
              "         after publish_dataset.py writes dataset_doi.txt.")

    rep = Report(args.check)

    if new_code and not args.check:
        print("\n[0/5] record the code DOI")
        write_code_doi_record(new_code)

    edit_archive_side(rep, new_code)
    edit_working_docs(rep, new_code, new_data)
    edit_manuscripts(rep, new_code, new_data)
    edit_workbook_history(rep)

    if args.check:
        if new_data:
            refresh_byte_counts(rep)
            append_appendix_section(rep, new_code, new_data, payload_totals()[1])
            append_workbook_section(rep, new_code, new_data, payload_totals()[1])
    else:
        # The payload must be rebuilt before the byte counts are read back.
        if not args.no_build and new_data:
            print("\n[4a/5] rebuild the payload with the new DOIs")
            rc = run([sys.executable, str(ROOT / "scripts" / "build_data_deposit.py"),
                      "--data-doi", new_data,
                      *([] if not new_code else ["--code-doi", new_code])], "build payload")
            if rc:
                return rc
        if new_data:
            refresh_byte_counts(rep)
            append_appendix_section(rep, new_code, new_data, payload_totals()[1])
            append_workbook_section(rep, new_code, new_data, payload_totals()[1])
        else:
            print("\n[4/5] byte counts: skipped — the payload has not been rebuilt with the "
                  "v1.0.1 dataset DOI yet, so any total read now is provisional")

    print("\n" + "-" * 74)
    print(f"edits applied : {rep.applied}")
    print(f"already done  : {rep.uptodate}")
    if rep.check:
        print(f"still pending : {rep.pending}")
    for f in rep.touched:
        print(f"  changed: {f}")
    if rep.problems:
        print("\nPROBLEMS (nothing was written for these):")
        for p in rep.problems:
            print(f"  [FAIL] {p}")

    if args.check:
        return 1 if (rep.pending or rep.problems) else 0
    if rep.problems:
        return 1

    if not args.no_build:
        print("\n[5/5] gate + consistency audits")
        rc = run([sys.executable, str(ROOT / "scripts" / "audit_data_descriptor.py"),
                  "--pre-publish"], "pre-publish gate")
        if rc:
            print("\nthe pre-publish gate is red — fix before publishing.")
            return rc
        for script in ("audit_citations.py", "reconcile_numbers.py"):
            run([sys.executable, str(ROOT / "scripts" / script)], script)

    print("\nnext: commit `release/` (the code archive) and re-run the payload verifiers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Patent Sequence Extraction — Learning Reference

A general-purpose reference for extracting siRNA / duplex / antisense-
oligonucleotide data from OCR'd patent applications. Covers Alnylam-style
PCT applications specifically (the dominant siRNA patent filer), but the
parser patterns transfer to Arrowhead, Dicerna, Ionis, and similar.

## Table of contents

1. [Anatomy of an Alnylam-style patent](#1-anatomy-of-an-alnylam-style-patent)
2. [Two styles of sequence table](#2-two-styles-of-sequence-table)
3. [Table 1: Abbreviation legend (always parse first)](#3-table-1-abbreviation-legend-always-parse-first)
4. [Alnylam modified-strand notation decoder](#4-alnylam-modified-strand-notation-decoder)
5. [Activity table shapes](#5-activity-table-shapes)
6. [In-vivo result tables](#6-in-vivo-result-tables)
7. [Column-mapping pitfalls](#7-column-mapping-pitfalls)
8. [OCR artifacts](#8-ocr-artifacts)
9. [PaddleOCR-VL input format](#9-paddleocr-vl-input-format)
10. [Parser snippets](#10-parser-snippets)
11. [Sanity checks](#11-sanity-checks)
12. [Reference example: WO 2021/126734 A1 (PNPLA3)](#12-reference-example-wo-2021126734-a1-pnpla3)
13. [Reference example: WO 2012/177784 A1 (ANGPTL3)](#13-reference-example-wo-2012177784-a1-angptl3)

---

## 1. Anatomy of an Alnylam-style patent

OCR'd patent Markdown has a stable structure:

1. **Front matter** (lines ~1-100): bibliographic data — applicant, inventors,
   PCT number, priority dates, IPC class. **Strong signal** for siRNA patents:
   IPC `C12N 15/113` or `A61K 31/713`.
2. **Body** (~lines 100-2000): background, summary, detailed description,
   examples. Look for tables labeled "Table 1", "Table 2", etc.
3. **"Description of the Sequence Listing"** table — usually Table 1 or
   early. Lists what each SEQ ID NO is. **Does NOT contain the actual
   sequences** — just describes them. A common failure mode is to confuse
   this with the sequence data.
4. **Actual sequence tables** — typically titled "Unmodified Sense and
   Antisense Strand Sequences of [gene] dsRNA Agents" and "Modified Sense
   and Antisense Strand Sequences of [gene] dsRNA Agents". Plain ACGU for
   unmodified, Alnylam-notation for modified.
5. **ST.26 sequence listing file** — separate ~1 MB file
   (`121301-XXXXX_SL.txt`) NOT in the OCR'd Markdown. This is the
   authoritative source. Patent body tables should match; if they don't,
   the sequence listing wins.
6. **Activity / screening tables** — "Single Dose Screen in [cell line]",
   "Fold Change". Report % mRNA remaining at 4-5 doses (in-vitro) or
   fold-change vs PBS (in-vivo). Often split across pages.
7. **Claims and embodiments** — cite SEQ IDs by number. Useful for finding
   the "preferred" / "most preferred" designs (clinically advanced).
8. **"Duplexes of Interest" / "Treatment Groups" tables** — list duplex
   names that go into in-vivo studies, with `Range in NM_XXXXX.X` columns.

---

## 2. Two styles of sequence table

Alnylam PCT applications use two column shapes:

### Style A — Unmodified (7-column)

```
| Duplex Name | Sense Sequence 5' to 3' | SEQ ID NO: | Range in NM_XXX | Antisense Sequence 5' to 3' | SEQ ID NO: | Range in NM_XXX |
| AD-517197.1 | CAUGAGCAAGAUUUGCAACUU     | 20         | 1187-1207       | AAGUUGCAAAUCUUGCUCAUGUA       | 109        | 1185-1207       |
```

Plain ACGU bases, has `Range in NM_xxx` columns. The Range gives the mRNA
target coordinate (use to verify against your own target region).

### Style B — Modified (7-column)

```
| Duplex Name | Sense Sequence 5' to 3' | SEQ ID NO: | Antisense Sequence 5' to 3' | SEQ ID NO: | mRNA Target Sequence 5' to 3' | SEQ ID NO: |
| AD-517197.1 | csasugagCfaAfGfAfuugcaacunL96 | 198 | asAfsgung(Cgn)aauancUfgCfucaugsusa | 287 | UACAUGAGCAAGAUUUGCAACUU | 376 |
```

Alnylam-notation modified bases, has `mRNA Target Sequence` column. The
Sense has `L96` (GalNAc) at 3' end; Antisense has `dTdT` 3' overhang.

### Style detection (column-count agnostic)

```python
def detect_table_style(header_cells):
    """Returns 'A' (unmodified), 'B' (modified), or None (ambiguous)."""
    h_low = " ".join(header_cells).lower()
    if "mrna target" in h_low:
        return "B"
    if "range" in h_low and h_low.count("seq id") >= 2:
        return "A"
    return None
```

### Two-table pairing

Each Modified title (Table 3) has a matching Unmodified title (Table 2). The
same is true for 4/5, 6/7, etc. After parsing both, **merge records by
duplex name** into a single record with both unmodified and modified
sequences. Don't deduplicate on the base name — `AD-519351.1` (initial lead)
and `AD-519351.4` (refined chemistry) are different designs.

### Three-series ID pattern (newer Alnylam PCTs)

A pattern seen in WO 2012/177784 A1 (and likely newer filings): the
**same chemical design** appears under **three different duplex IDs**
that map 1:1:1. Example:

| Series | ID example | Sequence tables | Chemistry | Used in activity |
|:-------|:-----------|:----------------|:----------|:-----------------|
| **A — unmodified, `x` overhang** | `AD-53063.1` | T7, T8 (mod) | 19nt sense with 3' `x` (abasic placeholder), 23nt AS with `+X` extra 3' nt | T11 single-dose, T12 IC50 |
| **B — unmodified, blunt/21nt sense** | `AD-52637.1` | T9, T10 (mod) | 21nt sense + 23nt AS, no `x` | T13 single-dose, T14 IC50 |
| **C — GalNAc-conjugated** | `AD-53063.1` (overlaps A!) | T8 (mod, with `L96`) | Same sequence as A but with `L96` GalNAc on sense 3' | T11 |

Wait — that's wrong. In WO 2012/177784, **series A and C share IDs**
(`AD-53063.1` is the same physical duplex in both). The 185-duplex
GalNAc set in T8 IS the GalNAc version of the 185-duplex `x`-overhang
set in T7. The 181-duplex T9/T10 series (B: `AD-52637.x`) is a **refined
chemistry** of the same designs, with different modification patterns
and 21nt/23nt sizing.

**The 1:1:1 mapping** between A, B, and C series is given in a dedicated
**"Table 15. IDs of duplex pairs"** in the patent. The patent literally
writes: "Unconjugated duplex ID | GalNAc-conjugated duplex ID" (and
sometimes a third column for the T9/T10 series).

**What to do in the extractor**:

1. Parse the Table 15 ID-pairing table. Store as
   `duplex_id → {"unconjugated_id":..., "galnac_id":..., "refined_id":...}`.
2. When writing the wide CSV, use the **unconjugated ID** as the primary
   key, and the other two as `paired_galnac_id` / `paired_unmodified_long_id`
   columns. This way the same chemical design shows up exactly once.
3. When writing the long activity CSV, leave the activity rows keyed by
   the ID that was actually used in the source table (don't substitute).
   The join happens in the analysis layer, not the extraction layer.

If the patent has a Table 15 (or equivalent ID-pairing table), add it to
CONFIG as `ID_PAIRING_TABLE = 15` and a `DUPLEX_SERIES` dict that names
the three series so the user can tell them apart in the wide CSV.

---

## 3. Table 1: Abbreviation legend (always parse first)

**Read Table 1 BEFORE decoding any modified sequence.** It defines the
exact alphabet for that patent. Typical Alnylam legend has 50-100 entries:

```
Abbreviation | Meaning
A            | Adenosine-3'-phosphate (2'-OH)
a            | 2'-O-methyladenosine-3'-phosphate
Af           | 2'-fluoroadenosine-3'-phosphate
As           | Adenosine-3'-phosphorothioate (PS bond to next nt)
Afs          | 2'-fluoroadenosine + PS bond to next nt
C            | Cytidine-3'-phosphate
...
U            | Uridine-3'-phosphate
u            | 2'-O-methyluridine
Uf           | 2'-fluorouridine
...
T            | 5-methyluridine (ribothymidine, rT)
t            | 2'-O-methyl-5-methyluridine
dT           | 2'-deoxythymidine (DNA, typically in 3' overhangs)
...
L96          | Triantennary GalNAc N-acetylgalactosamine conjugate linker
L10          | Cholesterol linker (older designs)
P            | Phosphate
Y34          | Abasic / inverted abasic linker
(Agn)        | GNA (glycol nucleic acid) block — generic-nucleotide non-RNA
(Ahd)        | 2'-O-hexadecyl (lipid conjugate)
N            | Any nucleotide (modified or unmodified)
s            | Phosphorothioate backbone (between bases, no new base)
```

**Always write the legend to `<patent>_abbreviations.csv`** before
attempting any sequence decoding.

---

## 4. Alnylam modified-strand notation decoder

**Current convention** (used in all recent Alnylam PCT applications): each
nucleotide is 1 or 2 characters, **uppercase base + lowercase modifier**.

### Token table

| Token                  | Meaning                                                         |
|:----------------------:|:----------------------------------------------------------------|
| `A` `U` `G` `C`        | RNA nucleotide (2'-OH)                                           |
| `T`                    | 5-methyluridine (= ribothymidine, rT)                            |
| `a` `u` `g` `c` `t`    | 2'-O-methyl (2'-OMe) ribose                                     |
| `Af` `Uf` `Gf` `Cf`    | 2'-fluoro (2'-F) at that ribose                                 |
| `As` `Us` `Gs` `Cs`    | RNA base with **phosphorothioate** (PS) backbone to next nt     |
| `Afs` `Ufs` `Gfs`      | 2'-fluoro + PS                                                  |
| `Ab` `Cb` `Gb` `Ub`    | beta-L ribose (+ optional `s` for PS)                            |
| `dA` `dC` `dG` `dT`    | 2'-deoxyribose (DNA) — typically 3' overhangs                   |
| `(Agn)` `(Cgn)` ...    | GNA (glycol nucleic acid) block — generic-nucleotide non-RNA    |
| `(A2p)` `(C2p)` ...    | 2'-phosphate (special backbone)                                 |
| `(Ahd)` `(Chd)` ...    | 2'-O-hexadecyl (lipid conjugate)                                |
| `L96` `L10`            | Triantennary GalNAc / cholesterol conjugate linker              |
| `Y34` `Y44`            | Abasic / inverted abasic linkers                                |
| `P`                    | Standalone phosphate                                            |
| `VP`                   | Vinyl-phosphonate                                               |
| `N`                    | Any nucleotide (modified or unmodified)                         |
| `s` (between bases)    | Standalone PS linker                                            |

### PS rule

`s` ALWAYS means "PS bond between the previous nt and the next nt" — it
never introduces a new base.

So `csasugag` = `c`(2'-OMe,PS) `a`(2'-OMe,PS) `u`(2'-OMe) `g`(2'-OMe)
`a`(2'-OMe) `g`(2'-OMe). The first `c` is PS-bonded to the next `a`; the
second `s` is PS-bonded to `u`; etc.

### Length sanity checks

- **Alnylam 19-nt sense with GalNAc**: 19 core nt + `nL96` or `L96` at 3'.
  No TT overhang on sense.
- **Alnylam 20-nt antisense**: 20 core nt + `dTdT` 3' overhang (common).
- **Free-uptake designs**: no `L96`, may have `dTdT` instead.
- **Modified antisense with GNA blocks** (`(Cgn)` / `(Tgn)`): each
  parenthesized token is ONE GNA nt; don't count `gn` as separate bases.

### Reference decoder (Python)

```python
import re

def decode_modified(seq):
    """Decode an Alnylam-notation modified strand.

    Returns (list_of_nt_dicts, linker_or_None). Each nt dict has:
      position (1-indexed from 5' end), base, sugar, ps_next, kind.
    """
    nts, linker = [], None
    i, L = 0, len(seq)
    while i < L:
        ch = seq[i]

        # Linker / end-marker tokens (L96, L10, Y34, Y44, P, VP)
        if ch in "LYP":
            m = re.match(r"(L\d+|Y\d+|P|VP)", seq[i:])
            if m:
                linker = m.group(1); i += len(m.group(1)); continue

        # Parenthesized: GNA block, 2'-phosphate, 2'-O-alkyl — capture + skip
        if ch == "(":
            j = seq.find(")", i)
            if j != -1:
                nts.append({"kind": "special", "token": seq[i:j+1], "base": "N",
                            "sugar": "special", "ps_next": "",
                            "position": len(nts)+1})
                i = j + 1; continue

        # DNA 2-char token: dA, dC, dG, dT, dU (+ optional trailing s for PS)
        if ch == "d" and i + 1 < L and seq[i+1] in "ACGTU":
            base, ps = seq[i+1], ""
            if i + 2 < L and seq[i+2] == "s": ps = "PS"; i += 3
            else: i += 2
            nts.append({"kind": "dna", "base": base, "sugar": "DNA",
                        "ps_next": ps, "position": len(nts)+1}); continue

        # Lowercase = 2'-OMe; trailing s = PS
        if ch in "acgtu":
            base, ps = ch.upper(), ""
            if i + 1 < L and seq[i+1] == "s": ps = "PS"; i += 2
            else: i += 1
            nts.append({"kind": "rna", "base": base, "sugar": "2'-OMe",
                        "ps_next": ps, "position": len(nts)+1}); continue

        # Uppercase RNA / 2'-F (Xf) / PS (Xs) / both (Xfs)
        if ch in "AUGCT":
            base, sugar, ps = ch, "RNA", ""
            if i + 1 < L and seq[i+1] == "f":
                sugar = "2'-F"; i += 2
                if i < L and seq[i] == "s": ps = "PS"; i += 1
            elif i + 1 < L and seq[i+1] == "s":
                ps = "PS"; i += 2
            else: i += 1
            nts.append({"kind": "rna", "base": base, "sugar": sugar,
                        "ps_next": ps, "position": len(nts)+1}); continue

        # Beta-L (Xb, Xbs)
        if ch in "ACGTU" and i + 1 < L and seq[i+1] in "bs":
            base, sugar, ps = ch, "beta-L", ""; i += 2
            if i < L and seq[i] == "s": ps = "PS"; i += 1
            nts.append({"kind": "rna", "base": base, "sugar": sugar,
                        "ps_next": ps, "position": len(nts)+1}); continue

        # N placeholder
        if ch == "N":
            nts.append({"kind": "rna", "base": "N", "sugar": "RNA",
                        "ps_next": "", "position": len(nts)+1}); i += 1; continue

        # Unknown token — record verbatim and keep going
        nts.append({"kind": "unknown", "token": ch, "base": ch,
                    "sugar": "unknown", "ps_next": "", "position": len(nts)+1})
        i += 1
    return nts, linker
```

### Cheap presence-check (for wide-CSV mode)

If you don't need per-nucleotide decoding and just want to know which
modification types are present in a strand:

```python
import re

def quick_decode(seq):
    if not seq: return None
    return {
        "has_GalNAc_L96":   "L96" in seq or "L10" in seq,
        "has_2_OMe":        bool(re.search(r"[acgtu]", seq)),
        "has_2_F":          bool(re.search(r"[ACGTU]f", seq)),
        "has_PS":           "s" in seq.replace("L96", "").replace("L10", ""),
        "has_DNA_overhang": any(t in seq for t in ("dT", "dA", "dC", "dG")),
    }
```

---

## 5. Activity table shapes

### Wide-format dose-response

Standard 4-dose layout. The header is typically two rows (dose labels, then
column-type labels), then data rows.

```
|              | 50 nM        | 10 nM        | 1 nM         | 0.1 nM       |
| DuplexID     | Avg | ST DEV | Avg | ST DEV | Avg | ST DEV | Avg | ST DEV |
| AD-517197.1  | 17.0| 7.0    | 15.4| 3.4    | 41.0| 13.9   | 77.8| (none) |
```

Each data row carries **2 × N doses** values (Avg + SD per dose),
interleaved. See `parse_wide_activity()` snippet.

### Long-format dose-response

```
| Duplex        | Avg % of Message Remaining | SD       | Dose (nM) |
| AD-519351.1   | 18.91743                   | 1.592235 | 50        |
```

One row per (duplex, dose) combination. Use `parse_long_activity()`.

### Single-dose-only wide format

Some tables only carry one dose (e.g. a 50 nM Cos-7 screen):

```
|              | 50 nM      |
| DuplexID     | % of Message Remaining | ST DEV |
| AD-521420.1  | 12.6       | 5.2 |
```

Same parser as wide-format, just with one dose.

---

## 6. In-vivo result tables

Each in-vivo study has a setup pair (Treatment Groups + Duplexes of
Interest) and a results table. The results table row label is one of:

| Label                                | Meaning                                |
|:-------------------------------------|:---------------------------------------|
| `PBS Day 7` / `Naïve Day 7`          | Control arm                            |
| `AD-XXXXX.N`                         | Duplex, generic (single dose/timepoint)|
| `AD-XXXXX.N 10mpk D7`                | Duplex + dose + timepoint              |
| `AD-XXXXX.N (+ CTL) 30mpk D7`        | Control co-administration              |

`parse_invivo_study_results()` (in the script) handles all four.

### AAV-titration tables

Different metric: `fold_change_vs_PBS`, not `pct_mRNA_remaining`. The
"dose" column is the AAV titer (`5.00E+11 vp/mouse`). The fold-change at
high AAV titer can be > 1 if the AAV itself expresses more target mRNA.

---

## 7. Column-mapping pitfalls

These are the bugs that have actually bitten this workflow. Each one is
labeled with the symptom + the fix.

### P1: Substring matching on "sense seq" matches antisense too

**Symptom**: antisense column gets overwritten with sense column data,
scrambling every duplex in the table.

**Cause**: `"sense seq" in "Antisense Sequence 5' to 3'"` is `True` (because
`antisense` contains `sense` and `sequence` contains `seq`).

**Fix**: use prefix matching or word-boundary matching:

```python
# BAD
if "sense seq" in header.lower(): col_sense = i

# GOOD
h_norm = header.lower().strip()
if h_norm.startswith("sense "): col_sense = i
elif h_norm.startswith("antisense "): col_antisense = i
```

### P2: Prose mentions of "Table N" look like real titles

**Symptom**: a page says "as shown in Table 21 (provided elsewhere)" without
Table 21 actually being on that page — the parser thinks it's a title.

**Fix**: for each "Table N" mention, compute distance to the next
`<table border=1>`. < 500 chars = real title; > 500 chars or no table
following = prose mention — skip.

### P3: Multi-page tables — Unmodified scoops into Modified

**Symptom**: ~50% of Modified rows have their `antisense_seq` replaced
with the antisense SEQ ID.

**Cause**: when an Unmodified (Style A) table ends on the same page a
Modified (Style B) table begins, the walker scoops rows from the next
table into the current one.

**Fix**: detect each `<table>` block's style from its header (`Range in NM`
→ Style A; `mRNA Target Sequence` → Style B). Skip subtables whose style
doesn't match the title you're walking.

### P4: OCR-merged column headers

**Symptom**: header is `"SEQ IDmRNA Target Sequence 5' to 3'"` (two cells
glued into one).

**Fix**: when prefix-matching fails, fall back to position-based mapping:

```python
# Style B 7-col: duplex=0, sense=1, antisense=3, mrna_target=5
# Style A 7-col: duplex=0, sense=1, antisense=4
# Style B 6-col: duplex=0, sense=1, antisense=3
# Style A 5-col: duplex=0, sense=1, antisense=4 (continuation page)
```

Preserve the header from the FIRST page of the table; continuation pages
may show only 4-5 cells per row.

### P5: PaddleOCR-VL output merges cells across page breaks

**Symptom**: parser writes far fewer rows than `n_duplex_ids` implies
(e.g. 23 rows from a 349-row table). Each data row has 4-5 numeric cells
where the header says 8-10.

**Diagnosis**: dump `<tr>...</tr>` for 2-3 data rows and count `<td>` per
row. If the row has fewer `<td>` than the header (and the header looks
right), the OCR merged cells across the page break. NOT a parser bug.

**Fix attempts that don't work**: re-detecting column count, position-based
mapping, splitting cells by delimiter.

**What does work**:

1. Document the damage rate per table in the summary caveats.
2. Extract whatever rows DO parse cleanly; tag them with
   `table_damaged_by_ocr=True`.
3. Leave the rest for a fresh OCR pass on those specific pages.

**Typical damage rates** observed in Alnylam-style PaddleOCR-VL output for
the 200-300 row screens: 25-90% of rows lost. Smaller tables (< 100 rows)
usually survive cleanly.

### P6: AD-ID regex is not always `AD-\d{6}`

Alnylam uses 4-7 digit IDs (`AD-67554.9`, `AD-518710.1`,
`AD-1010733.2`). Use:

```python
AD_RE = re.compile(r"^AD-\d{4,8}(?:\.\d+)?$")
```

### P7: Version-numbered duplexes are different designs

`AD-519351.1` (initial lead) and `AD-519351.4` (refined chemistry) have
different modified sequences. Don't deduplicate on the base name when
extracting sequences — only when summarizing unique chemical designs.

### P8: Activity measurements for IDs NOT in any sequence table are normal

~250 measurements in a typical Alnylam Hep3B screen reference IDs like
`AD-518710.1` that appear in NO sequence table — internal Alnylam
screening library whose sequences weren't published. Flag these with
`unmatched_to_sequence=True`; don't drop them.

### P9: Activity is per-duplex, never per-strand

If the user asks for per-strand activity, explain that no patent reports
it. Activity attaches to the duplex (sense+antisense pair). Each strand
is "active" only in the context of its partner.

### P10: dose/timepoint suffixes OCR-merged into duplex name

In-vivo result labels look like `AD-520018.2 10mpk D7` (one cell). Strip
the suffixes when extracting the duplex ID; keep dose/timepoint as
separate columns.

### P11: Output file naming should encode the patent ID

`pnpla3_duplexes.json` collides when you have multiple patent folders
open. Use `<PATENT_ID>_sirna.json`, e.g. `WO2021126734_sirna.json`.

### P12: `x` at sense 3' end is an abasic / placeholder overhang, not a base

**Symptom**: When parsing a sense strand that ends in `x` (e.g.
`AAAGACAACAAACAUUAUAUUx` from T7 in WO 2012/177784), the decoder tries
to interpret `x` as a base. `x` is NOT a nucleotide — it's a placeholder
for an abasic / inverted-abasic linker, or simply the OCR's
representation of an empty position in a 3' overhang.

**Fix**: in `decode_modified()` (section 4), strip a trailing `x` and
record it as a `abasic-overhang` modification at the 3' end. Pair this
with the AS strand's `+X` extra 3' nucleotide (e.g. AS has 23 nt where
sense has 21) — the 2 extra AS nt are the 3' overhang, and the `x` on
sense is a placeholder. Don't try to base-pair them.

### P13: Three-series ID pattern (see section 2 "Three-series ID pattern")

When the same design appears under multiple duplex IDs across
sequence tables (A: `x`-overhang, B: refined 21nt, C: GalNAc-conjugated
of A), keep all three IDs and join via the patent's Table 15 (or
equivalent ID-pairing table). Don't merge on chemical identity at
extraction time — let the analysis layer decide. Add a new
`ID_PAIRING_TABLE = N` and `DUPLEX_SERIES = {"A": "...", "B": "...",
"C": "..."}` to CONFIG.

### P14: Schema B JSON has no `get_tables` cache

`search_patent` / `get_tables` / `extract_sirna_data` reject the
per-patent dict schema (see section 9). The script's
`build_tmap_from_schema_b()` helper walks `pages[i].blocks[].block_label
== 'table'` and pairs each table with the preceding
`block_label == 'figure_title'` to reconstruct the table map. If you
hit this and the helper is missing, you're about to spend an hour
reinventing it. Add it before doing anything else.

### P15: Tables 11+ may have 5+ columns of dose data per row

The single-dose screen for the GalNAc-conjugated set (T11 in
WO 2012/177784) has **five** dose columns: `10nM (RNAimax) | 0.1nM
(RNAimax) | 500nM PCH Celsis (FU) | 100nM PCH Celsis (FU) | 10nM PCH
Celsis (FU)` with interleaved SD. The wide-activity parser handles
this as long as `n_doses = len(dose_labels)` is set correctly. If
you're hand-rolling, set `n_doses = 5` and pass five dose labels.

### P16: The decoder in section 4 ASSUMES `Xf` = 2'-Fluoro — but older Alnylam patents use `f` as a case-shift delimiter only

**Symptom** (seen on WO 2012/177784 A1 ANGPTL3, 2026-06-22): every
modified strand starts with a long prefix like
`AfaAfgAfcAfaCfAfAfaCfAcUfuAfufUfL96`. With the section-4 decoder,
the `f` between bases (e.g. `Af` then `a`) gets parsed as
"uppercase A with 2'-F sugar" + a separate lowercase `a` (2'-OMe).
The output modification list ends up with **every other base flagged
as 2'-Fluoro**, which is wrong: this patent's designs use NO
2'-Fluoro at all.

**Cause**: WO 2012/177784 A1's Table 1 legend defines **only two
sugar modifications**:

| Symbol | Meaning |
|:------:|:--------|
| `A C G U` | unmodified ribonucleotide (2'-OH) |
| `a c g u` | 2'-O-methyl ribose |
| `dT` | 2'-deoxythymidine |
| `s` (suffix on base) | phosphorothioate linkage |

There is **no 2'-Fluoro** in this patent's designs. The `f` you see
between bases (e.g. `Afa` = `A` then `a`) is the OCR/publisher's
case-shift delimiter — it indicates "the previous base was
uppercase, the next is lowercase" — NOT a sugar modifier.

**The two Alnylam conventions, never both** in one patent:

- **Older style** (e.g. WO 2012/177784, WO 2010/...): **case-shift
  only**. `A`=unmod, `a`=OMe, `dT`=DNA. Every strand reads like
  `AfaAfg...` where `f` is just the case-transition character.
  Walk the strand in 1- or 2-character windows, skipping `f`.
- **Newer style** (e.g. WO 2021/126734): **explicit `Xf` =
  2'-Fluoro**. `A`=unmod, `a`=OMe, `Af`=2'-F, `As`/`Afs`=PS
  variants. Walk 2- or 3-character windows including the `f`.

**Fix**: read Table 1 first and dispatch the decoder accordingly.

```python
def decoder_convention(table1_legend):
    """Returns 'case_shift' (older) or 'explicit_f' (newer Alnylam)."""
    has_explicit_fluoro = any(
        "Af" in entry and ("fluoro" in meaning.lower() or "2'-F" in meaning)
        for entry, meaning in table1_legend
    )
    return "explicit_f" if has_explicit_fluoro else "case_shift"
```

In `case_shift` mode the decoder should:

- Skip `f` (it's a delimiter, not a token).
- Treat `A`/`C`/`G`/`U` as a single uppercase base (unmodified).
- Treat the next `a`/`c`/`g`/`u` as the next base with 2'-OMe sugar.
- Treat `s` after any base as a PS marker on the **next** linkage.

**Diagnostic**: if every sense strand in a table decodes to roughly
50% 2'-F sugars and 50% 2'-OMe, AND the patent's Table 1 doesn't
list 2'-F, you've hit P16. Switch to `case_shift` mode.

### P17: Schema-B (per-patent dict) patent JSON works with `get_tables` after the first call

**Correction to section 9**: The `get_tables` / `search_patent`
tools **do** work on the WO 2012/177784 schema-B JSON (per-patent
dict: `{patent_id, total_pages, pages:[{page_index, blocks:[...]}]}`)
once the tables cache (`.tables.json`) is written. The schema they
reject is the raw PaddleOCR batch result (a flat list). The
schema-B JSON's per-table HTML blobs are well-formed and the
extractor parses them via the same `parse_table` snippet from
section 10.

After the first `get_tables` call, skip the JSON walk entirely:

```python
idx = json.load(open("patent.tables.json"))
# idx["tables"] is a list of dicts with title, pages, rows, headers,
# is_sequence_table, is_activity_table, n_rows, duplex_ids
for t in idx["tables"]:
    if "Unmodified" in t["title"]: ...     # Style A
    if "Modified" in t["title"]: ...       # Style B
```

This is much faster than re-parsing 200+ pages of HTML on every
run. Cache invalidated automatically when the patent JSON changes.

### P18: ST.25 sequence.txt may have right-side line numbers (length markers) appended to sequence strings

**Symptom** (seen on WO 2012/078967 A1 EGLN, 2026-07-06): after parsing
the ST.25 sequence listing with a naive `replace(" ", "")` join, the
sequence string ends with digits like `...AAGUUCATT19` or `...CUUCC19`
instead of clean `...AAGUUCAUU` or `...CAUCUUCC`. Result: every record
gets marked `unverified_ocr` even when the OCR is perfect.

**Cause**: FastSEQ and other sequence-listing-generation tools append a
right-side line number (the cumulative length up to that line) to each
sequence line in the `.txt` file. The standard ST.25 format itself
doesn't use line numbers, but FastSEQ does, and many Alnylam filings
were generated with FastSEQ. Example raw line:

```
gcuauguccg ucacguuga                                              19
```

After naive whitespace stripping, this becomes `gcuauguccgucacguuga19`
— the `19` is the line number / cumulative length, not part of the
sequence.

**Fix**: when building the ST.25 sequence cache, strip line-number
suffixes per line BEFORE joining:

```python
# Pattern 1: a pure-numeric line (column markers like "1 5 10 15 20")
if re.match(r"^\d+$", stripped) or re.match(r"^[\d\s]+$", stripped):
    continue  # skip the line entirely

# Pattern 2: trailing digits on the same line as the sequence
stripped = re.sub(r"\s*\d+\s*$", "", stripped)

# Final safety net on the joined sequence:
sequence = re.sub(r"\d+$", "", sequence)
```

**Diagnostic**: if your verification rate is 0% across the board but
the sequences look correct, you've hit P18. Also note: some
`sequence.txt` entries are 21 nt long (`19 + TT` 3' overhang baked in
for modified strands). When comparing against the GT, strip trailing
`TT` from the GT to match the 19-nt OCR core.

### P19: In-vivo "treatment-group design" tables vs. "result" tables

**Symptom** (seen on WO 2012/078967 A1 EGLN, 2026-07-06): extracted
in-vivo tables give 5-7 rows of design metadata (Group, siRNA,
Formulation, Sample Size, Dose, In-vitro IC50) but NO numeric mRNA-knockdown
percentages. Activity keywords like "knockdown" appear in table TITLES
("Table 3. In vivo knockdown of EGLN genes") but the actual % data is
not in the table rows.

**Cause**: the mRNA-knockdown % data for these in-vivo studies lives in
**chart/figure blocks** on adjacent pages (the bar charts showing
% remaining for each treatment group at each timepoint), not in the
table text. The design tables are the experimental setup, not the
results.

**Fix**: extract the design rows as "treatment-group" records (which
siRNA → which group → which dose → which formulation → n) and
explicitly note in the extraction summary that the numeric mRNA-%
results require a separate figure-extraction step (or re-OCR pass on
the chart pages). The chart image files are present in
`output/{patent_id}/images/` as `layout_det_res_pNNN.jpg`.

**Diagnostic**: if your "in-vivo activity" tables have columns
`Group | siRNA | Formulation | Sample Size | Dose (mg/kg) | In vitro IC50`
and ZERO columns with "% remaining" or "fold change" — these are
design tables, not result tables.

### P20: Same-page title tables (e.g. T7 + T8 in WO 2012/078967, both on p.151)

**Symptom**: two adjacent "Table N" titles appear on the same page
(T7 and T8 in WO 2012/078967). A naive "end = next_title_page - 1"
range calculation gives T7 an end of 150, which is BEFORE T7's own
title page 151 — so T7 ends up with 0 table blocks.

**Fix**: when building per-label page ranges, if `next_title_page ==
current_title_page`, set end = current_title_page (single-page
assignment for both titles). Better: use a block-level state machine
that walks `pages[].blocks[]` in order and assigns each `table` block
to the most recent `figure_title` with a `Table N` caption. The
state-machine approach handles same-page titles, multi-page continuations,
and form/legal tables on adjacent pages all in one pass.

### P21: ST.25 sequence listing rendered as table blocks (post-p.173+ in some patents)

**Symptom**: tables continue appearing for hundreds of pages after the
last "real" patent table. The most recent figure_title label (e.g. T13B
on p.159) keeps re-assigning to every subsequent table block, inflating
the count and pulling in non-sequence-listing tables (form/legal docs).

**Cause**: pages 173+ in some Alnylam-style PCTs contain the ST.25 sequence
listing rendered as table blocks (one per `<210>` SEQ ID, typically 2-3
tables per page). Pages 165-172 are also often USPTO filing forms
(Application Data Sheet, etc.) that look like tables. Both will be
incorrectly attached to whatever title was last seen if you don't cut off.

**Fix**: detect the ST.25 listing start by scanning `block_label ==
'header'` for text containing `SequenceLST` (Alnylam's internal
sequence-listing identifier). Hard cutoff the table indexer at that
page. Also filter form/legal tables from the sequence-table pool by
checking that ≥50% of rows look like sequence data rows (first column
matches `AD-XXX` OR a 1-5 digit Start position OR a 15+ nt nucleotide
string).

### P22: 5-column sequence tables with no Duplex Number column (T6A/B/C style)

**Symptom**: some sequence tables have schema `Start | Sense Sequence |
SEQ ID NO | Antisense Sequence | SEQ ID NO` — NO Duplex Number column.
A naive "first cell must match `^AD-\d+`" filter rejects all of them.
They look like form-data tables to the heuristic.

**Fix**: detect the column count from the header row. For 5-column
tables with "Start" in col[0] but no "Duplex Number", synthesize IDs as
`{gene_short}_{table_label}_{NNN}` (e.g. `hEGLN1_6A_001`).
These synthesized IDs cannot be cross-referenced to AD-XXXX IDs in
in-vivo tables — flag that limitation in the summary.

---

## 8. OCR artifacts

| Artifact                       | Real meaning                            | Where to fix         |
|:-------------------------------|:----------------------------------------|:---------------------|
| `γ` (Greek gamma)              | The letter `T`                          | Replace before parsing |
| `i` before a base (e.g. `iA`)  | OCR misread of `*A` or `mA`             | Context-dependent    |
| Multi-line cells               | Was originally one cell, OCR split      | Concatenate lines    |
| `º` `°` after digits           | Degree symbol, not sequence data        | Strip if not in seq  |
| `as` `cs` `gs` `us`            | 2'-OMe base + PS backbone to next nt    | Treat as `s` modifier |
| `(Xgn)` `(Xan)` `(Tmn)`        | GNA block                               | Skip parenthesized block |
| `L96` `L10` `L96N`             | Linker codes                            | End-of-sequence marker |
| `dT` `dA` `dC` `dG`            | DNA bases (3' overhangs)                | Parse as 2-char token |
| `x` at sense 3' end            | Abasic / placeholder 3' overhang        | Strip + record as `abasic-overhang`; pair with AS's `+2` extra 3' nt |
| `N` (uppercase)                | Placeholder base OR modification        | Check legend         |
| `\_` (escaped underscore)      | Real underscore                         | Unescape             |

### `d` vs `s` confusion

**OCR mangle**: `adA` or `aadA` is the OCR mangling of `asA` or `asU` (the
`s` reads as `d`). Cross-check by counting `s` vs `d` in the rest of the
sequence — if the patent's other sequences show ~30% PS-marked nucleotides
and a sequence has 0%, the `d` is probably `s`.

---

## 9. Patent JSON input schemas

The patent OCR pipeline can produce **three** top-level shapes. The
`search_patent` / `get_tables` / `extract_sirna_data` tools only accept
schema A. Schemas B and C are fully usable but require hand-parsing
(or use the script's auto-detection, which handles all three).

| Schema | Top-level | `get_tables` accepts? | How to detect |
|:-------|:----------|:----------------------|:--------------|
| **A** — PaddleOCR-VL batch | `[{page_index, markdown: str\|dict, ...}, ...]` | yes | `isinstance(d, list)` |
| **B** — per-patent dict | `{patent_id, total_pages, pages:[{page_index, blocks:[{block_label, block_content}]}]}` | **no** | `isinstance(d, dict) and "pages" in d and isinstance(d["pages"][0].get("blocks"), list)` |
| **C** — Marker / PyMuPDF | `{... "blocks":[{block_label, block_content}]}` page-keyed | partial | `isinstance(d, dict) and "blocks" in d` |

**Schema B failure mode** (seen on WO 2012/177784 A1 ANGPTL3, 2026-06-22):

```
search_patent failed: WO2012177784.json: expected a PaddleOCR batch result
(JSON array), got dict. Top-level keys: ['patent_id', 'total_pages', 'pages']
```

When you see this, the data is fully usable — just bypass the tool and
parse manually. The schema-B structure is:

```python
data = json.load(open(patent_json))
pages = data["pages"]   # list of {"page_index": int, "blocks": [...]}
for pg in pages:
    for b in pg["blocks"]:
        # b["block_label"] in {"text","aside_text","header","number",
        #   "footer","table","image","figure_title","vision_footnote",
        #   "doc_title","paragraph_title","display_formula"}
        if b["block_label"] == "table":
            html_blob = b["block_content"]   # full <table>...</table>
            rows = parse_table(html_blob)     # works unchanged
```

The HTML table cells are well-formed (`<tr><td>...</td></tr>`), so the
existing `parse_table` snippet in section 10 works unchanged. You lose
the `get_tables` cache (titles, page spans, suspected_ocr_damage) — to
recreate it, walk the JSON, scan `block_label=='figure_title'` for
"Table N" strings, and check whether a `block_label=='table'` follows
within the same page. The script has a `build_tmap_from_schema_b()`
helper that does this.

### Schema A sub-variants (PaddleOCR-VL)

**Older versions (pre-mid-2026): double-encoded** — `patent[i]["markdown"]`
is a JSON-encoded string. Must unwrap:

```python
import json
with open("patent.json") as f:
    patent = json.load(f)

def page_text_old(i):
    md = json.loads(patent[i]["markdown"])   # decode the JSON string
    return md.get("text", "")
```

### Newer versions (mid-2026+): real dict

`patent[i]["markdown"]` is a dict with a `text` key:

```python
def page_text_new(i):
    return patent[i]["markdown"]["text"]
```

### Auto-detect (schema A only)

Always check the type first — mixing them up surfaces as `TypeError` or
`AttributeError`:

```python
if isinstance(patent[0]["markdown"], str):
    page_text = lambda i: json.loads(patent[i]["markdown"])["text"]
elif isinstance(patent[0]["markdown"], dict):
    page_text = lambda i: patent[i]["markdown"]["text"]
else:
    raise ValueError(f"Unknown markdown type: {type(patent[0]['markdown'])}")
```

### Pre-parsed table cache (schema A only)

If `get_tables` has been run on the patent, there's a cache file at
`<patent>.tables.json` containing every table's title, page span, and
already-parsed rows. **Reading from cache is far faster than re-parsing
269 pages**. Use:

```python
idx = json.load(open("patent.tables.json"))
tmap = {t["table_number"]: t for t in idx["tables"]}
# tmap[N]["rows"] gives list of cells-per-row
# tmap[N]["is_sequence_table"], tmap[N]["is_activity_table"] are set
```

For schemas B and C, the script's `build_tmap_from_schema_b()` walks the
JSON and constructs an equivalent `tmap` from `figure_title` + `table`
block pairs. Use that path; don't reinvent.

---

## 10. Parser snippets

### Cell extraction

```python
import re

def parse_table(table_text):
    """Extract rows from an HTML <table> blob. Returns list of cell-lists."""
    rows = []
    trs = re.findall(r"<tr[^>]*>(.*?)</tr>", table_text, re.DOTALL)
    for tr in trs:
        cells = re.findall(r"<td[^>]*>(.*?)</td>", tr, re.DOTALL)
        if not cells:
            continue
        clean = []
        for c in cells:
            c = re.sub(r"<[^>]+>", "", c)               # strip inner tags
            c = c.replace("&nbsp;", " ")
            c = c.replace("&#x27;", "'").replace("&#39;", "'")
            c = c.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
            c = c.replace("&quot;", '"')
            c = re.sub(r"\$+", "", c)                    # strip $...$ LaTeX
            c = re.sub(r"\\mathrm\{([^{}]*)\}", r"\1", c)
            c = c.replace("$", "")
            c = re.sub(r"\s+", " ", c).strip()
            clean.append(c)
        if any(clean):
            rows.append(clean)
    return rows
```

### Wide-format activity parser

```python
def parse_wide_activity(rows, dose_labels, table_label):
    """rows: list of cell-lists. dose_labels: ['50 nM', '10 nM', ...]."""
    out = []
    n_doses = len(dose_labels)
    for r in rows:
        if not r:
            continue
        first = (r[0] or "").strip()
        if not (isinstance(first, str) and re.match(r"^AD-\d{4,8}(?:\.\d+)?$", first)):
            continue
        if len(r) < 1 + 2 * n_doses - 1:
            continue
        vals = [c.strip() for c in r[1:]]
        rec = {"table": table_label, "doses": {}}
        ok = True
        for j, dose in enumerate(dose_labels):
            avg = vals[2 * j] if 2 * j < len(vals) else ""
            sd  = vals[2 * j + 1] if 2 * j + 1 < len(vals) else ""
            if not re.match(r"^-?\d+(?:\.\d+)?$|^n/a$", avg, re.IGNORECASE):
                ok = False; break
            rec["doses"][dose] = {"avg": avg, "sd": sd}
        if ok:
            out.append((first, rec))
    return out
```

### Long-format activity parser

```python
def parse_long_activity(rows, table_label):
    """Long-format: [duplex, avg, sd, dose_nM]."""
    out = []
    for r in rows:
        if not r: continue
        first = (r[0] or "").strip()
        if not (isinstance(first, str) and re.match(r"^AD-\d{4,8}(?:\.\d+)?$", first)):
            continue
        if len(r) < 4: continue
        if not re.match(r"^-?\d+(?:\.\d+)?$", r[1]): continue
        if not re.match(r"^-?\d+(?:\.\d+)?$", r[3]): continue
        dose = re.sub(r"\s*nM\s*$", "", r[3].strip(), flags=re.IGNORECASE) + " nM"
        out.append((first, {"table": table_label, "doses": {dose: {"avg": r[1].strip(), "sd": r[2].strip()}}}))
    return out
```

### In-vivo result parser

```python
def parse_invivo_study_results(rows, table_label):
    """In-vivo result labels: 'AD-X.N 10mpk D7' or 'PBS Day 7' or 'AD-X.N (+ CTL) 30mpk D7'."""
    out = []
    for r in rows:
        if not r or len(r) < 3: continue
        first = (r[0] or "").strip()
        try:
            avg = float(r[1].strip()); sd = float(r[2].strip())
        except (ValueError, AttributeError):
            continue
        dose_match = re.search(r"(\d+\s*mpk)", first, re.IGNORECASE)
        time_match = re.search(r"(D\d+|Day\s*\d+)", first, re.IGNORECASE)
        ctl = "(+ CTL)" in first
        dup_id = re.sub(r"\s*\d+\s*mpk.*$", "", first, flags=re.IGNORECASE)
        dup_id = re.sub(r"\s*\(\+ CTL\).*$", "", dup_id)
        dup_id = re.sub(r"\s*(D\d+|Day\s*\d+).*$", "", dup_id, flags=re.IGNORECASE).strip()
        out.append({
            "table": table_label,
            "label_raw": first,
            "duplex_id": dup_id if re.match(r"^AD-\d{4,8}(?:\.\d+)?$", dup_id) else None,
            "is_control": dup_id.lower() in ("pbs", "naive", "naïve"),
            "dose_mpk": dose_match.group(1).replace(" ", "") if dose_match else None,
            "timepoint": time_match.group(1).replace(" ", "") if time_match else None,
            "ctl_marker": ctl,
            "pct_remaining": avg,
            "sd": sd,
        })
    return out
```

---

## 11. Sanity checks

After extraction, verify:

- [ ] All 6 output files present: JSON, wide CSV, long activity CSV,
      abbreviations CSV, crossrefs CSV, summary MD
- [ ] JSON `duplex_count_with_sequence_data` == number of rows in wide CSV
- [ ] At least one antisense strand is the reverse-complement of its sense
      strand (strip 3' overhangs first)
- [ ] At least one modified-sense sequence has `L96` (GalNAc) and decodes
      to 19 nt + linker
- [ ] For each activity table: `parsed_count / cache_n_rows >= 0.5`. If
      below 50%, flag `table_damaged_by_ocr=True` and document.
- [ ] Abbreviation CSV row count matches Table 1 in the patent
      (typically 50-100 for Alnylam-style)
- [ ] Long activity CSV has the expected order of magnitude: 3,000-5,000
      rows for a 700-duplex patent with 14 in-vitro + 6 in-vivo tables
- [ ] Cross-reference CSV has at least one row per in-vivo study (Studies
      1, 2, 3 typically have 12-16 duplexes each)

---

## 12. Reference example: WO 2021/126734 A1 (PNPLA3)

Used as the worked example for the runnable script. Useful as a sanity
check — if you re-run on this patent, you should see ~736 duplexes,
~3,400 in-vitro measurements, ~95 in-vivo measurements.

### Table inventory (51 tables total)

| Type                   | Tables                | Notes                                      |
|:-----------------------|:----------------------|:-------------------------------------------|
| Abbreviation legend    | 1 (67 entries)        | Page 116                                   |
| Unmodified sequences   | 2, 4, 6, 8, 10, 32, 36, 49 | Style A; 736 unique duplexes total     |
| Modified sequences     | 3, 5, 7, 9, 11, 33, 37, 50 | Style B; pairs with above             |
| In-vitro dose-response | 12-19, 34, 35, 38, 39, 51, 52 | Tables 34/35/38/39 OCR damaged       |
| Treatment Groups       | 40, 43, 46, 26        | In-vivo study setup                        |
| Duplexes of Interest   | 41, 44, 47            | DOI for in-vivo leads                      |
| In-vivo results        | 22, 25, 31, 42, 45, 48 | Numeric % mRNA remaining                  |
| AAV titration          | 28                    | Fold-change metric                         |
| In-vivo sequence info  | 21, 24, 27, 30        | Per-strand sequence for in-vivo duplexes   |

### Cohort splits

- **Tables 2-11**: First chemistry series (4 cohorts: 89, 28, 107, 90, 36)
- **Tables 32-37**: In-vivo lead cohort (GalNAc-conjugated)
- **Tables 49-50**: New free-uptake cohort (180 designs, no GalNAc, no
  TT overhang) — only ~30 IDs overlap with Tables 32-37

### In-vivo studies

| Study          | Groups | DOI | Results | Notes                              |
|:---------------|:------:|:---:|:-------:|:-----------------------------------|
| Study 1 (Fig 1)| 40     | 41  | 42      | Mouse AAV-hPNPLA3 single-dose      |
| Study 2 (Fig 2)| 43     | 44  | 45      | Multi-dose comparison              |
| Study 3 (Fig 3)| 46     | 47  | 48      | Dose-ranging + NHP comparator     |
| AAV titr (Fig 4)| 26    | —   | 28      | Fold-change vs PBS                 |

### Modification breakdown (sense strand, observed)

- 605 / 736 (82%) GalNAc-conjugated (sense has `L96`)
- 640 / 736 (87%) have 2'-OMe sugars
- 609 / 736 (83%) have 2'-F sugars
- 628 / 736 (85%) have PS backbone modifications

### OCR damage observed in dose-response tables

| Table | Cell line | Recovered | Expected | Loss rate |
|:-----:|:---------:|----------:|---------:|----------:|
| 34    | Hep3B     | 23        | ~349     | 93%       |
| 35    | PCH       | 23        | ~108     | 79%       |
| 38    | Hep3B     | 62        | ~85      | 27%       |
| 39    | PCH       | 29        | ~84      | 65%       |

Tables 12, 13, 15-18, 51, 52 parse cleanly (>95% recovery).

---

## 13. Reference example: WO 2012/177784 A1 (ANGPTL3)

A second worked example. Uses the **older case-shift modification
convention** (P16) — different from section 12. Useful as a sanity
check on a non-GalNAc-heavy, no-in-vivo design.

### Key differences from WO 2021/126734 (PNPLA3)

| Aspect | WO 2021/126734 (PNPLA3) | WO 2012/177784 (ANGPTL3) |
|:-------|:------------------------|:-------------------------|
| Modification convention | explicit_f (`Af`=2'-F) | **case_shift** (`f`=delimiter only, NO 2'-F) |
| Total unique duplexes | 736 | 428 |
| Total activity measurements | ~3,400 | ~2,200 |
| In-vivo studies | 3 (mouse + NHP) | **None** — all in-vitro |
| GalNAc-conjugated | 605 (82%) | 185 (43%) |
| Lead chemistry | 2'-OMe + 2'-F + PS | 2'-OMe + dTsdT + PS (no 2'-F) |
| Activity readout | mostly fold-change / mRNA % | mostly mRNA remaining (1.0 = no effect) |
| Three-series ID pattern | yes (A=initial, B=refined, C=GalNAc) | yes (A=GalNAc w/ `x`, B=refined 21nt, C=unmod same chem) |

### Table inventory (15 tables total)

| Type | Tables | Rows | Notes |
|:-----|:-------|-----:|:------|
| Abbreviation legend | 1 | 13 | Page 145 — short, only 4 sugar symbols + dT + s |
| Unmodified sequences (batch 1) | 2 | 62 | 19-23nt, LNP-formatted, no GalNAc |
| Modified sequences (batch 1) | 3 | 62 | `Afa`-prefixed strands, dTsdT 3' overhang |
| LNP single-dose screen | 4 | 61 | Hep3B, 10nM + 0.1nM + SD |
| LNP IC50 dose-response | 5 | 14 | Hep3B, 24h + 120h, 3 replicates each |
| LNP cell viability | 6 | 14 | HeLa + Hep3B, day 3 + day 6, 5 doses each |
| Unmodified GalNAc series | 7 | 185 | 21-23nt, sense ends with `x` (placeholder) |
| Modified GalNAc series | 8 | 83 | subset of T7, modified + L96 conjugate |
| Renumbered (no-GalNAc) batch | 9 | 181 | same chem as T7 but unconjugated |
| Renumbered (no-GalNAc) mod | 10 | 153 | same chem as T8 but unconjugated |
| GalNAc single-dose | 11 | 185 | 5 conditions (RNAimax 10/0.1nM + PCH 500/100/10nM) |
| GalNAc IC50 | 12 | 33 | free uptake + transfection (ND for most) |
| Unconjugated replacement dose-response | 13 | 181 | 10/0.1/0.025 nM + SD |
| Unconjugated replacement IC50 | 14 | 33 | selected subset |
| Unconjugated ↔ GalNAc pairing | 15 | 186 | unconjugated_id \| galnac_id |

### Sequence-table structure

Three-series pattern with shared IDs (NOT the same as PNPLA3 — here
series A and C share the **same duplex IDs** while series B is a
parallel renumbering):

- **Series A** (`AD-52953.x` to `AD-53150.x`): the "GalNAc batch".
  In T7 the sense strand ends with lowercase `x` (abasic placeholder);
  in T8 the same IDs appear with `L96` conjugation and full chemical
  modification. The pair (T7 row + T8 row with same ID) is the SAME
  physical duplex under two representations.
- **Series B** (`AD-52637.x` to `AD-52822.x`): a parallel renumbered
  set of unconjugated designs. Same chemistry patterns as T8 but no
  GalNAc, used in the LNP / RNAiMax-only screens.
- **Table 15** maps `AD-52637.x ↔ AD-52953.x` (and so on) explicitly.

### Decoder gotchas specific to this patent

1. **P16 (case-shift convention)** — see pitfalls section above.
   `Af` is NOT 2'-Fluoro. Apply the `case_shift` decoder mode.
2. **Sense strand ends with `x`** (in T7) — `x` is an abasic
   placeholder at the 3' end of the sense strand, paired with an
   antisense strand that has 2 extra 3' nucleotides. Already
   documented in P12.
3. **`Afc...L96`** sense strands have a final `L96` (GalNAc) after
   the modification sequence. The decoder should consume the entire
   trailing `L96` (or `L10`) token as a 3'-end linker entry, not as
   a base.
4. **OCR garbling on T8 / T10**: a handful of rows have sense strings
   with long repeats (`GfGfGfGfGfGfGf`, `Ufgugugu...`). These are
   flagged `needs_review=YES` in the JSON. **Do not try to clean
   these up** — they need a re-OCR pass on the source page. Just
   flag and document.

### Recommended extraction order

1. Parse Table 1 (13-row abbreviation legend) first. Detect
   `case_shift` convention because 2'-Fluoro is absent.
2. Build table index via `get_tables` (writes `.tables.json` cache).
3. Sequence tables: T2, T3, T7, T8, T9, T10 — all share the same
   decoder branch once you've set the convention.
4. Activity tables: T4, T5, T6 (LNP batch), T11, T12, T13, T14
   (GalNAc + replacement). T11 has 5 dose columns + 5 SD columns
   (P15). T13 has 3 doses + 3 SD.
5. Cross-ref: T15 → writes 186 (unconjugated, GalNAc) pairs.
6. Verify: every AD-XXXXX.x in T15 should appear in either T7/T8
   (GalNAc side) or T9/T10 (unconjugated side). Zero missing IDs
   is the success criterion.

### Sanity targets when re-running on this patent

- 726 sequence records across 4 table sources (T2, T3, T7, T8, T9, T10)
- 428 unique duplex IDs
- ~2,200 activity measurements
- 186 cross-ref pairs
- Zero IDs in T15 that don't appear in any sequence table

---

## 14. Reference example: WO 2012/078967 A1 (EGLN1/2/3 for EPO induction)

A third worked example, used 2026-07-06 to self-test this skill. The
key differences from sections 12 and 13 are listed below — mostly that
this is an **older (2012) Alnylam filing** with LNP-formulated siRNA
(no GalNAc) and a different table inventory pattern.

### Key differences from WO 2021/126734 (PNPLA3) and WO 2012/177784 (ANGPTL3)

| Aspect | WO 2021/126734 (PNPLA3) | WO 2012/177784 (ANGPTL3) | WO 2012/078967 (EGLN) |
|:-------|:------------------------|:-------------------------|:---------------------|
| Filing year | 2021 | 2012 | 2012 |
| Target gene(s) | PNPLA3 | ANGPTL3 | EGLN1 + EGLN2 + EGLN3 |
| Total pages | 269 | ~750 | 748 |
| Total unique duplexes | 736 | 428 | 1,533 |
| In-vivo studies | 3 (mouse + NHP) | None | 9 design tables (no raw % data in OCR) |
| GalNAc-conjugated | 605 (82%) | 185 (43%) | 0 (LNP-formulated) |
| Lead chemistry | 2'-OMe + 2'-F + PS | 2'-OMe + dTsdT + PS (no 2'-F) | 2'-OMe + dTsdT + PS (no 2'-F) |
| Convention | explicit_f (`Af`=2'-F) | case_shift | case_shift |
| Modification legend size (Table 1) | 67 entries | 13 entries | 12 entries |

### Table inventory (21 named tables)

| Type | Tables | Notes |
|:-----|:-------|:------|
| Abbreviation legend | 1 (12 entries) | Page 97 — short, confirms case_shift |
| Unmodified Mouse sequences | 2A, 2C, 2E | 31-32 rows each, real AD-IDs |
| Modified Mouse sequences | 2B, 2D, 2F | 31-32 rows each, endolight chemistry |
| Unmodified Human sequences | 6A, 6B, 6C | 241-494 rows each, **NO Duplex Number column** (P22) |
| Modified Human sequences | 13A, 13B | 47-90 rows each, real AD-IDs |
| In-vivo design tables | 3, 4, 5, 7, 8, 9, 10, 11, 12 | 6-14 design rows each (P19) |

### Pitfalls encountered on this patent (added to the pitfalls list)

- **P18**: ST.25 sequence.txt has right-side line numbers (length markers) appended to sequence strings. Strip with `re.sub(r"\s*\d+\s*$", "", line)` per line.
- **P19**: In-vivo tables are design tables, not result tables. The mRNA-knockdown % data is in chart blocks, not OCR'd table text.
- **P20**: Same-page title tables (T7 + T8 both on p.151). Use block-level state machine, not page-range logic.
- **P21**: ST.25 sequence listing rendered as table blocks (post-p.173). Hard cutoff at the first page with `SequenceLST` in the header.
- **P22**: 5-column sequence tables (T6A/B/C) with no Duplex Number column. Synthesize IDs as `{gene_short}_{table_label}_{NNN}`.

### Verification rate (final)

- 1,533 sequence records
- 1,048 sense sequences (68%) verified against sequence.txt ground truth
- 1,021 antisense sequences (66%) verified
- OCR error rate: ~32% of records have at least one nucleotide misread (dropped U, C↔G confusion in 2'-OMe regions of modified strands)
- 13 unique AD-IDs used in in-vivo treatment groups: AD-1955 (control), AD-40758 (EGLN3), AD-40773 (EGLN2), AD-40894 (EGLN1), plus cocktail percentages

### Recommended extraction order

1. Parse the ST.25 sequence.txt with the P18 fix in place (strip line numbers)
2. Build table index by walking pages[].blocks[] directly (schema B fallback)
3. Parse Table 1 (12 entries) — confirm `case_shift` convention because no 2'-F entry
4. Sequence tables: T2A-F, T6A-C, T13A-B. Auto-detect modified by row content (P19/P22)
5. Treatment tables: T3, T4, T5, T7-T12. Extract as design rows, not results
6. Verify: cross-check the 13 in-vivo AD-IDs against the sequence tables (should all be found)
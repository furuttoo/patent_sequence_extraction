# Patent layout judgement — what the code cannot see

The scripts in `scripts/` parse cells, decode notation, strip ST.25 line numbers
and validate every record. This file is the other half: the judgements that need
a reader who can see the page. Nothing here restates what the code already
enforces; each item names the evidence to look for and what to do when it is
missing.

Coverage markers (`[P2]`, `[TABLE-STYLE]`, ...) index the items below; a gap is
listed here only when no script can decide it.

---

## 1. Which table is this?  [TABLE-STYLE] [CONTINUATION]

The cache hands over rows and headers, not intent. Decide from the page:

- **Style A (unmodified, 7 columns)**: duplex id, sense seq, sense length,
  antisense seq, antisense length, target/gene, notes.
- **Style B (modified, 7 columns)**: duplex id, sense seq, antisense seq,
  target/gene. Lengths are absent because the notation carries them.
- 6- and 5-column variants exist (`[P22]`): the duplex column may be missing
  entirely and the id has to be synthesised from table label + row index.

**Continuation vs new table** `[CONTINUATION]`: a table running across a page
break usually repeats the header but not the title, and its row numbering
continues. A new table restarts numbering and carries its own title. Getting
this wrong duplicates rows; the code cannot tell them apart because both look
like "a few rows with a header".

## 2. Identifiers  [P6] [P7] [P10] [P13]

- `[P6]` `AD-\d{6}` is the common shape, not the only one. `AD-12115.1`,
  `AD-1955`, `AL-DP-4440`, `AV00999.um` all appear. Never normalise or truncate
  an identifier to fit a regex; keep what is printed.
- `[P7]` Version suffixes (`.1`, `.2`) are different designs, not typos of one
  duplex. Do not merge them.
- `[P10]` Dose or timepoint text is sometimes OCR-merged into the identifier
  (`AD-40758 30mpk 4wk`). Split it: the id is the compound, the rest is the
  measurement condition.
- `[P13]` Newer filings carry three parallel series (SEQ ID, duplex, compound).
  They are not interchangeable; record which one the table actually prints and
  keep the other two in `source_quote`.

## 3. Activity table shapes  [ACT-WIDE] [ACT-LONG] [DESIGN-VS-RESULT]

- `[ACT-WIDE]` One row per duplex, a dose per column, a cell line per column
  group. The column header carries the condition, not the row.
- `[ACT-LONG]` One row per measurement: duplex, cell line, dose, metric, value.
- `[DESIGN-VS-RESULT]` Treatment-group tables (which animals got which compound)
  are **design**, not measured activity. Only result tables produce records.
  The two can sit on the same page and share a title prefix.

## 4. Whose activity is it?  [ACT-PER-DUPLEX] [ORPHAN-ACTIVITY] [DUP-MEASUREMENTS]

- `[ACT-PER-DUPLEX]` Activity belongs to the duplex, never to one strand. Both
  strand records of that duplex may carry it, but the value must not be counted
  twice when collapsing.
- `[ORPHAN-ACTIVITY]` An activity table routinely lists compounds that appear in
  no sequence table (controls, prior art, screening hits). That is normal: write
  the measurement only if the sequence is also available somewhere in the same
  patent, otherwise leave the row out of the duplex export.
- `[DUP-MEASUREMENTS]` The same identifier can appear several times: different
  concentrations, timepoints, cell lines, or genuine repeats. Only identical
  conditions are duplicates; the rest are distinct measurements.

## 5. OCR damage: recognise it, do not repair it  [P4] [P5] [P12] [OCR-SHAPE]

- `[P4]` Merged column headers: two header cells collapse into one
  (`Sense sequence / Antisense sequence`), shifting every column by one.
- `[P5]` PaddleOCR-VL sometimes merges cells across a page break, so a row ends
  in the middle of a sequence and continues on the next page.
- `[P12]` A terminal `x` is an abasic or placeholder overhang, not a base. It
  must not enter `sequence`.
- `[OCR-SHAPE]` Known damage shapes: dropped `U`, `C`↔`G` in 2'-OMe runs,
  `d`↔`s` confusion. When a sequence looks broken, keep it, set
  `needs_review`, and quote the printed row; never correct a base from context.

## 6. Boundaries between tables and prose  [P2] [P3] [P20] [P21]

- `[P2]` Prose that mentions "Table 4" is not a table title. Titles are bold or
  centred and sit adjacent to the header row.
- `[P3]` A multi-page unmodified table can visually absorb the beginning of a
  modified table on its last page. Check whether the column count changes
  mid-page.
- `[P20]` Two tables can share one page (e.g. T7 and T8 on p.151). Track the
  title block, not the page range.
- `[P21]` A rendered ST.25 listing looks like a table block but is a sequence
  listing: bases only, no duplex ids, no activity. Prefer the listing parser and
  stop treating those pages as sequence tables at the first `SequenceLST`
  heading.

## 7. Precedence and hand-off  [PRECEDENCE] [NO-GUESS]

- `[PRECEDENCE]` The patent's own legend (`abbrev_map`, `step6_legends.json`)
  wins over every default. Defaults apply only to symbols the patent does not
  define.
- `[NO-GUESS]` Two things are never inferred: the meaning of letter case
  (`detect_case_convention` returns `unknown` rather than assuming 2'-O-Methyl)
  and the bases behind a damaged cell. Both end in `flagged` / `needs_review`
  with the printed evidence in `source_quote`.

---

## What is deliberately absent

Everything the code already owns: notation token tables, position arithmetic,
ST.25 line-number stripping, sequence alphabet, field validation, CSV layout.
Look those up in `modification_extractor.py`, `parse_sequence_listing_st25.py`
and `schemas_pydantic.py` instead of duplicating them here.

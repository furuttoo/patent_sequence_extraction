# Modification Decoding Reference

Agent reference for interpreting siRNA chemical modification shorthand found
in OCR'd patent sequence tables and text.

This reference is used **after** the abbreviation map has been built from the
patent's own legend pages (see `abbrev_map` in the page cache). When a symbol
appears in `abbrev_map`, that mapping takes precedence over the defaults below.
Use these defaults only for symbols the patent does not define locally.

---

## Decoding priority

1. **`abbrev_map` in the page cache** — built by reading the patent's own
   abbreviation legend table (pages flagged `ABBREV`). Always consult this
   first. The patent is the ground truth for its own notation.
2. **Style A defaults below** — apply when a modifier suffix (`m`, `f`, `n`,
   `s`) appears in a sequence but the patent has no legend, or the legend does
   not cover that symbol.
3. **Style B convention** — lowercase letter = ribonucleotide modified per the
   case convention stated in the patent's footnote. Read the footnote on the
   sequence page to confirm before assuming.

---

## Style A — inline suffix codes

Modifier characters appear directly appended to the nucleotide letter within
the sequence string.

### Sugar (2') modifications

| Symbol pattern | Decoded type | Notes |
|---|---|---|
| `Am`, `Cm`, `Gm`, `Um`, `Tm` | 2'-OMe | Most common Alnylam pattern; lowercase variants `am`, `cm` etc. identical |
| `mA`, `mG`, `mC`, `mU` | 2'-OMe | Prefix form; same chemistry as suffix `m` |
| `Af`, `Cf`, `Gf`, `Uf`, `Tf` | 2'-Fluoro | Uppercase base + `f` suffix |
| `fA`, `fG`, `fC`, `fU` | 2'-Fluoro | Prefix form |
| `mg` | 2'-OMe G variant | Seen in some filings; treat as 2'-OMe unless patent specifies otherwise |
| `dT`, `dA`, `dC`, `dG` | 2'-deoxy (DNA) | H at 2′, no hydroxyl; `dTdT` or `TT` at 3′ end = deoxythymidine overhang |

### Backbone modifications

| Symbol pattern | Decoded type | Notes |
|---|---|---|
| `*` between residues | phosphorothioate | Linkage at that position |
| `s` suffix (`as`, `us`, `Gs`, `Cs`) | phosphorothioate | Linkage after that nucleotide, 1-indexed |
| no marker | phosphodiester | Standard unmodified backbone |

### Ambiguous / multi-use symbols

| Symbol | Default interpretation | Flag? |
|---|---|---|
| `n` | any nucleotide (degenerate position) | Flag if in a modification context |
| `I` | inosine | Flag; confirm from context |
| `T` / `t` | thymidine or ribothymidine | Confirm from context; at 3′ overhang = `dT` |

---

## Style B — case convention

Some patents encode modifications via letter case rather than suffix codes.
The exact convention is always stated in a footnote on the page or at the end
of the table. Read that footnote before applying any case rule.

Common conventions seen in Alnylam filings:

- **lowercase** = 2'-OMe ribonucleotide; **uppercase** = unmodified RNA
- **lowercase** = ribonucleotide (RNA); **uppercase** = deoxyribonucleotide (DNA)

Do not assume which convention applies. If the page has no footnote and no
`abbrev_map` entry covering case, mark the modification as `flagged: true` and
`needs_review: true`.

---

## Conjugates and other modifications

| Symbol / term | Decoded type | Notes |
|---|---|---|
| `L96` suffix | GalNAc | Hepatocyte-targeting ligand; not a sugar or backbone mod |
| `GalNAc` | GalNAc | Same; may appear as text description rather than inline code |
| `LNA` | LNA (locked nucleic acid) | Bicyclic sugar; described in text, rarely as inline code |
| `2'-MOE` / `2'-O-Methoxyethyl` | 2'-MOE | Described in text |
| `2'-NH2` / `2'-Amino` | 2'-Amino | Described in text |
| Cholesterol / `Chol` | cholesterol conjugate | End-conjugate; note strand and position (3′/5′) |

---

## Extraction rules

- **Do not invent data.** Only emit a record for a siRNA the text explicitly
  describes with at least an identifier or a sequence.
- **One `modifications[]` entry per modification type per position.** Positions
  are 1-indexed from the 5′ end of the relevant strand.
- **Ambiguous symbol:** put your best-guess decode in `type`, the original
  symbol in `raw_symbol`, set `flagged: true` on that modification entry, and
  `needs_review: true` on the record.
- **Missing data:** omit the field rather than guessing. Never synthesize a
  sequence, position, or activity value not stated in the source.
- **`needs_review: true`** whenever any field is uncertain, ambiguous, or
  absent but expected.
- **Source provenance:** populate `source_quote` with the verbatim snippet
  (table row, sentence) from which the record was drawn. This is the primary
  audit trail.
- **OCR damage:** if a sequence looks corrupted (wrong length, non-IUPAC
  characters, broken mid-run), note it in `source_quote`, set
  `needs_review: true`, and do not attempt to auto-correct. Flag as
  `ocr_damage_suspected` in the record's notes if that field exists.
- **Sequence ground truth:** always prefer the `sequence.txt` / ST.25 / ST.26
  listing over OCR'd table sequences when both are available. Join on SEQ ID NO.
  See `sequence-listing-parsing-guide.md` for ST format details.

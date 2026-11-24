---
name: patent-sequence-extraction
description: Extract patent siRNA sequences and quantitative activity from OCR, EPO descriptions and ST.25/ST.26 listings.
version: 3.1.0
---
# Patent sequence extraction
Use the seven steps below. The prepared interpreter, inputs and output paths are supplied by the launcher. Do not rediscover the environment. Use only this patent’s sources and the skill; other trials and patents are not evidence. Use the old PS7 stages: Steps 1–3, Steps 4–5, Steps 6–7. With prepared page and sequence caches, skip the first conversation entirely; do not reconfirm or audit these inputs. Within Steps 4–5, save each table and continue in the SAME conversation; never stop just because one table or batch is saved.

## Ground rule: the sequence listing is authoritative

A printed ST.25/ST.26 listing beats OCR'd body tables, always. Spot-testing Alnylam filings put OCR body sequences at only 35–41% exact agreement with the listing; the dominant failures were dropped `U` and G↔U confusion. So:

- When a printed SEQ ID maps explicitly into the listing, take the bases verbatim from the listing and set `sequence_source: "sequence_txt"`.
- Without a listing, take the printed sequence, set `sequence_source: "ocr"`, and keep the OCR notation for chemistry. An absent listing is normal, not an error.
- Never repair an ambiguous base. Flag it and move on; a listing is not proof that chemistry is resolved.

## Reading order

When a sequence listing is present, read it before the body tables: **3 → 1 → 2 → 4 → 5 → 7** (Step 6 still runs when chemistry is unresolved). Without a listing, follow Steps 1 → 7.

1. Start from the pages the prepared cache marks SEQ and ACTIVITY; that routing is the route to the tables. Read a page the cache does not pool only when a printed table or paragraph points there. Do not re-scan the whole patent or audit the cache.
2. Legend tables are routed to output/activity/step6_legends.json. Do not inspect or validate them during Steps 4–5; Step 6 handles them after data is saved.
3. Use the prepared ST.25/ST.26 seq_cache. Listing bases are authoritative when the printed SEQ ID maps explicitly. Without listing, use source sequence text and mark OCR uncertainty. Never repair ambiguous bases. EPO body text is fulltext, not sequence_txt. A listing is not proof that chemistry is resolved.
4. Read the sequence tables, including continuations and footnotes. Build schemas_pydantic.SiRNARecord objects: one strand per record, not a duplex dictionary. Consult model_json_schema() once for field names; this is the only output contract. Preserve printed identifiers and the printed notation in source_quote; `sequence` itself holds bases only (A/C/G/T/U/N) and chemistry never goes there. Column position may help identify strands only when the table supports it. Do not automatically strip suffixes or merge modified/unmodified compounds. Save reliable rows immediately through the provided commit helper; uncertain rows go to pending. Defer chemistry to Step 6.

   **Modifications are grouped, not per-nucleotide rows.** Write `modifications` as a list of groups; one group is one chemistry applied at one or more positions of one strand:

   ```json
   {"canonical_name": "Phosphorothioate", "positions": [1, 2], "raw_names": ["*m", "*f"]}
   {"canonical_name": "2-O-Methyl",       "positions": [1, 2, 3], "raw_names": ["mU", "mG", "mG"]}
   ```

   - `canonical_name` is one name from `references/modification_vocabulary.json`; that file is the only authority on spelling and aliases (162 canonical names, 1,022 aliases). Write `2-O-Methyl`, not `2'-OMe`; `2-F`, not `2'-Fluoro`.
   - `positions` are 1-based within that strand, ascending, one entry per nucleotide the chemistry covers.
   - `raw_names` are the printed symbols, one per position in the same order, so the evidence survives even when the canonical name is `other`.
   - Resolve the name from the patent's own `abbrev_map` first, then from the defaults implemented in `scripts/modification_extractor.py`. If neither supports a name, use `other` and keep `raw_names`.
   - A strand with no chemistry has `modifications: []`. Do not emit bare base letters as modifications.

5. Read measured activity even on pages that also contain sequences. Resolve IDs using explicit table/paragraph evidence. Save a separate record for each distinct measurement/condition. Percent knockdown and percent mRNA remaining are different metrics; a ratio becomes percent only when its normalization is explicit. Preserve precision, inequality/uncertainty in notes and raw evidence. Do not treat dose, assay instructions or control definitions as measured activity. Mixture activity belongs to the mixture, not automatically each component. Missing mappings are local pending items, not a reason to fail the patent. First commit one confirmed compound and measurement with final:false. Then continue incrementally; do not wait for other compounds or repeated source verification. Commit each finished table before the next.
6. Only now resolve chemistry from the patent legend and original source_quote; use scripts/modification_extractor.py and references/modification_extraction.md when relevant. Letter case is a convention, not a fact: confirm it from the legend (lowercase = 2'-OMe, or lowercase = RNA with uppercase = DNA); with no stated convention keep the record flagged instead of assuming 2'-O-Methyl. Preserve already saved sequence/activity/identity. A sequence listing does not remove this step. Unresolved chemistry remains flagged, without blocking valid activity.

   **Legend drag.** When a SEQ/ACTIVITY page uses a modifier symbol that is not in `abbrev_map`, search forward from that page for the next page whose `abbrev_footnotes` or `abbrev_tables` define it, add that page to the working set, update `abbrev_map`, and re-decode the affected records from their `source_quote`. Drag pages only for a missing legend definition.

7. Export one row per compound with schemas_pydantic.write_duplex_csv: sense and antisense side by side, modification names and positions as the aligned ` * ` strings, `0` when a strand carries no chemistry, and an empty cell when the patent only prints one strand. Columns follow the dataset layout through `knockdown`. Every row is validated by the DuplexRow model before it is written; record pending evidence separately. A written CSV is not certification of full coverage.

## Inputs
The launcher prepares OCR page JSON and seq_cache. EPO evidence is in output/activity/epo_context.json, with original paragraph IDs/order and XML path. Read only matching table/example/ID context. EPO segment indexes are not PDF page numbers. Embedded images or absent table cells require the matching OCR/PDF; never invent missing values. The raw XML stays accessible. Use the supplied PDF renderer only if text cannot resolve the local evidence.

## References
- `references/sequence_listing_parsing_guide.md`: ST.25/ST.26 formats.
- `references/modification_vocabulary.json`: canonical modification names and aliases; `canonical_name` must come from here.
- `references/patent_layout_judgement.md`: the judgements no script can make — which table this is, continuation pages, id variants, activity attribution, OCR damage shapes. Read it when a table looks ambiguous.
- `references/modification_extraction.md`: the chemistry reconciliation contract for chemistry that stays unresolved.
- `scripts/modification_extractor.py`: the default symbol conventions (residue prefixes/suffixes, linkage marks, bracketed conjugates); patent definitions still win.
- `scripts/modification_extractor.py`: modification tokenization, position assignment, canonical naming and grouped serialization.
- Tm-only requests retain the separate `references/tm_extraction.md` workflow.

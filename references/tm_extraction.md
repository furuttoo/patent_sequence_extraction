# Sequence–Tm extraction

Use this mode only when the requested output is sequence, standardized
modification/position groups, and melting temperature.

1. Run `scripts/tm_page_filter.py`. An isolated `Tm` modification token is not
   a melting-temperature hit. If the patent has no contextual Tm hit, skip it
   and do not create a CSV.
2. Prefer ST.25/ST.26 sequence-listing data when a table SEQ ID resolves.
   Otherwise retain a plain OCR nucleotide sequence verbatim. Never repair a
   sequence. Normalize OCR/LaTeX Tm headers such as `$T_m$`, `$T_{m}$`, and
   `$T_{{m}}$`. When Tm is published in a side table keyed by `SEQ ID pair` or
   `SEQ ID mero`, join it to the companion strand table on that exact key. Emit
   one validated duplex row carrying the pair's scalar Tm. Keep sense and
   antisense sequences and their modification/position groups in separate
   columns. Concatenate ordered sense mero fragments into the sense strand; do
   not infer a join when the key is missing or ambiguous.
3. Accept only a scalar Tm. Normalize it with decimal arithmetic to exactly two
   decimal places. Reject inequalities, ranges, approximations, uncertainty,
   and missing-value labels into the audit JSON.
4. Validate every final row with `SequenceTmOutputRecord`, then write only via
   `write_sequence_tm_csv`. The writer validates the complete batch and uses an
   atomic replacement.
5. The final column contract is:
   `patent_id_origin,sense_sequence,antisense_sequence,modification_sense,`
   `sense_position,modification_antisense,antisense_position,tm_value`.
6. Use `merge_sequence_tm_csv` to combine outputs. It revalidates every input,
   removes exact duplicates, sorts deterministically, and writes atomically.

The Tm extractor decodes known modified shorthand through
`scripts/modification_extractor.py`. Unknown chemistry is retained in the
rejection audit with reason `unresolved_chemistry`. Review unique unknowns,
write an auditable chemistry resolution map, and rerun with `--chemistry-map`;
never label unresolved modified notation as unmodified.

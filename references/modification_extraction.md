# Modification extraction and standardization

Use `scripts/modification_extractor.py` for both ordinary sequence extraction
and Tm-only extraction.

1. Build the patent abbreviation map first. Patent definitions override the
   default vocabulary.
2. Tokenize the OCR notation. Only nucleotide tokens increment the one-indexed
   position; a phosphorothioate marker is assigned to the preceding nucleotide.
3. Use the sequence listing as the final nucleotide sequence when its SEQ ID
   resolves, but always derive chemistry and positions from OCR notation.
4. Normalize known chemistry through
   `references/modification_vocabulary.json` and serialize groups as aligned
   ` * `-separated modification and position strings.
5. A token that resolves to nothing is written as `other` with its
   `raw_symbol` kept and `flagged=true`. It is never dropped, collapsed to `0`,
   or guessed from the sequence.
6. Unresolved chemistry is revisited in bulk, not per record: the redo pass
   (`scripts/redo_step6_entries.py`) groups a patent's unresolved entries, sends
   them together with that patent's legend, and merges every accepted verdict
   back onto all measurement rows of the same `(sirna_id, strand)`. A rejected
   verdict leaves the stored value untouched. `load_resolution_map()` is there
   when a hand-curated `{raw_name: canonical_name}` map exists, but the pipeline
   does not require one.

The unmodified representation is `modification=0, position=0`. A core sequence
may omit a two-nucleotide terminal overhang while the position data uses the
next two positions; the Pydantic output contract permits that dataset
convention while rejecting positions beyond it.

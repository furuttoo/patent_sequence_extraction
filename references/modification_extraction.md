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
5. Collect every unknown token as `UnresolvedChemistry`. Group unique raw terms
   across the patent and reason over them once using the patent definition,
   notation footnote, source row, and vocabulary.
6. Write an auditable chemistry map containing `raw_name`, `canonical_name`,
   resolution type, evidence, confidence, and `needs_review`. Rerun extraction
   with that map. Never silently convert an unknown chemistry to `other` or `0`.

The unmodified representation is `modification=0, position=0`. A core sequence
may omit a two-nucleotide terminal overhang while the position data uses the
next two positions; the Pydantic output contract permits that dataset
convention while rejecting positions beyond it.

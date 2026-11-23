"""
Parse an ST.25 PatentIn-format sequence listing file into a
{seq_id: {sequence, length, mol_type, organism, description, features}} dict.

ST.25 lines look like:
    <210> 1
    <211> 19
    <212> RNA
    <213> Artificial Sequence
    <220>
    <223> iRNA agent
    <400>
    gcuauguccg ucacguuga

Important pitfalls handled (see references/alnylam_pct_anatomy.md P18):
- FastSEQ appends right-side line numbers to each sequence line. Strip them
  per line BEFORE joining. Patterns to handle:
    1. Pure-numeric lines (column markers like "1 5 10 15 20")
    2. Trailing digits on the same line as the sequence (e.g. "...guuga 19")
- Some entries are 21 nt long (19 + TT 3' overhang baked in for modified strands).
  Leave as-is — the caller decides whether to strip TT for comparison.

Usage as a library:
    from parse_sequence_listing_st25 import parse_sequence_listing
    summary, records = parse_sequence_listing("WO2012078967_sequence.txt")
    # records is a dict: {seq_id_int: {sequence, length, mol_type, ...}}

Usage as a CLI:
    python parse_sequence_listing_st25.py WO2012078967_sequence.txt
"""
import re
import sys
import json
from pathlib import Path


def parse_sequence_listing(path):
    """Parse an ST.25 .txt sequence listing.

    Returns:
        summary: dict with counts (n_records, n_rna, n_dna, n_protein, n_other)
        records: dict[int, dict] keyed by SEQ ID NO
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        raw = f.read()

    # Split on <210> block boundaries; keep only chunks starting with <210>
    parts = re.split(r"(?=<210>\s*\d+\s*\n)", raw)

    records = {}
    for chunk in parts:
        chunk = chunk.strip()
        if not chunk.startswith("<210>"):
            continue

        seq_id = None
        length = None
        mol_type = None
        organism = None
        description = None
        features = []
        in_seq = False
        seq_lines = []

        for line in chunk.split("\n"):
            line = line.rstrip("\r")

            if line.startswith("<210>"):
                m = re.match(r"<210>\s*(\d+)", line)
                if m:
                    seq_id = int(m.group(1))
                in_seq = False
            elif line.startswith("<211>"):
                m = re.match(r"<211>\s*(\d+)", line)
                if m:
                    length = int(m.group(1))
            elif line.startswith("<212>"):
                m = re.match(r"<212>\s*(\S+)", line)
                if m:
                    mol_type = m.group(1)
            elif line.startswith("<213>"):
                m = re.match(r"<213>\s*(.+)", line)
                if m:
                    organism = m.group(1).strip()
            elif line.startswith("<220>"):
                pass
            elif line.startswith("<223>"):
                m = re.match(r"<223>\s*(.+)", line)
                if m:
                    description = m.group(1).strip()
            elif line.startswith("<221>"):
                m = re.match(r"<221>\s*(.+)", line)
                if m:
                    features.append(("name", m.group(1).strip()))
            elif line.startswith("<222>"):
                m = re.match(r"<222>\s*(.+)", line)
                if m:
                    features.append(("location", m.group(1).strip()))
            elif line.startswith("<400>"):
                in_seq = True
            elif in_seq:
                stripped = line.strip()
                if not stripped:
                    continue
                # P18 fix: skip pure-numeric lines (FastSEQ column markers)
                if re.match(r"^\d+$", stripped):
                    continue
                if re.match(r"^[\d\s]+$", stripped):
                    continue
                # P18 fix: strip trailing digits on the same line as the sequence
                stripped = re.sub(r"\s*\d+\s*$", "", stripped)
                seq_lines.append(stripped)

        if seq_id is not None:
            sequence = "".join(seq_lines).replace(" ", "").replace("\t", "").upper()
            # P18 final safety net
            sequence = re.sub(r"\d+$", "", sequence)
            records[seq_id] = {
                "seq_id": seq_id,
                "length": length,
                "mol_type": mol_type,
                "organism": organism,
                "description": description,
                "features": features,
                "sequence": sequence,
            }

    # Build summary
    from collections import Counter
    by_mol = Counter(r["mol_type"] for r in records.values() if r["mol_type"])
    summary = {
        "n_records": len(records),
        "n_rna": by_mol.get("RNA", 0),
        "n_dna": by_mol.get("DNA", 0),
        "n_protein": by_mol.get("PRT", 0),
        "n_other": sum(v for k, v in by_mol.items() if k not in ("RNA", "DNA", "PRT")),
        "by_mol_type": dict(by_mol),
    }
    return summary, records


def build_cache(path, force=False):
    """Parse sequence listing and write a JSON cache next to it.

    Returns: (cache_path, status) where status is "cache_built" or "cache_hit".
    """
    src = Path(path)
    cache_path = src.with_name(f"{src.stem}_seqcache.json")
    if cache_path.exists() and not force:
        return str(cache_path), "cache_hit"
    summary, records = parse_sequence_listing(path)
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump({"summary": summary, "records": records}, f, indent=2)
    return str(cache_path), "cache_built"


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python parse_sequence_listing_st25.py <sequence.txt> [--force]")
        sys.exit(1)
    force = "--force" in sys.argv
    path = sys.argv[1]
    cache_path, status = build_cache(path, force=force)
    summary, _ = parse_sequence_listing(path)
    print(f"Status: {status}")
    print(f"Cache:  {cache_path}")
    print(f"Records: {summary['n_records']}  (RNA={summary['n_rna']}  DNA={summary['n_dna']}  PRT={summary['n_protein']}  other={summary['n_other']})")

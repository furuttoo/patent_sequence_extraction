"""
parse_sequence_listing_st26.py

Parser + cache builder for ST.26-format patent sequence listings — XML
files (WIPO DTD ST26SequenceListing, post-2022 filings). Use this file
when the source is a `.xml` sequence listing with <SequenceData
sequenceIDNumber="N"> / <INSDSeq> structure.

For pre-2022 filings that use the flat-tag ST.25 `.txt` format instead,
use `parse_sequence_listing_st25.py` — do not extend this file to handle
that format.

Structural differences vs. ST.25 that this parser is built around:
  - SEQ ID lives in the <SequenceData sequenceIDNumber="N"> attribute —
    a single, schema-enforced source of truth. There is no <210>/<400>
    redundancy to cross-check here, unlike ST.25, because the XML schema
    itself guarantees a SequenceData block can't declare a mismatched ID.
  - Every position-level modification is its OWN <INSDFeature> element
    with key=modified_base and a <location> + <note> qualifier — no
    range-compressed annotations like ST.25's "<222> 1..19" covering many
    positions under one description. Expect one feature per modified base.
  - <INSDFeature key=misc_feature> blocks separately mark which position
    ranges are RNA vs. DNA within a single strand — relevant for chimeric
    backbones (e.g. GNA-containing oligos), and has no ST.25 equivalent.
  - Organism lives inside the `source` feature's qualifiers, not a
    top-level tag.
  - No OCR risk: XML sequence listings are machine-generated and
    schema-validated at submission time, so there's no character-level
    corruption to guard against the way there is with ST.25 `.txt` files
    (which were often typed/edited by hand and preceded structured
    validation tooling).

Output record schema matches parse_sequence_listing_st25.py's schema
(same field names) so downstream code doesn't need to branch on format,
only the two ST.25-specific fields (`parse_anomaly` from the <210>/<400>
cross-check) don't apply here and are simply absent.

Usage (single file):
    python parse_sequence_listing_st26.py path/to/SL.xml

Usage (cache a whole tree of patent folders):
    python parse_sequence_listing_st26.py --scan /path/to/patents_root

Caching: writes `<stem>_seqcache.json` next to each source file. Skips
re-parsing if the cache is newer than the source (mtime check). --force
rebuilds regardless.
"""

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

# ---- tunables (kept identical to the ST.25 parser for consistent roles) --

OLIGO_LEN_MIN = 15
OLIGO_LEN_MAX = 30

PRIMER_HINTS = ("primer", "rt-pcr", "pcr primer")
REFERENCE_HINTS = ("mrna", "genomic", "transcript", "reference sequence")
# ST.26 commonly uses "synthetic construct" where ST.25 used "artificial
# sequence" — both mean the same thing (a designed/synthesized molecule,
# not a natural isolate), so both are treated as the "designed oligo"
# organism signal.
ARTIFICIAL_ORGANISM_HINTS = ("artificial", "synthetic construct")


def _parse_location(raw):
    """Normalize an <INSDFeature_location> string into a list of
    (start, end) int tuples. Same syntaxes as ST.25's <222>: '1', '1..21'."""
    if not raw:
        return []
    raw = raw.strip()
    out = []
    for piece in raw.split(","):
        piece = piece.strip()
        nums = re.findall(r"\d+", piece)
        if not nums:
            continue
        if len(nums) >= 2:
            out.append((int(nums[0]), int(nums[1])))
        else:
            out.append((int(nums[0]), int(nums[0])))
    return out


def _guess_role(mol_type, organism, length, features, sequence):
    """Soft hint only — never used to drop or alter a record. Mirrors the
    ST.25 classifier but recognizes ST.26's "synthetic construct" organism
    label in addition to "artificial sequence"."""
    text_blob = " ".join(
        (f.get("qualifier") or "") for f in features
    ).lower()
    organism_l = (organism or "").lower()

    if any(h in text_blob for h in PRIMER_HINTS):
        return "primer"
    if any(h in text_blob for h in REFERENCE_HINTS) or length > OLIGO_LEN_MAX * 4:
        return "reference_sequence"
    if any(h in organism_l for h in ARTIFICIAL_ORGANISM_HINTS) and OLIGO_LEN_MIN <= length <= OLIGO_LEN_MAX:
        return "oligo_candidate"
    if length > OLIGO_LEN_MAX:
        return "reference_sequence"
    return "unclassified"


def parse_record(sd_element):
    """Parse one <SequenceData> element into the shared record schema."""
    sid_raw = sd_element.get("sequenceIDNumber")
    if sid_raw is None:
        return None
    sid = int(sid_raw)

    insd = sd_element.find("INSDSeq")
    if insd is None:
        return {"seq_id": sid, "parse_anomaly": "no_INSDSeq_element"}

    length = int(insd.findtext("INSDSeq_length") or 0)
    mol_type = insd.findtext("INSDSeq_moltype", "")
    sequence = (insd.findtext("INSDSeq_sequence") or "").upper()

    # organism lives inside the 'source' feature's qualifiers
    organism = ""
    for feat in insd.findall(".//INSDFeature"):
        if feat.findtext("INSDFeature_key") == "source":
            for q in feat.findall(".//INSDQualifier"):
                if q.findtext("INSDQualifier_name") == "organism":
                    organism = q.findtext("INSDQualifier_value", "")

    # extract every non-source feature: modified_base (per position) and
    # misc_feature (RNA/DNA span annotations)
    features = []
    for feat in insd.findall(".//INSDFeature"):
        key = feat.findtext("INSDFeature_key", "")
        if key == "source":
            continue

        loc_raw = feat.findtext("INSDFeature_location", "")
        location = _parse_location(loc_raw)

        quals = {}
        for q in feat.findall(".//INSDQualifier"):
            qname = q.findtext("INSDQualifier_name", "")
            qval = q.findtext("INSDQualifier_value", "")
            quals[qname] = qval

        qualifier = quals.get("note") or quals.get("mod_base") or ""

        features.append({
            "key": key,          # "modified_base" or "misc_feature"
            "location": location,
            "qualifier": qualifier,
            # ST.26-specific: RNA/DNA span info from misc_feature, no
            # ST.25 equivalent — kept as a separate field rather than
            # overloading `qualifier` so downstream code can tell them apart
            "misc_note": quals.get("note") if key == "misc_feature" else None,
        })

    role = _guess_role(mol_type, organism, length, features, sequence)

    record = {
        "seq_id": sid,
        "length_declared": length,
        "length_actual": len(sequence),
        "mol_type": mol_type,
        "organism": organism,
        "features": features,
        "sequence": sequence,
        "likely_role": role,
        "format": "ST26",
    }
    if length != len(sequence):
        record["length_mismatch"] = True
    return record


def parse_sequence_listing(path):
    tree = ET.parse(str(path))
    root = tree.getroot()

    declared_total_str = root.findtext("SequenceTotalQuantity")
    declared_total = int(declared_total_str) if declared_total_str else None

    records = {}
    parse_errors = []
    duplicate_ids = []
    anomalies = []

    for sd in root.findall("SequenceData"):
        try:
            rec = parse_record(sd)
        except Exception as e:
            parse_errors.append(f"seq {sd.get('sequenceIDNumber', '?')}: {e}")
            continue
        if rec is None:
            continue
        sid = rec.get("seq_id")
        if sid in records:
            duplicate_ids.append(sid)
        if rec.get("parse_anomaly"):
            anomalies.append(f"seq_id {sid}: {rec['parse_anomaly']}")
        records[sid] = rec

    gaps = []
    if records:
        present = set(records.keys())
        expected_max = declared_total if declared_total else max(present)
        gaps = sorted(set(range(1, expected_max + 1)) - present)

    summary = {
        "source_file": str(path),
        "format": "ST26",
        "declared_total_sequences": declared_total,
        "parsed_total_sequences": len(records),
        "patentin_version": None,  # not applicable to ST.26
        "parse_errors": parse_errors,
        "count_mismatch": (
            declared_total is not None and declared_total != len(records)
        ),
        "duplicate_seq_ids": duplicate_ids,
        "missing_seq_ids": gaps,
        # <210>/<400> cross-check has no ST.26 equivalent — the schema
        # enforces ID integrity structurally, so this stays empty here
        "id_anomalies": anomalies,
    }
    return summary, records


def build_cache(xml_path, force=False):
    xml_path = Path(xml_path)
    cache_path = xml_path.parent / f"{xml_path.stem}_seqcache.json"

    if cache_path.exists() and not force:
        src_mtime = xml_path.stat().st_mtime
        cache_mtime = cache_path.stat().st_mtime
        if cache_mtime >= src_mtime:
            return cache_path, "cache_hit"

    summary, records = parse_sequence_listing(xml_path)
    payload = {"summary": summary, "records": records}
    cache_path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return cache_path, "cache_built"


def scan_and_cache(root_dir, force=False):
    root = Path(root_dir)
    results = []
    for xml_file in sorted(root.rglob("*.xml")):
        if "seqcache" in xml_file.stem.lower():
            continue
        cache_path, status = build_cache(xml_file, force=force)
        summary = json.loads(cache_path.read_text())["summary"]
        flags = []
        if summary["count_mismatch"]:
            flags.append("COUNT_MISMATCH")
        if summary["duplicate_seq_ids"]:
            flags.append(f"DUPLICATES({len(summary['duplicate_seq_ids'])})")
        if summary["missing_seq_ids"]:
            flags.append(f"GAPS({len(summary['missing_seq_ids'])})")
        if summary["id_anomalies"]:
            flags.append(f"ANOMALIES({len(summary['id_anomalies'])})")
        flagstr = "  <-- " + ", ".join(flags) if flags else ""
        print(f"{status:12s} {summary['parsed_total_sequences']:>5}/"
              f"{str(summary['declared_total_sequences']):>5}  {xml_file}{flagstr}")
        results.append({
            "xml_file": str(xml_file),
            "cache": str(cache_path),
            "status": status,
            "declared_total_sequences": summary["declared_total_sequences"],
            "parsed_total_sequences": summary["parsed_total_sequences"],
            "count_mismatch": summary["count_mismatch"],
        })
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", nargs="?", help="single ST.26 .xml sequence listing to parse")
    ap.add_argument("--scan", metavar="ROOT_DIR",
                     help="scan a directory tree and cache every *.xml sequence listing found")
    ap.add_argument("--force", action="store_true",
                     help="rebuild cache even if a newer one already exists")
    args = ap.parse_args()

    if args.scan:
        results = scan_and_cache(args.scan, force=args.force)
        print(f"\n{len(results)} .xml file(s) cached under {args.scan}")
        return

    if not args.path:
        ap.error("provide a path or --scan ROOT_DIR")

    cache_path, status = build_cache(args.path, force=args.force)
    summary = json.loads(cache_path.read_text())["summary"]
    print(json.dumps(summary, indent=2))
    print(f"\n[{status}] cache written to: {cache_path}")


if __name__ == "__main__":
    main()

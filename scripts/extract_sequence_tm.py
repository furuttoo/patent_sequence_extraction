"""Extract merge-safe sequence/Tm rows from OCR patent JSON.

The extractor handles plain sequences, dataset-standard modification/position
columns, and known chemically encoded notation through the shared modification
extractor. Unknown chemistry is retained in the rejection audit until an
agent-authored resolution map is supplied; it is never mislabeled unmodified.
"""

from __future__ import annotations

import argparse
import html
from html.parser import HTMLParser
import importlib.util
import json
import os
import re
from pathlib import Path
import sys
from typing import Any, Mapping, Optional

# macOS' system Python does not include this skill's Pydantic dependency. When
# an agent invokes this script with plain ``python3``, transparently restart it
# under the Hermes environment before importing any Pydantic-backed modules.
if importlib.util.find_spec("pydantic") is None:
    hermes_python = Path(
        os.environ.get(
            "HERMES_PYTHON",
            Path.home() / ".hermes/hermes-agent/venv/bin/python",
        )
    ).expanduser()
    if hermes_python.is_file() and Path(sys.executable).resolve() != hermes_python.resolve():
        os.execv(
            str(hermes_python),
            [str(hermes_python), str(Path(__file__).resolve()), *sys.argv[1:]],
        )
    raise ModuleNotFoundError(
        "Pydantic is required; set HERMES_PYTHON to a Python environment containing pydantic"
    )

from pydantic import ValidationError

from modification_extractor import extract_modifications, load_resolution_map
from schemas_pydantic import (
    ModificationGroup,
    SequenceTmOutputRecord,
    serialize_modification_groups,
    write_sequence_tm_csv,
)
from tm_page_filter import build_cache


_SEQ_HEADER = re.compile(r"sequence|oligonucleotide|oligo", re.I)
_SEQ_ID_HEADER = re.compile(r"SEQ\s*ID", re.I)
_PAIR_HEADER = re.compile(r"SEQ\s*ID.*(?:pair|mero)|(?:pair|mero).*SEQ\s*ID", re.I)
_MOD_HEADER = re.compile(r"modification", re.I)
_POSITION_HEADER = re.compile(r"position", re.I)
_PLAIN_SEQUENCE = re.compile(r"[ACGTUNacgtun\s]+")
_PAIR_VALUE = re.compile(r"(?P<pair>\d+(?:\s*/\s*\d+)+)")
_SCALAR_PREFIXED_PAIR = re.compile(
    r"^\s*(?P<tm>[+-]?(?:\d+(?:[.,]\d+)?|\.\d+))\s+"
    r"(?P<pair>\d+(?:\s*/\s*\d+)+)\s*$"
)


class _TableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.rows: list[list[str]] = []
        self._row: Optional[list[str]] = None
        self._cell: Optional[list[str]] = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, Optional[str]]]) -> None:
        if tag.lower() == "tr":
            self._row = []
        elif tag.lower() in {"td", "th"} and self._row is not None:
            self._cell = []

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._cell.append(data)

    def handle_endtag(self, tag: str) -> None:
        lower = tag.lower()
        if lower in {"td", "th"} and self._cell is not None and self._row is not None:
            self._row.append(re.sub(r"\s+", " ", html.unescape("".join(self._cell))).strip())
            self._cell = None
        elif lower == "tr" and self._row is not None:
            if any(self._row):
                self.rows.append(self._row)
            self._row = None


def parse_html_table(value: str) -> list[list[str]]:
    parser = _TableParser()
    parser.feed(value)
    return parser.rows


def _patent_id(cache: Mapping[str, Any]) -> str:
    return str(cache.get("patent_id") or "")


def _compact_header(value: str) -> str:
    """Flatten OCR/LaTeX header markup for semantic matching."""
    text = html.unescape(str(value or "")).replace("\\n", " ")
    return re.sub(r"[^A-Za-z0-9]+", "", text).casefold()


def _is_tm_header(value: str) -> bool:
    compact = _compact_header(value)
    if "meltingtemperature" in compact:
        return True
    if compact in {"tm", "tmc", "tmoc", "tmdegc"}:
        return True
    # OCR commonly preserves TeX punctuation in headers such as $T_m$,
    # $T_{m}$, $T_{{m}}$, or appends T_m to a merged sequence header.
    return ("$" in value or "_" in value or "{" in value) and compact.endswith("tm")


def _pair_key(value: str) -> Optional[str]:
    match = _PAIR_VALUE.search(str(value or ""))
    if not match:
        return None
    return "/".join(re.findall(r"\d+", match.group("pair")))


def _split_sequence_fragments(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"(?:\\n|\n)+", value or "") if part.strip()]


def _header_indexes(rows: list[list[str]]) -> Optional[tuple[int, dict[str, int]]]:
    for row_index, row in enumerate(rows[:5]):
        indexes: dict[str, int] = {}
        for index, cell in enumerate(row):
            if _SEQ_ID_HEADER.search(cell) and not _PAIR_HEADER.search(cell):
                indexes.setdefault("seq_id", index)
            elif _SEQ_HEADER.search(cell):
                indexes.setdefault("sequence", index)
            if _is_tm_header(cell):
                indexes.setdefault("tm", index)
            if _MOD_HEADER.search(cell):
                indexes.setdefault("modification", index)
            if _POSITION_HEADER.search(cell):
                indexes.setdefault("position", index)
        if "sequence" in indexes and "tm" in indexes:
            return row_index, indexes
    return None


def _pair_table_header(
    rows: list[list[str]],
) -> Optional[tuple[int, int, Optional[int], list[int]]]:
    """Return header row, pair column, optional Tm column, and sequence columns."""
    for row_index, row in enumerate(rows[:6]):
        pair_columns = [i for i, cell in enumerate(row) if _PAIR_HEADER.search(cell)]
        if not pair_columns:
            continue
        tm_columns = [i for i, cell in enumerate(row) if _is_tm_header(cell)]
        sequence_columns = [
            i for i, cell in enumerate(row)
            if _SEQ_HEADER.search(cell) and not _PAIR_HEADER.search(cell)
        ]
        return row_index, pair_columns[-1], tm_columns[0] if tm_columns else None, sequence_columns
    return None


def _build_pair_sequence_map(
    cache: Mapping[str, Any],
) -> dict[str, list[tuple[str, str]]]:
    """Index OCR strand notation by a slash-delimited SEQ ID pair/mero key."""
    pairs: dict[str, list[tuple[str, str]]] = {}
    for page in cache.get("pages", []):
        for block in page.get("blocks", []):
            if block.get("block_label") != "table":
                continue
            rows = parse_html_table(str(block.get("block_content", "")))
            header = _pair_table_header(rows)
            if header is None:
                continue
            header_index, pair_column, _tm_column, sequence_columns = header
            if not sequence_columns:
                continue
            data_rows = rows[header_index + 1 :]
            for offset, row in enumerate(data_rows):
                if pair_column >= len(row):
                    continue
                key = _pair_key(row[pair_column])
                if key is None:
                    continue
                seq_ids = key.split("/")
                raw_sequences = [
                    fragment
                    for column in sequence_columns
                    if column < len(row)
                    for fragment in _split_sequence_fragments(row[column])
                ]
                # Some OCR tables put a third mero strand on the following
                # short row as ``SEQ ID | sequence``.
                lookahead = offset + 1
                while len(raw_sequences) < len(seq_ids) and lookahead < len(data_rows):
                    continuation = data_rows[lookahead]
                    if len(continuation) >= 2 and continuation[0].strip() == seq_ids[len(raw_sequences)]:
                        raw_sequences.extend(_split_sequence_fragments(continuation[1]))
                        lookahead += 1
                        continue
                    break
                if len(raw_sequences) < len(seq_ids):
                    continue
                pairs[key] = list(zip(seq_ids, raw_sequences[: len(seq_ids)]))
    return pairs


def _pair_tm_value(row: list[str], tm_column: int, pair_column: int) -> Optional[str]:
    if tm_column < len(row):
        candidate = row[tm_column].strip()
        if re.fullmatch(r"[+-]?(?:\d+(?:[.,]\d+)?|\.\d+)", candidate):
            return candidate
    # Column boundaries sometimes collapse so the last cell becomes
    # ``46.3 50/65/82`` while T_m remains attached to the sequence header.
    if pair_column < len(row):
        merged = _SCALAR_PREFIXED_PAIR.fullmatch(row[pair_column])
        if merged:
            return merged.group("tm")
    return row[tm_column] if tm_column < len(row) else None


def _clean_plain_sequence(raw: str) -> Optional[str]:
    if not raw or not _PLAIN_SEQUENCE.fullmatch(raw):
        return None
    sequence = re.sub(r"\s+", "", raw)
    # Mixed case is a modification signal in supported patent notation.
    # Uniform lowercase is allowed as an unmodified OCR rendering.
    if re.search(r"[ACGUTN]", sequence) and re.search(r"[acgutn]", sequence):
        return None
    return sequence.upper() or None


def _lookup_sequence(seq_cache: Mapping[Any, Any], raw_seq_id: str) -> Optional[str]:
    match = re.search(r"\d+", raw_seq_id or "")
    if not match:
        return None
    seq_id = int(match.group(0))
    record = seq_cache.get(seq_id) or seq_cache.get(str(seq_id))
    if isinstance(record, Mapping):
        sequence = record.get("sequence")
        return str(sequence).upper() if sequence else None
    return None


def load_sequence_listing(path: str | Path | None) -> dict[Any, Any]:
    if path is None:
        return {}
    listing = Path(path)
    if listing.suffix.lower() == ".xml":
        from parse_sequence_listing_st26 import parse_sequence_listing
    else:
        from parse_sequence_listing_st25 import parse_sequence_listing
    _, records = parse_sequence_listing(listing)
    return records


def _append_sequence_tm_record(
    records: list[SequenceTmOutputRecord],
    rejected: list[dict[str, Any]],
    *,
    cache: Mapping[str, Any],
    seq_cache: Mapping[Any, Any],
    chemistry_resolutions: Mapping[str, str],
    page_index: Any,
    source_row_index: int,
    source_row: list[str],
    raw_sequence: str,
    raw_seq_id: str = "",
    tm_value: Any,
    explicit_modification: Optional[str] = None,
    explicit_position: Optional[str] = None,
) -> None:
    ocr_sequence = _clean_plain_sequence(raw_sequence)
    has_explicit_modification_columns = (
        explicit_modification is not None and explicit_position is not None
    )
    decoded_modification = None
    if ocr_sequence is None:
        decoded_modification = extract_modifications(
            raw_sequence,
            abbrev_map=cache.get("abbrev_map", {}),
            resolutions=chemistry_resolutions,
        )
        if decoded_modification.unresolved and not has_explicit_modification_columns:
            rejected.append(
                {
                    "page_index": page_index,
                    "row_index": source_row_index,
                    "source_row": source_row,
                    "raw_sequence": raw_sequence,
                    "seq_id": raw_seq_id or None,
                    "reason": "unresolved_chemistry",
                    "unresolved": [
                        item.model_dump() for item in decoded_modification.unresolved
                    ],
                }
            )
            return
        ocr_sequence = decoded_modification.base_sequence

    sequence = _lookup_sequence(seq_cache, raw_seq_id) if raw_seq_id else None
    sequence = sequence or ocr_sequence
    if sequence is None:
        rejected.append(
            {
                "page_index": page_index,
                "row_index": source_row_index,
                "source_row": source_row,
                "raw_sequence": raw_sequence,
                "seq_id": raw_seq_id or None,
                "reason": "sequence_decode_failed",
            }
        )
        return

    if decoded_modification is not None and not has_explicit_modification_columns:
        modification, position = serialize_modification_groups(decoded_modification.groups)
    else:
        modification = explicit_modification or "0"
        position = explicit_position or "0"
    payload = {
        "patent_id_origin": _patent_id(cache),
        "sequence": sequence,
        "modification": modification,
        "position": position,
        "tm_value": tm_value,
    }
    try:
        records.append(SequenceTmOutputRecord.model_validate(payload))
    except ValidationError as exc:
        rejected.append(
            {
                "page_index": page_index,
                "row_index": source_row_index,
                "source_row": source_row,
                "raw_sequence": raw_sequence,
                "seq_id": raw_seq_id or None,
                "reason": "pydantic_validation_failed",
                "error": str(exc),
            }
        )


def _append_duplex_tm_record(
    records: list[SequenceTmOutputRecord],
    rejected: list[dict[str, Any]],
    *,
    cache: Mapping[str, Any],
    seq_cache: Mapping[Any, Any],
    chemistry_resolutions: Mapping[str, str],
    page_index: Any,
    source_row_index: int,
    source_row: list[str],
    components: list[tuple[str, str]],
    tm_value: Any,
) -> None:
    """Build one duplex row: first component antisense, remaining components sense."""
    if len(components) < 2:
        rejected.append({
            "page_index": page_index, "row_index": source_row_index,
            "source_row": source_row, "reason": "missing_duplex_partner",
        })
        return

    decoded_components: list[tuple[str, list[ModificationGroup]]] = []
    unresolved_items = []
    for seq_id, raw_sequence in components:
        plain = _clean_plain_sequence(raw_sequence)
        decoded = None
        if plain is None:
            decoded = extract_modifications(
                raw_sequence,
                abbrev_map=cache.get("abbrev_map", {}),
                resolutions=chemistry_resolutions,
            )
            unresolved_items.extend(decoded.unresolved)
            plain = decoded.base_sequence
        sequence = _lookup_sequence(seq_cache, seq_id) or plain
        if not sequence:
            rejected.append({
                "page_index": page_index, "row_index": source_row_index,
                "source_row": source_row, "seq_id": seq_id,
                "reason": "sequence_decode_failed",
            })
            return
        decoded_components.append((sequence, decoded.groups if decoded else []))

    if unresolved_items:
        rejected.append({
            "page_index": page_index, "row_index": source_row_index,
            "source_row": source_row, "reason": "unresolved_chemistry",
            "unresolved": [item.model_dump() for item in unresolved_items],
        })
        return

    def combine(parts: list[tuple[str, list[ModificationGroup]]]) -> tuple[str, str, str]:
        sequence = "".join(part[0] for part in parts)
        grouped: dict[str, set[int]] = {}
        offset = 0
        for part_sequence, groups in parts:
            for group in groups:
                grouped.setdefault(group.canonical_name, set()).update(
                    offset + position for position in group.positions
                )
            offset += len(part_sequence)
        merged_groups = [
            ModificationGroup(canonical_name=name, positions=tuple(sorted(positions)))
            for name, positions in grouped.items()
        ]
        modification, position = serialize_modification_groups(merged_groups)
        return sequence, modification, position

    antisense_sequence, modification_antisense, antisense_position = combine(
        decoded_components[:1]
    )
    sense_sequence, modification_sense, sense_position = combine(decoded_components[1:])
    payload = {
        "patent_id_origin": _patent_id(cache),
        "sense_sequence": sense_sequence,
        "antisense_sequence": antisense_sequence,
        "modification_sense": modification_sense,
        "sense_position": sense_position,
        "modification_antisense": modification_antisense,
        "antisense_position": antisense_position,
        "tm_value": tm_value,
    }
    try:
        records.append(SequenceTmOutputRecord.model_validate(payload))
    except ValidationError as exc:
        rejected.append({
            "page_index": page_index, "row_index": source_row_index,
            "source_row": source_row, "reason": "pydantic_validation_failed",
            "error": str(exc),
        })


def extract_from_cache(
    cache: Mapping[str, Any],
    *,
    seq_cache: Optional[Mapping[Any, Any]] = None,
    chemistry_resolutions: Optional[Mapping[str, str]] = None,
) -> tuple[list[SequenceTmOutputRecord], list[dict[str, Any]]]:
    seq_cache = seq_cache or {}
    chemistry_resolutions = chemistry_resolutions or {}
    records: list[SequenceTmOutputRecord] = []
    rejected: list[dict[str, Any]] = []
    pair_sequences = _build_pair_sequence_map(cache)

    for page in cache.get("pages", []):
        page_index = page.get("page_index")
        for block in page.get("blocks", []):
            if block.get("block_label") != "table":
                continue
            rows = parse_html_table(str(block.get("block_content", "")))
            pair_header = _pair_table_header(rows)
            if pair_header is not None and pair_header[2] is not None:
                header_index, pair_column, tm_column, _sequence_columns = pair_header
                required_max = max(pair_column, tm_column)
                for source_row_index, row in enumerate(
                    rows[header_index + 1 :], start=header_index + 1
                ):
                    if len(row) <= required_max:
                        continue
                    key = _pair_key(row[pair_column])
                    if key is None:
                        continue
                    components = pair_sequences.get(key)
                    if not components:
                        rejected.append(
                            {
                                "page_index": page_index,
                                "row_index": source_row_index,
                                "source_row": row,
                                "seq_id_pair": key,
                                "reason": "pair_sequence_join_failed",
                            }
                        )
                        continue
                    tm_value = _pair_tm_value(row, tm_column, pair_column)
                    _append_duplex_tm_record(
                        records,
                        rejected,
                        cache=cache,
                        seq_cache=seq_cache,
                        chemistry_resolutions=chemistry_resolutions,
                        page_index=page_index,
                        source_row_index=source_row_index,
                        source_row=row,
                        components=components,
                        tm_value=tm_value,
                    )
                continue

            header = _header_indexes(rows)
            if header is None:
                continue
            header_index, indexes = header
            required_max = max(indexes.values())

            for source_row_index, row in enumerate(rows[header_index + 1 :], start=header_index + 1):
                if len(row) <= required_max:
                    rejected.append(
                        {
                            "page_index": page_index,
                            "row_index": source_row_index,
                            "source_row": row,
                            "reason": "short_table_row",
                        }
                    )
                    continue
                rejected.append({
                    "page_index": page_index,
                    "row_index": source_row_index,
                    "source_row": row,
                    "reason": "missing_duplex_partner",
                })
    return records, rejected


def aggregate_unresolved_chemistry(rejected: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group unresolved chemistry so an agent reasons about each raw term once."""
    grouped: dict[str, dict[str, Any]] = {}
    for rejected_row in rejected:
        for item in rejected_row.get("unresolved", []):
            raw_name = str(item.get("raw_name") or item.get("raw_symbol") or "").strip()
            if not raw_name:
                continue
            entry = grouped.setdefault(
                raw_name,
                {
                    "raw_name": raw_name,
                    "raw_symbols": set(),
                    "patent_definitions": set(),
                    "positions": set(),
                    "occurrences": 0,
                    "examples": [],
                },
            )
            if item.get("raw_symbol"):
                entry["raw_symbols"].add(item["raw_symbol"])
            if item.get("patent_definition"):
                entry["patent_definitions"].add(item["patent_definition"])
            entry["positions"].update(item.get("positions", []))
            entry["occurrences"] += 1
            if len(entry["examples"]) < 3:
                entry["examples"].append(
                    {
                        "sequence_ocr": item.get("sequence_ocr"),
                        "page_index": rejected_row.get("page_index"),
                        "source_row": rejected_row.get("source_row"),
                    }
                )

    return [
        {
            "raw_name": entry["raw_name"],
            "raw_symbols": sorted(entry["raw_symbols"]),
            "patent_definitions": sorted(entry["patent_definitions"]),
            "positions": sorted(entry["positions"]),
            "occurrences": entry["occurrences"],
            "examples": entry["examples"],
        }
        for entry in sorted(grouped.values(), key=lambda item: item["raw_name"].casefold())
    ]


def extract_patent(
    patent_json: str | Path,
    output_csv: str | Path,
    *,
    sequence_listing: str | Path | None = None,
    chemistry_map: str | Path | None = None,
    rejected_json: str | Path | None = None,
) -> dict[str, Any]:
    cache = build_cache(patent_json)
    if cache["skipped_no_tm"]:
        return {"status": "skipped_no_tm", "patent_id": cache["patent_id"], "rows": 0}

    records, rejected = extract_from_cache(
        cache,
        seq_cache=load_sequence_listing(sequence_listing),
        chemistry_resolutions=load_resolution_map(chemistry_map),
    )
    write_sequence_tm_csv(records, output_csv)

    unresolved = aggregate_unresolved_chemistry(rejected)
    if rejected_json is not None:
        destination = Path(rejected_json)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with open(destination, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "patent_id_origin": cache["patent_id"],
                    "unresolved_chemistry": unresolved,
                    "rejected_rows": rejected,
                },
                handle,
                ensure_ascii=False,
                indent=2,
            )

    return {
        "status": "written",
        "patent_id": cache["patent_id"],
        "rows": len(records),
        "rejected": len(rejected),
        "unresolved_chemistry": len(unresolved),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patent_json")
    parser.add_argument("output_csv")
    parser.add_argument("--sequence-listing")
    parser.add_argument("--chemistry-map")
    parser.add_argument("--rejected-json")
    args = parser.parse_args()
    result = extract_patent(
        args.patent_json,
        args.output_csv,
        sequence_listing=args.sequence_listing,
        chemistry_map=args.chemistry_map,
        rejected_json=args.rejected_json,
    )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()

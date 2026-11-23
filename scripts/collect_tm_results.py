#!/usr/bin/env python3
"""Collect completed patent Tm CSVs into one validated, merge-safe dataset.

The extractor writes each final CSV atomically, so this collector can safely run
while extraction workers are active. Files that are absent, temporary, malformed,
or fail SequenceTmOutputRecord validation are never included.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Any


LEGACY_TM_CSV_COLUMNS = [
    "patent_id_origin",
    "sequence",
    "modification",
    "position",
    "tm_value",
]

_BROAD_TM_EVIDENCE = re.compile(
    r"melting|denaturation|(?<![A-Za-z0-9])T\s*\$?\s*(?:[_-]\s*\{*\s*)?m\s*\}*\s*\$?(?![A-Za-z0-9])",
    re.I,
)


def ensure_pydantic_runtime() -> None:
    """Re-exec with Hermes' Python when the caller's Python lacks Pydantic."""
    if importlib.util.find_spec("pydantic") is not None:
        return
    hermes_python = Path.home() / ".hermes/hermes-agent/venv/bin/python"
    if hermes_python.is_file() and Path(sys.executable).resolve() != hermes_python.resolve():
        os.execv(str(hermes_python), [str(hermes_python), str(Path(__file__).resolve()), *sys.argv[1:]])
    raise SystemExit(
        "ERROR: Pydantic is unavailable and the Hermes Python runtime was not found at "
        f"{hermes_python}"
    )


def resolve_skill_dir(explicit: str | None) -> Path:
    candidates = [
        Path(explicit).expanduser() if explicit else None,
        Path.home() / ".hermes/skills/patent-sequence-exraction",
        Path.home() / ".hermes/skills/patent-sequence-extraction",
        Path(__file__).resolve().parent / "patent_extraction_skill",
    ]
    for candidate in candidates:
        if candidate and (candidate / "scripts/schemas_pydantic.py").is_file():
            return candidate.resolve()
    raise SystemExit("ERROR: patent extraction skill not found; pass --skill-dir PATH")


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        os.replace(temporary, path)
    except Exception:
        if temporary and temporary.exists():
            temporary.unlink()
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate and merge completed WO*_sequence_tm.csv files."
    )
    parser.add_argument(
        "--root",
        default="/Volumes/home/output",
        help="Root directory to search recursively (default: /Volumes/home/output).",
    )
    parser.add_argument(
        "--output",
        default="/Volumes/home/output/tm_results_merged.csv",
        help="Merged CSV destination.",
    )
    parser.add_argument(
        "--job-list",
        help=(
            "Patent-folder list used for fast direct lookups. When omitted, "
            "ROOT/joblist_tm.txt is used if present; otherwise recursive discovery is used."
        ),
    )
    parser.add_argument(
        "--manifest",
        help="Successful-file list (default: OUTPUT with .files.txt suffix).",
    )
    parser.add_argument(
        "--report",
        help="JSON collection report (default: OUTPUT with .report.json suffix).",
    )
    parser.add_argument("--skill-dir", help="Patent extraction skill directory override.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit nonzero if any candidate CSV is invalid; valid output is still written.",
    )
    parser.add_argument(
        "--audit-skips",
        action="store_true",
        help=(
            "Re-open OCR JSON for skipped patents and flag possible false negatives. "
            "This is slower on network filesystems."
        ),
    )
    return parser.parse_args()


def sibling_with_suffix(path: Path, suffix: str) -> Path:
    base = path.with_suffix("") if path.suffix else path
    return base.parent / f"{base.name}{suffix}"


def candidates_from_job_list(
    job_list: Path,
    output: Path,
) -> tuple[list[Path], list[Path], list[dict[str, str]], int, int]:
    """Classify terminal artifacts without counting no-Tm skips as pending."""
    candidates: set[Path] = set()
    skip_markers: set[Path] = set()
    conflicts: list[dict[str, str]] = []
    listed_jobs = 0
    pending_jobs = 0
    with job_list.open(encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            listed_jobs += 1
            patent_folder = Path(line).expanduser()
            if not patent_folder.is_absolute():
                patent_folder = job_list.parent / patent_folder
            patent_folder = patent_folder.resolve()
            patent_id = patent_folder.name
            if not re.fullmatch(r"WO\d+", patent_id):
                pending_jobs += 1
                continue
            result = (patent_folder / "output" / f"{patent_id}_sequence_tm.csv").resolve()
            skip_marker = (
                patent_folder / "output" / f"{patent_id}_tm_skipped.json"
            ).resolve()
            has_result = result.is_file() and result != output
            has_skip = skip_marker.is_file()
            if has_result and has_skip:
                conflicts.append(
                    {
                        "patent_id": patent_id,
                        "csv_path": str(result),
                        "skip_marker_path": str(skip_marker),
                        "error": "both result CSV and no-Tm skip marker are present",
                    }
                )
            elif has_result:
                candidates.add(result)
            elif has_skip:
                skip_markers.add(skip_marker)
            else:
                pending_jobs += 1
    return sorted(candidates), sorted(skip_markers), conflicts, listed_jobs, pending_jobs


def validate_skip_marker(path: Path) -> dict[str, Any]:
    """Validate a terminal no-Tm marker without claiming it is a true negative."""
    match = re.fullmatch(r"(WO\d+)_tm_skipped\.json", path.name)
    if not match:
        raise ValueError(f"unexpected skip-marker filename: {path.name}")
    expected_patent = match.group(1)
    with path.open(encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("skip marker must contain a JSON object")
    if payload.get("patent_id_origin") != expected_patent:
        raise ValueError(
            "patent_id_origin "
            f"{payload.get('patent_id_origin')!r} does not match filename {expected_patent!r}"
        )
    if payload.get("status") != "skipped_no_tm":
        raise ValueError(f"expected status 'skipped_no_tm', got {payload.get('status')!r}")
    return {
        "patent_id": expected_patent,
        "path": str(path),
        "source": payload.get("source"),
        "timestamp": payload.get("timestamp"),
        "audit_status": "not_audited_for_false_negative",
    }


def locate_patent_json(skip_marker: Path, patent_id: str) -> Path | None:
    patent_folder = skip_marker.parent.parent
    candidates = [
        patent_folder / f"{patent_id}.json",
        patent_folder / "output" / f"{patent_id}.json",
        patent_folder / "output" / patent_id / f"{patent_id}.json",
    ]
    return next((path.resolve() for path in candidates if path.is_file()), None)


def audit_skipped_patent(item: dict[str, Any], build_cache: Any) -> dict[str, Any]:
    """Apply a contextual scan plus a deliberately broad recall-oriented check."""
    patent_id = str(item["patent_id"])
    marker_path = Path(str(item["path"]))
    patent_json = locate_patent_json(marker_path, patent_id)
    if patent_json is None:
        item["audit_status"] = "audit_error_json_not_found"
        return item
    try:
        cache = build_cache(patent_json, patent_id=patent_id)
        strict_hits = int(cache.get("tm_hit_count", 0))
        raw_text = patent_json.read_text(encoding="utf-8", errors="replace")
        broad_match = _BROAD_TM_EVIDENCE.search(raw_text)
        item["audit_json"] = str(patent_json)
        item["current_filter_hit_count"] = strict_hits
        item["broad_lexical_hit"] = bool(broad_match)
        if strict_hits:
            item["audit_status"] = "possible_false_negative_current_filter_hit"
        elif broad_match:
            item["audit_status"] = "needs_agent_review_broad_tm_evidence"
            item["broad_evidence"] = re.sub(r"\s+", " ", broad_match.group(0))[:160]
        else:
            item["audit_status"] = "no_tm_lexical_evidence"
    except Exception as exc:
        item["audit_status"] = "audit_error"
        item["audit_error"] = str(exc)
    return item


def main() -> int:
    args = parse_args()
    ensure_pydantic_runtime()

    root = Path(args.root).expanduser().resolve()
    output = Path(args.output).expanduser().resolve()
    manifest = (
        Path(args.manifest).expanduser().resolve()
        if args.manifest
        else sibling_with_suffix(output, ".files.txt")
    )
    report = (
        Path(args.report).expanduser().resolve()
        if args.report
        else sibling_with_suffix(output, ".report.json")
    )
    if not root.is_dir():
        raise SystemExit(f"ERROR: search root not found: {root}")

    explicit_job_list = Path(args.job_list).expanduser().resolve() if args.job_list else None
    if explicit_job_list and not explicit_job_list.is_file():
        raise SystemExit(f"ERROR: job list not found: {explicit_job_list}")
    automatic_job_list = root / "joblist_tm.txt"
    job_list = explicit_job_list or (automatic_job_list if automatic_job_list.is_file() else None)

    skill_dir = resolve_skill_dir(args.skill_dir)
    sys.path.insert(0, str(skill_dir / "scripts"))
    from schemas_pydantic import (  # type: ignore[import-not-found]
        SequenceTmOutputRecord,
        TM_CSV_COLUMNS,
        write_sequence_tm_csv,
    )
    if args.audit_skips:
        from tm_page_filter import build_cache as build_tm_cache  # type: ignore[import-not-found]
    else:
        build_tm_cache = None

    filename_pattern = re.compile(r"^(WO\d+)_sequence_tm\.csv$")
    if job_list:
        (
            candidates,
            skip_marker_paths,
            artifact_conflicts,
            listed_jobs,
            pending_jobs,
        ) = candidates_from_job_list(job_list, output)
        discovery_mode = "job_list_direct_lookup"
    else:
        candidates = sorted(
            path.resolve()
            for path in root.rglob("WO*_sequence_tm.csv")
            if path.is_file() and path.resolve() != output
        )
        listed_jobs = None
        pending_jobs = None
        skip_marker_paths = []
        artifact_conflicts = []
        discovery_mode = "recursive_fallback"
    successful: list[dict[str, Any]] = []
    header_only: list[dict[str, Any]] = []
    legacy_files: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    skipped_no_tm: list[dict[str, Any]] = []
    invalid_skip_markers: list[dict[str, str]] = []
    validated_records: list[Any] = []

    for path in skip_marker_paths:
        try:
            skipped_no_tm.append(validate_skip_marker(path))
        except Exception as exc:
            invalid_skip_markers.append({"path": str(path), "error": str(exc)})
    if build_tm_cache is not None:
        skipped_no_tm = [
            audit_skipped_patent(item, build_tm_cache) for item in skipped_no_tm
        ]
    skip_audit_candidates = [
        item
        for item in skipped_no_tm
        if item["audit_status"]
        in {
            "possible_false_negative_current_filter_hit",
            "needs_agent_review_broad_tm_evidence",
        }
    ]
    skip_audit_errors = [
        item for item in skipped_no_tm if str(item["audit_status"]).startswith("audit_error")
    ]

    for path in candidates:
        match = filename_pattern.fullmatch(path.name)
        if not match:
            continue
        expected_patent = match.group(1)
        file_records: list[Any] = []
        try:
            with path.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames == LEGACY_TM_CSV_COLUMNS:
                    legacy_files.append(
                        {
                            "patent_id": expected_patent,
                            "path": str(path),
                            "header": reader.fieldnames,
                            "status": "legacy_five_column_output_requires_force_rerun",
                        }
                    )
                    continue
                if reader.fieldnames != TM_CSV_COLUMNS:
                    raise ValueError(
                        f"expected header {TM_CSV_COLUMNS}, got {reader.fieldnames}"
                    )
                for line_number, row in enumerate(reader, start=2):
                    try:
                        record = SequenceTmOutputRecord.model_validate(dict(row))
                    except Exception as exc:
                        raise ValueError(f"line {line_number}: {exc}") from exc
                    if record.patent_id_origin != expected_patent:
                        raise ValueError(
                            f"line {line_number}: patent_id_origin "
                            f"{record.patent_id_origin!r} does not match filename {expected_patent!r}"
                        )
                    file_records.append(record)
            file_result = {
                "patent_id": expected_patent,
                "path": str(path),
                "rows": len(file_records),
            }
            if file_records:
                validated_records.extend(file_records)
                successful.append(file_result)
            else:
                header_only.append(file_result)
        except Exception as exc:
            failures.append({"path": str(path), "error": str(exc)})

    # Exact row deduplication keeps reruns or duplicate patent trees merge-safe.
    unique: dict[tuple[str, ...], Any] = {}
    for record in validated_records:
        row = record.to_csv_row()
        unique[tuple(row[column] for column in TM_CSV_COLUMNS)] = record
    ordered = sorted(
        unique.values(),
        key=lambda record: (
            record.patent_id_origin,
            record.sense_sequence,
            record.antisense_sequence,
            record.tm_value,
            record.modification_sense,
            record.modification_antisense,
        ),
    )

    # The skill writer revalidates the complete batch and atomically publishes it.
    write_sequence_tm_csv(ordered, output)
    atomic_text(
        manifest,
        "".join(f"{item['path']}\n" for item in successful),
    )
    report_payload = {
        "root": str(root),
        "output": str(output),
        "discovery_mode": discovery_mode,
        "job_list": str(job_list) if job_list else None,
        "listed_jobs": listed_jobs,
        "results_not_yet_present": pending_jobs,
        "true_pending_jobs": pending_jobs,
        "candidate_files": len(candidates),
        "valid_files": len(successful) + len(header_only),
        "successful_files": len(successful),
        "header_only_files": len(header_only),
        "legacy_files": len(legacy_files),
        "skipped_no_tm": len(skipped_no_tm),
        "skip_audit_enabled": args.audit_skips,
        "skip_audit_candidates": len(skip_audit_candidates),
        "skip_audit_errors": len(skip_audit_errors),
        "artifact_conflicts": len(artifact_conflicts),
        "invalid_skip_markers": len(invalid_skip_markers),
        "invalid_files": len(failures),
        "input_rows": len(validated_records),
        "merged_unique_rows": len(ordered),
        "successful": successful,
        "header_only": header_only,
        "legacy": legacy_files,
        "skipped": skipped_no_tm,
        "skip_false_negative_policy": (
            "Skip markers are terminal for pending-count purposes but are not treated "
            "as verified negatives. Use --audit-skips for a slower contextual and broad "
            "lexical recall check; broad-only hits still require agent review."
        ),
        "skip_audit_candidate_details": skip_audit_candidates,
        "skip_audit_error_details": skip_audit_errors,
        "conflicts": artifact_conflicts,
        "invalid_skip_marker_details": invalid_skip_markers,
        "failures": failures,
    }
    atomic_text(report, json.dumps(report_payload, indent=2, ensure_ascii=False) + "\n")

    print(f"Successful files : {len(successful)}")
    print(f"Header-only files: {len(header_only)}")
    print(f"Legacy files     : {len(legacy_files)}")
    skip_label = "audited" if args.audit_skips else "not false-negative audited"
    print(f"Skipped no Tm    : {len(skipped_no_tm)} ({skip_label})")
    if args.audit_skips:
        print(f"Skip audit flags : {len(skip_audit_candidates)}")
        print(f"Skip audit errors: {len(skip_audit_errors)}")
    print(f"Artifact conflicts: {len(artifact_conflicts)}")
    print(f"Invalid skip marks: {len(invalid_skip_markers)}")
    print(f"Invalid files    : {len(failures)}")
    print(f"Merged rows      : {len(ordered)}")
    print(f"Discovery        : {discovery_mode}")
    if job_list:
        print(f"Job list         : {job_list}")
        print(f"True pending     : {pending_jobs}")
    print(f"Merged CSV       : {output}")
    print(f"File manifest    : {manifest}")
    print(f"Report           : {report}")
    strict_errors = failures or invalid_skip_markers or artifact_conflicts
    return 1 if args.strict and strict_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

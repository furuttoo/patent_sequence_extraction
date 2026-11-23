"""Filter OCR patent pages for genuine melting-temperature evidence.

An isolated ``Tm`` nucleotide-modification token is not a hit. A page must
contain explicit melting/thermal language or a Tm table/header with sequence
or Celsius context.
"""

from __future__ import annotations

import html
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable


_DIRECT_TM = re.compile(
    r"melting\s+(?:temperature|point)|thermal\s+(?:melting|stability|denaturation)",
    re.I,
)
_TM_TOKEN = re.compile(
    r"(?<![A-Za-z0-9])T\s*\$?\s*(?:[_-]\s*\{*\s*)?m\s*\}*\s*\$?(?![A-Za-z0-9])",
    re.I,
)
_CELSIUS = re.compile(r"(?:°\s*C|degrees?\s*C|Celsius)", re.I)
_SEQUENCE_CONTEXT = re.compile(
    r"sequence|oligonucleotide|oligo|duplex|sense|antisense|SEQ\s+ID",
    re.I,
)
_TM_VALUE_CONTEXT = re.compile(r"\bTm\s+(?:value|measurement|was|were|of)\b", re.I)
_TM_HEADER_CELL = re.compile(r"^\s*T\s*[_-]?m(?:\s*\([^)]*C[^)]*\)|\s*°?\s*C)?\s*$", re.I)


@dataclass
class TmPageResult:
    page_index: int
    is_tm_hit: bool
    has_table: bool
    evidence: list[str] = field(default_factory=list)


def _decode(value: Any) -> str:
    return html.unescape(str(value or ""))


def _has_tm_header_cell(table_html: str) -> bool:
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", table_html, flags=re.I | re.S)[:5]:
        cells = re.findall(r"<t[hd][^>]*>(.*?)</t[hd]>", row, flags=re.I | re.S)
        for cell in cells:
            plain = re.sub(r"<[^>]+>", "", cell)
            decoded = html.unescape(plain)
            compact = re.sub(r"[^A-Za-z0-9]+", "", decoded).casefold()
            if (
                _TM_HEADER_CELL.fullmatch(decoded)
                or compact in {"tm", "tmc", "tmoc", "tmdegc"}
                or "meltingtemperature" in compact
            ):
                return True
    return False


def _schema_a_blocks(page: dict[str, Any]) -> list[dict[str, str]]:
    markdown = page.get("markdown", "")
    if isinstance(markdown, str):
        try:
            decoded = json.loads(markdown)
        except (json.JSONDecodeError, TypeError):
            decoded = {"text": markdown}
    elif isinstance(markdown, dict):
        decoded = markdown
    else:
        decoded = {}
    return [{"block_label": "text", "block_content": _decode(decoded.get("text", ""))}]


def iter_pages(data: Any) -> Iterable[dict[str, Any]]:
    """Yield normalized ``{page_index, blocks}`` pages for schemas A/B/C."""
    if isinstance(data, list):
        for index, page in enumerate(data):
            yield {
                "page_index": page.get("page_index", index),
                "blocks": page.get("blocks") or _schema_a_blocks(page),
            }
        return
    if not isinstance(data, dict):
        raise ValueError("patent JSON must be a list or object")
    if isinstance(data.get("pages"), list):
        for index, page in enumerate(data["pages"]):
            yield {
                "page_index": page.get("page_index", index),
                "blocks": page.get("blocks") or _schema_a_blocks(page),
            }
        return
    if isinstance(data.get("blocks"), list):
        yield {"page_index": data.get("page_index", 0), "blocks": data["blocks"]}
        return

    page_values = [value for value in data.values() if isinstance(value, dict) and "blocks" in value]
    if page_values:
        for index, page in enumerate(page_values):
            yield {"page_index": page.get("page_index", index), "blocks": page["blocks"]}
        return
    raise ValueError("unsupported patent JSON schema: no pages or blocks found")


def score_page(page: dict[str, Any]) -> TmPageResult:
    blocks = page.get("blocks", [])
    decoded_blocks = [
        (str(block.get("block_label", "text")), _decode(block.get("block_content", "")))
        for block in blocks
    ]
    has_table = any(label == "table" for label, _ in decoded_blocks)
    evidence: list[str] = []

    for label, content in decoded_blocks:
        direct = _DIRECT_TM.search(content)
        if direct:
            evidence.append(direct.group(0))
            continue
        token = _TM_TOKEN.search(content)
        if not token:
            continue
        local_context = content[max(0, token.start() - 40) : token.end() + 40]
        contextual = (
            _CELSIUS.search(local_context)
            or _TM_VALUE_CONTEXT.search(content)
            or (
                label == "table"
                and _has_tm_header_cell(content)
                and _SEQUENCE_CONTEXT.search(content)
            )
        )
        if contextual:
            start = max(0, token.start() - 40)
            end = min(len(content), token.end() + 80)
            evidence.append(re.sub(r"\s+", " ", content[start:end]).strip())

    return TmPageResult(
        page_index=int(page.get("page_index", 0)),
        is_tm_hit=bool(evidence),
        has_table=has_table,
        evidence=evidence,
    )


def scan(patent_json: str | Path) -> tuple[Any, list[dict[str, Any]], list[TmPageResult]]:
    with open(patent_json, encoding="utf-8") as handle:
        data = json.load(handle)
    pages = list(iter_pages(data))
    results = [score_page(page) for page in pages]
    return data, pages, results


def build_cache(
    patent_json: str | Path,
    *,
    patent_id: str = "",
    include_table_continuations: bool = True,
) -> dict[str, Any]:
    data, pages, results = scan(patent_json)
    by_index = {int(page["page_index"]): page for page in pages}
    result_by_index = {result.page_index: result for result in results}
    selected = {result.page_index for result in results if result.is_tm_hit}

    if include_table_continuations:
        for page_index in list(selected):
            result = result_by_index[page_index]
            next_index = page_index + 1
            next_page = by_index.get(next_index)
            if result.has_table and next_page:
                next_result = result_by_index[next_index]
                next_text = " ".join(
                    _decode(block.get("block_content", "")) for block in next_page.get("blocks", [])
                )
                if next_result.has_table and (
                    re.search(r"\bcontinued\b", next_text, re.I)
                    or _SEQUENCE_CONTEXT.search(next_text)
                ):
                    selected.add(next_index)

    patent_from_data = data.get("patent_id", "") if isinstance(data, dict) else ""
    pages_out = []
    for page_index in sorted(selected):
        page = by_index[page_index]
        result = result_by_index[page_index]
        pages_out.append(
            {
                "page_index": page_index,
                "direct_tm_hit": result.is_tm_hit,
                "evidence": result.evidence,
                "blocks": [
                    {
                        "block_label": block.get("block_label", "text"),
                        "block_content": _decode(block.get("block_content", "")),
                    }
                    for block in page.get("blocks", [])
                ],
            }
        )

    return {
        "patent_id": patent_id or patent_from_data or Path(patent_json).stem,
        "total_pages": len(pages),
        "tm_hit_count": sum(result.is_tm_hit for result in results),
        "skipped_no_tm": not any(result.is_tm_hit for result in results),
        "abbrev_map": {},
        "pages": pages_out,
    }

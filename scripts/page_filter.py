"""
page_filter.py  —  Pass 1: regex pool builder
==============================================
Scans an OCR'd siRNA patent JSON and returns the subset of pages worth
caching for downstream extraction.  No scoring, no heuristics: a page is
either in the pool or it isn't.

Five binary signals, applied in order:

  ABBREV   — abbreviation / legend table (style E)
  SEQ      — ≥2 nucleotide runs of ≥15 chars
  MOD_A    — inline suffix modifier (m/f/n) inside a nucleotide run
  MOD_B    — mixed-case nucleotide run (lowercase = modified convention)
  ACTIVITY — activity keyword AND (table block OR duplex name on same page)

MOD_A and MOD_B never pull in pages that SEQ hasn't already caught, but they
annotate *how* sequences on a page are encoded so the extraction agent knows
which decoding path to take without re-reading the page.

ABBREV pages are handled specially: the agent (not this module) builds the
actual letter-code → chemistry mapping by reading those pages.  This module
extracts the raw table HTML and any notation footnotes verbatim so the agent
has everything it needs in one place.

Pass 2 (abbreviation map building) and Pass 3 (drag / forward-search for
missing notation keys) are separate concerns handled by the agent.
"""

from __future__ import annotations

import html
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Compiled patterns
# ---------------------------------------------------------------------------

# Nucleotide run: base chars + modifier chars, minimum 15 long.
# Character class is broad; length is the noise gate. 'o' is the alkoxy
# modifier: WO2026131490's primary 表1 (p111) is written GosGosCoAo..., and
# without 'o' the page scored zero runs, so the one table holding the patent's
# siRNA sequences was never pooled and the agent re-scanned the whole source.
_NUC_RUN = re.compile(r'[ACGUTacgut][ACGUTacgutmfnsdpo]{14,}')

# Style A: modifier letter sandwiched between nucleotide-class chars.
# Applied only inside confirmed nucleotide runs, never against full page text.
_MOD_A = re.compile(r'[ACGUTacgut][mfn](?=[ACGUTacgut])')

# Abbreviation legend table: both header forms seen across Alnylam filings.
_ABBREV = re.compile(r'Abbreviation\s*\$?\s*\^|Nucleotide\(s\)', re.I)

# Activity: assay-type keyword. The English terms alone missed every Chinese
# filing: WO2026131490's "抑制率" table (p152) and its in-vivo pages 153-155,
# and WO2026098677's "mRNA 剩余表达水平%" table (p57-73) scored zero, so the
# cache routed activity to nothing (or to the search report) and the agent went
# hunting through the full source and the 547 KB EPO context instead. Those
# pages are the measurements; without them Step 5 has no input.
_ACTIVITY = re.compile(
    r'inhibit|knockdown|silenc|IC50|EC50|ED50|mRNA\s+reduction'
    r'|%\s*inhib|Luc\s+reporter'
    r'|抑制率|抑制活性|敲低率|敲低效率|沉默率|沉默效率'
    r'|剩余表达|残留表达|表达量|表达水平|半数抑制|半抑制浓度',
    re.I,
)

# The international search report quotes the claims and repeats "inhibit" /
# "抑制" in prose and tables. It is not a measurement page: on WO2026098677 it
# was the ONLY page flagged ACTIVITY (p87/93), which is exactly the wrong place
# to send Step 5.
_SEARCH_REPORT = re.compile(r'INTERNATIONAL\s+SEARCH\s+REPORT|国际检索报告', re.I)

# Gate: activity only kept when a table block or duplex name is also present.
# Cuts ~29 prose-only pages that mention "inhibition" in passing.
_DUPLEX = re.compile(r'(?:AL[-\u2013]DP[-\u2013]|AD[-\u2013]|AL\d{3})\d+', re.I)

# A data cell: one number, optionally signed/approximate, optionally a percent.
# Used only to recognise the data-row pages of a table that starts earlier.
# Every signal above keys on a word -- a metric name, a nucleotide run, a legend
# header -- and the second and third pages of a long table carry none of them:
# they open mid-table with rows like <td>BC-100042</td><td>88.6</td> under a
# running header. So they scored zero and were dropped, taking the measurements
# with them. On WO2024173593 that was 11 pages holding 1,314 numeric cells,
# including the 88.6 / 90.5 / 86.4 / 66.3 percent-knockdown columns of pages
# 74-78 and the 100 nM dose-response headers of pages 84-86. Three separate
# extraction runs reported it before it was believed: "Phase 1's filter skipped
# continuation pages (74-76, 78, 84-87, 90, 93, 126) that hold table rows".
_NUMERIC_CELL = re.compile(r'<td[^>]*>\s*[<>~\u2248\u2264\u2265]?\s*\d+(?:[.,]\d+)?\s*%?\s*</td>')
# Three rows of numbers: enough to be a data table, few enough to catch a short
# continuation. A page with one or two stray numeric cells is prose.
_MIN_DATA_CELLS = 3


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PageResult:
    page_index: int
    in_pool:    bool

    # which signals fired
    sig_abbrev:   bool = False   # style E: abbreviation legend table
    sig_seq:      bool = False   # ≥2 nucleotide runs ≥15 nt
    sig_mod_a:    bool = False   # style A: inline suffix modifier in run
    sig_mod_b:    bool = False   # style B: mixed-case run
    sig_activity: bool = False   # gated activity signal
    sig_table_cont: bool = False # data-row page of a table that began earlier

    has_table: bool = False
    numeric_cells: int = 0       # data cells, for continuation detection

    # for SEQ pages: representative runs so agent can see encoding at a glance
    sequence_samples: list[str] = field(default_factory=list)

    # for ABBREV pages: raw table HTML(s) + any notation footnotes on the page,
    # verbatim and decoded — the agent reads these to build the letter→chemistry map
    abbrev_tables:    list[str] = field(default_factory=list)   # one entry per <table>
    abbrev_footnotes: list[str] = field(default_factory=list)   # vision_footnote / text blocks
                                                                 # containing notation keys

    def signals(self) -> list[str]:
        out = []
        if self.sig_abbrev:   out.append('ABBREV')
        if self.sig_seq:      out.append('SEQ')
        if self.sig_mod_a:    out.append('MOD-A')
        if self.sig_mod_b:    out.append('MOD-B')
        if self.sig_activity: out.append('ACTIVITY')
        # Last, and additive: a continuation inherits the signals of the page
        # whose table it continues, so every consumer that routes on SEQ or
        # ACTIVITY keeps working, and TABLE-CONT only tells the agent that this
        # page opens mid-table and its columns are defined on an earlier page.
        if self.sig_table_cont: out.append('TABLE-CONT')
        return out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Notation-key phrases that appear in text/footnote blocks on or near
# abbreviation pages.  Kept deliberately broad — the agent interprets content.
_NOTATION_KEY = re.compile(
    r'lower\s*case|upper\s*case|capital\s+letter|denoted\s+by|'
    r'represent\s+(?:ribo|deoxy)|standard\s+nomenclature|abbreviations?\s+of',
    re.I,
)


def _decode(raw: str) -> str:
    return html.unescape(raw)


def _page_blocks(page: dict[str, Any]) -> list[tuple[str, str]]:
    """Return list of (block_label, decoded_content) for a page."""
    return [
        (b['block_label'], _decode(b.get('block_content', '')))
        for b in page.get('blocks', [])
    ]


def _nuc_runs(text: str) -> list[str]:
    return [m.group() for m in _NUC_RUN.finditer(text) if len(m.group()) >= 15]


def _mod_b(run: str) -> bool:
    """True if run is mixed-case after stripping the common dT/TT overhang."""
    core = re.sub(r'(?:dT|TT)$', '', run)
    return (
        len(core) >= 10
        and bool(re.search(r'[acgut]', core))
        and bool(re.search(r'[ACGU]', core))
    )


def _extract_abbrev(blocks: list[tuple[str, str]]) -> tuple[list[str], list[str]]:
    """
    Pull raw table HTML and notation footnotes from an ABBREV page.

    Returns:
        tables    — list of raw <table>...</table> strings (HTML-decoded)
        footnotes — list of text/footnote block contents that contain
                    notation-key phrases (e.g. "lower case = ribonucleotide")
    """
    tables:    list[str] = []
    footnotes: list[str] = []

    for label, content in blocks:
        if label == 'table':
            tables.append(content)
        elif label in ('vision_footnote', 'text', 'figure_title', 'paragraph_title'):
            if _NOTATION_KEY.search(content):
                footnotes.append(content)

    return tables, footnotes


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def score_page(page: dict[str, Any]) -> PageResult:
    """Run all signals on one page dict and return a PageResult."""
    pi     = page['page_index']
    blocks = _page_blocks(page)
    text   = ' '.join(c for _, c in blocks)

    has_table = any(label == 'table' for label, _ in blocks)
    runs      = _nuc_runs(text)

    sig_abbrev   = bool(_ABBREV.search(text))
    sig_seq      = len(runs) >= 2
    sig_mod_a    = any(_MOD_A.search(r) for r in runs)
    sig_mod_b    = any(_mod_b(r) for r in runs)
    sig_activity = (
        bool(_ACTIVITY.search(text))
        and (has_table or bool(_DUPLEX.search(text)))
        and not _SEARCH_REPORT.search(text)
    )

    in_pool = any([sig_abbrev, sig_seq, sig_mod_a, sig_mod_b, sig_activity])

    result = PageResult(
        page_index   = pi,
        in_pool      = in_pool,
        sig_abbrev   = sig_abbrev,
        sig_seq      = sig_seq,
        sig_mod_a    = sig_mod_a,
        sig_mod_b    = sig_mod_b,
        sig_activity = sig_activity,
        has_table    = has_table,
        numeric_cells = len(_NUMERIC_CELL.findall(text)),
        sequence_samples = runs[:4],
    )

    if sig_abbrev:
        result.abbrev_tables, result.abbrev_footnotes = _extract_abbrev(blocks)

    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def link_table_continuations(
    results:   list[PageResult],
    min_cells: int = _MIN_DATA_CELLS,
) -> list[PageResult]:
    """Pool the data-row pages that follow a pooled page's table.

    Walks forward in page_index order, so a table spanning five pages is
    carried the whole way: each page admitted becomes the anchor for the next.
    A page joins only when all three hold -- it is not already pooled, the page
    immediately before it is, and it carries at least ``min_cells`` numeric
    table cells of its own. The adjacency requirement is what keeps this from
    sweeping in every table in the patent: a data page with no pooled page in
    front of it is not a continuation of anything.

    The admitted page inherits its anchor's signals verbatim rather than being
    labelled activity-by-assumption. That matters downstream: round_kind() in
    run_activity_extraction treats a round as activity-only when no page in it
    carries SEQ/MOD-A/MOD-B/ABBREV, so a sequence table's continuation has to
    keep saying SEQ or its round would be told not to create records.

    A continuation of an ABBREV legend keeps the signal but not the extracted
    abbrev_tables; its rows are in ``blocks``, which is what the agent reads.
    """
    for index, result in enumerate(results):
        if index == 0 or result.in_pool:
            continue
        anchor = results[index - 1]
        if not anchor.in_pool or anchor.page_index != result.page_index - 1:
            continue
        if not result.has_table or result.numeric_cells < min_cells:
            continue
        result.sig_abbrev = anchor.sig_abbrev
        result.sig_seq = anchor.sig_seq
        result.sig_mod_a = anchor.sig_mod_a
        result.sig_mod_b = anchor.sig_mod_b
        result.sig_activity = anchor.sig_activity
        result.sig_table_cont = True
        result.in_pool = True
    return results


def scan(patent_json: str | Path) -> list[PageResult]:
    """
    Scan every page in a patent JSON.
    Returns one PageResult per page, in page_index order.
    Check .in_pool on each to identify cache members.
    """
    with open(patent_json) as f:
        data = json.load(f)
    return link_table_continuations([score_page(p) for p in data['pages']])


def pool_pages(results: list[PageResult]) -> list[PageResult]:
    """Return only in_pool pages, preserving page_index order."""
    return [r for r in results if r.in_pool]


def build_cache(
    results:     list[PageResult],
    patent_json: str | Path,
    patent_id:   str = '',
) -> dict[str, Any]:
    """
    Serialise the pool to a cache dict (write as <patent>_page_cache.json).

    Structure
    ---------
    {
      patent_id, total_pages, pool_size, pages_skipped, detected_signals,

      abbrev_map_pending: true,   # agent must populate this before extraction
                                  # by reading the abbrev pages listed below

      pages: [
        {
          page_index, signals, has_table, sequence_samples,

          # ABBREV pages only:
          abbrev_tables:    [...],   # raw <table> HTML for agent to parse
          abbrev_footnotes: [...],   # notation-key sentences for agent to read

          blocks: [{ block_label, block_content }, ...]
        },
        ...
      ]
    }

    All block_content values are HTML-decoded.  The agent never needs to
    re-read the source patent JSON.
    """
    with open(patent_json) as f:
        raw_pages = {p['page_index']: p for p in json.load(f)['pages']}

    pool = pool_pages(results)

    detected: set[str] = set()
    for r in pool:
        detected.update(r.signals())

    pages_out = []
    for r in pool:
        entry: dict[str, Any] = {
            'page_index':       r.page_index,
            'signals':          r.signals(),
            'has_table':        r.has_table,
            'sequence_samples': r.sequence_samples,
            'blocks': [
                {
                    'block_label':   b['block_label'],
                    'block_content': _decode(b.get('block_content', '')),
                }
                for b in raw_pages[r.page_index].get('blocks', [])
            ],
        }

        if r.sig_abbrev:
            entry['abbrev_tables']    = r.abbrev_tables
            entry['abbrev_footnotes'] = r.abbrev_footnotes

        pages_out.append(entry)

    return {
        'patent_id':        patent_id or Path(patent_json).stem,
        'total_pages':      len(results),
        'pool_size':        len(pool),
        'pages_skipped':    len(results) - len(pool),
        'detected_signals': sorted(detected),
        'abbrev_map':       {},      # agent writes letter→chemistry pairs here
        'pages':            pages_out,
    }

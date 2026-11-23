"""
schemas.py — siRNA record shape, validation, and CSV writer

Three classes define the data contract:

    Modification  — one chemical modification at one position on one strand
    Activity      — one activity measurement (in-vitro or in-vivo)
    SiRNARecord   — one strand record (sense OR antisense); the agent builds
                    two per duplex and passes them to write_activity_csv()

The agent should:
    1. Import SiRNARecord, Modification, Activity from this file.
    2. Build records by constructing SiRNARecord(**data) — pydantic validates
       field names and types at construction; extra="forbid" rejects any field
       not declared below.
    3. Call write_activity_csv(records, path) to flatten and write the CSV.
       Never write the CSV directly — to_csv_row() is the only flattening path.

Modification decoding rules: modification_decoding_reference.md
Patent-specific symbol→chemistry map: abbrev_map in the page cache
"""

from __future__ import annotations

import csv
import json
import os
import re
import tempfile
import unicodedata
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path
from typing import Any, Literal, Mapping, Optional
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from pydantic import BaseModel, Field, PrivateAttr, ValidationInfo, field_serializer, field_validator, model_validator


# ---------------------------------------------------------------------------
# Normalisation helpers
# ---------------------------------------------------------------------------

_WHITESPACE_RE = re.compile(r"\s+")
_NUMBER_PATTERN = r"(?:\d+(?:[.,]\d+)?|\.\d+)(?:[eE][+-]?\d+)?"
_UNIT_PATTERN = r"(?:pmol/L|nmol/L|umol/L|mmol/L|mol/L|pM|nM|uM|mM|M|%)"

_CONCENTRATION_FACTORS_TO_NM = {
    "M": Decimal("1000000000"),
    "mM": Decimal("1000000"),
    "uM": Decimal("1000"),
    "nM": Decimal("1"),
    "pM": Decimal("0.001"),
}

_CONCENTRATION_UNIT_ALIASES = {
    "m": "M",
    "mol/l": "M",
    "mm": "mM",
    "mmol/l": "mM",
    "um": "uM",
    "umol/l": "uM",
    "nm": "nM",
    "nmol/l": "nM",
    "pm": "pM",
    "pmol/l": "pM",
}

_PERCENT_METRICS = {
    "knockdown_percent",
    "mRNA_remaining_percent",
    "cell_viability_percent",
}
_CONCENTRATION_METRICS = {"IC50", "EC50"}


def _clean_text(value: Any) -> Optional[str]:
    """Coerce scalar input to a clean, Unicode-normalised string."""
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError("boolean values are not valid scientific measurements")
    text = unicodedata.normalize("NFKC", str(value))
    text = text.replace("μ", "u").replace("µ", "u")
    text = _WHITESPACE_RE.sub(" ", text).strip()
    return text or None


def _decimal(value: str) -> Decimal:
    value = value.replace(",", ".")
    try:
        result = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f"invalid numeric value: {value!r}") from exc
    if not result.is_finite():
        raise ValueError(f"non-finite numeric value: {value!r}")
    return result


def _format_decimal(value: Decimal) -> str:
    if value == 0:
        return "0"
    rendered = format(value.normalize(), "f")
    return rendered.rstrip("0").rstrip(".") if "." in rendered else rendered


def _canonical_concentration_unit(unit: Optional[str]) -> Optional[str]:
    if unit is None:
        return None
    cleaned = unit.replace("μ", "u").replace("µ", "u").replace(" ", "")
    canonical = _CONCENTRATION_UNIT_ALIASES.get(cleaned.casefold())
    if canonical is None:
        raise ValueError(f"unsupported concentration unit: {unit!r}")
    return canonical


@dataclass(frozen=True)
class _NormalisedQuantity:
    value: str
    warning: Optional[str] = None


def _normalise_quantity(
    value: Any,
    *,
    quantity_type: Literal["concentration", "percent"],
    default_unit: Optional[str] = None,
) -> _NormalisedQuantity:
    """Normalise a point, interval, inequality, or value with uncertainty."""
    text = _clean_text(value)
    if text is None:
        return _NormalisedQuantity("")

    text = re.sub(r"(?i)\bpercent\b", "%", text)
    text = re.sub(r"(?i)^\s*(?:approximately|approx\.?|about)\s*", "~", text)
    text = text.replace("≤", "<=").replace("≥", ">=").replace("−", "-")
    text = text.replace("—", "–")

    relation = ""
    relation_match = re.match(r"^(<=|>=|<|>|~)\s*", text)
    if relation_match:
        relation = relation_match.group(1)
        text = text[relation_match.end():].strip()

    range_match = re.fullmatch(
        rf"(?P<low>{_NUMBER_PATTERN})\s*(?P<low_unit>{_UNIT_PATTERN})?\s*"
        rf"(?:-|–|(?i:to))\s*"
        rf"(?P<high>{_NUMBER_PATTERN})\s*(?P<high_unit>{_UNIT_PATTERN})?",
        text,
        flags=re.IGNORECASE,
    )
    uncertainty_match = re.fullmatch(
        rf"(?P<point>{_NUMBER_PATTERN})\s*(?:±|\+/-)\s*"
        rf"(?P<uncertainty>{_NUMBER_PATTERN})\s*(?P<unit>{_UNIT_PATTERN})?",
        text,
        flags=re.IGNORECASE,
    )
    point_match = re.fullmatch(
        rf"(?P<point>{_NUMBER_PATTERN})\s*(?P<unit>{_UNIT_PATTERN})?",
        text,
        flags=re.IGNORECASE,
    )

    expected_unit = "%" if quantity_type == "percent" else None

    def resolve_unit(raw_unit: Optional[str]) -> Optional[str]:
        candidate = raw_unit or default_unit or expected_unit
        if quantity_type == "percent":
            if candidate not in (None, "%"):
                raise ValueError(f"percentage value cannot use unit {candidate!r}")
            return "%"
        if candidate == "%":
            raise ValueError("concentration cannot be expressed as a percentage")
        return _canonical_concentration_unit(candidate)

    def convert(number: Decimal, unit: Optional[str]) -> Decimal:
        if quantity_type == "concentration" and unit is not None:
            return number * _CONCENTRATION_FACTORS_TO_NM[unit]
        return number

    def validate_bounds(*numbers: Decimal) -> None:
        if any(number < 0 for number in numbers):
            raise ValueError("measurement cannot be negative")
        if quantity_type == "percent" and any(number > 100 for number in numbers):
            raise ValueError("percentage values must be between 0 and 100")

    if range_match:
        if relation:
            raise ValueError("a range cannot also have an inequality relation")
        low_unit = resolve_unit(range_match.group("low_unit") or range_match.group("high_unit"))
        high_unit = resolve_unit(range_match.group("high_unit") or range_match.group("low_unit"))
        low = convert(_decimal(range_match.group("low")), low_unit)
        high = convert(_decimal(range_match.group("high")), high_unit)
        validate_bounds(low, high)
        if low > high:
            raise ValueError("range lower bound exceeds upper bound")
        unit = "%" if quantity_type == "percent" else ("nM" if low_unit else None)
        suffix = unit or ""
        warning = None if unit else "concentration_unit_missing"
        return _NormalisedQuantity(
            f"{_format_decimal(low)}–{_format_decimal(high)}{suffix if suffix == '%' else (' ' + suffix if suffix else '')}",
            warning,
        )

    if uncertainty_match:
        unit = resolve_unit(uncertainty_match.group("unit"))
        point = convert(_decimal(uncertainty_match.group("point")), unit)
        uncertainty = convert(_decimal(uncertainty_match.group("uncertainty")), unit)
        validate_bounds(point, uncertainty)
        output_unit = "%" if quantity_type == "percent" else ("nM" if unit else None)
        suffix = output_unit or ""
        warning = None if output_unit else "concentration_unit_missing"
        return _NormalisedQuantity(
            f"{relation}{_format_decimal(point)} ± {_format_decimal(uncertainty)}"
            f"{suffix if suffix == '%' else (' ' + suffix if suffix else '')}",
            warning,
        )

    if point_match:
        unit = resolve_unit(point_match.group("unit"))
        point = convert(_decimal(point_match.group("point")), unit)
        validate_bounds(point)
        output_unit = "%" if quantity_type == "percent" else ("nM" if unit else None)
        suffix = output_unit or ""
        warning = None if output_unit else "concentration_unit_missing"
        return _NormalisedQuantity(
            f"{relation}{_format_decimal(point)}"
            f"{suffix if suffix == '%' else (' ' + suffix if suffix else '')}",
            warning,
        )

    return _NormalisedQuantity(_clean_text(value) or "", "unparsed_quantity")


def _normalise_percentage(value: Any, *, remaining: bool = False) -> _NormalisedQuantity:
    """Accept an explicit percentage, preserving the source decimal precision."""
    text = _clean_text(value)
    if text is None:
        return _NormalisedQuantity("")

    match = re.fullmatch(r"(?P<value>[+-]?\d+(?:[.,]\d+)?)\s*%", text)
    if not match:
        raise ValueError(
            "percentage must be a number followed by '%'; "
            "ranges, inequalities, approximations, and uncertainty are not allowed"
        )

    number = _decimal(match.group("value"))
    if (remaining and number < 0) or (not remaining and number > Decimal("100")):
        raise ValueError("remaining percentage must be nonnegative; knockdown cannot exceed 100%")
    return _NormalisedQuantity(f"{_format_decimal(number)}%")


# What `metric="other"` may hold. "other" exists for a genuine readout that none
# of the named metrics fits -- a fold-change, a log reduction, a ratio to
# control. It is not a free-text field, and treating it as one had a measurable
# cost: of 2,047 rows exported across the first 12 patents, 136 carried
# something that is not an activity magnitude at all -- 20 rows of "18 weeks"
# (a culture duration), 6 of "t1/2 = 32 days" (serum stability), 9 of
# "effective at inhibiting angiogenesis" (a conclusion), and 91 with a metric
# and no value. All were dropped downstream, but only after being written,
# counted, and then read back as activity data the merge filter had "lost".
_TIME_UNIT = (r'(?:s|sec|secs|second|seconds|min|mins|minute|minutes|h|hr|hrs|'
              r'hour|hours|d|day|days|wk|wks|week|weeks|mo|month|months|yr|yrs|'
              r'year|years)')
# A half-life is a stability readout -- how long the oligo survives, not how
# much target it silenced. Named explicitly because it carries a number and a
# unit and would pass every other check here.
_HALF_LIFE_RE = re.compile(r't\s*(?:1\s*/\s*2|½)|half[-\s]?li(?:fe|ves)', re.I)
# A duration standing alone: one number, a time unit, and nothing else numeric.
# "72h post-dose, 40% knockdown" keeps its second number and is not matched.
_DURATION_ONLY_RE = re.compile(rf'\D*\d+(?:[.,]\d+)?\s*{_TIME_UNIT}\b\D*', re.I)
_DIGIT_RE = re.compile(r'\d')


def _non_measurement_kind(value: Optional[str]) -> Optional[str]:
    """Why an `other` value is not an activity magnitude; None if it is one.

    Order matters: a half-life is checked before the duration rule because
    "t1/2 = 15 seconds (control)" satisfies both, and "half_life" is the more
    useful label for anyone who later wants those rows back.
    """
    if value is None:
        return "absent"
    if _HALF_LIFE_RE.search(value):
        return "half_life"
    if not _DIGIT_RE.search(value):
        return "qualitative"
    if _DURATION_ONLY_RE.fullmatch(value):
        return "duration"
    return None


@dataclass(frozen=True)
class GeneResolution:
    symbol: str
    warning: Optional[str] = None
    source: Optional[str] = None


class GeneSymbolResolver:
    """
    Resolve gene aliases to canonical symbols without trusting fuzzy ranking.

    HGNC is authoritative for human genes. MyGene.info is an exact-field
    fallback and supports other species. Patent-specific ambiguous terms can
    be supplied through overrides, e.g. {"KSP": "KIF11"} when the patent
    explicitly defines KSP as kinesin spindle protein.
    """

    def __init__(
        self,
        *,
        overrides: Optional[Mapping[str, str]] = None,
        allow_network: bool = True,
        timeout_seconds: float = 5.0,
    ) -> None:
        self.overrides = {
            _clean_text(key).casefold(): _clean_text(symbol)
            for key, symbol in (overrides or {}).items()
            if _clean_text(key) and _clean_text(symbol)
        }
        self.allow_network = allow_network
        self.timeout_seconds = timeout_seconds
        self._cache: dict[tuple[str, int], GeneResolution] = {}

    def resolve(self, value: str, *, taxid: int = 9606) -> GeneResolution:
        cleaned = _clean_text(value)
        if not cleaned:
            return GeneResolution("")
        key = (cleaned.casefold(), taxid)
        if key in self._cache:
            return self._cache[key]

        override = self.overrides.get(cleaned.casefold())
        if override:
            result = GeneResolution(override, source="override")
            self._cache[key] = result
            return result

        result = self._resolve_exact(cleaned, taxid=taxid)
        if result.warning == "gene_not_resolved" and "/" in cleaned:
            result = self._resolve_composite(cleaned, taxid=taxid)
        self._cache[key] = result
        return result

    def _resolve_exact(self, value: str, *, taxid: int) -> GeneResolution:
        if not self.allow_network:
            return GeneResolution(value, "gene_not_resolved")

        services_available = False
        if taxid == 9606:
            hgnc_result, hgnc_available = self._resolve_hgnc(value)
            services_available |= hgnc_available
            if hgnc_result is not None:
                return hgnc_result

        mygene_result, mygene_available = self._resolve_mygene(value, taxid=taxid)
        services_available |= mygene_available
        if mygene_result is not None:
            return mygene_result

        warning = "gene_not_resolved" if services_available else "gene_service_unavailable"
        return GeneResolution(value, warning)

    def _resolve_composite(self, value: str, *, taxid: int) -> GeneResolution:
        parts = [part.strip() for part in value.split("/") if part.strip()]
        resolved = [self.resolve(part, taxid=taxid) for part in parts]
        symbols = {
            item.symbol
            for item in resolved
            if item.warning not in {"gene_not_resolved", "gene_service_unavailable", "gene_ambiguous"}
        }
        if len(symbols) == 1:
            warning = None if all(item.warning is None for item in resolved) else "gene_partially_resolved"
            return GeneResolution(symbols.pop(), warning, "composite")
        if len(symbols) > 1:
            return GeneResolution(value, "gene_ambiguous")
        return GeneResolution(value, "gene_not_resolved")

    def _get_json(self, url: str, *, accept: str = "application/json") -> Optional[dict]:
        request = Request(url, headers={"Accept": accept, "User-Agent": "sirna-schema/1.0"})
        try:
            with urlopen(request, timeout=self.timeout_seconds) as response:
                return json.load(response)
        except (OSError, ValueError, json.JSONDecodeError):
            return None

    def _resolve_hgnc(self, value: str) -> tuple[Optional[GeneResolution], bool]:
        available = False
        candidates: set[str] = set()
        for field in ("symbol", "alias_symbol", "prev_symbol"):
            search_url = f"https://rest.genenames.org/search/{field}/{quote(value, safe='')}"
            payload = self._get_json(search_url)
            if payload is None:
                # Treat the first transport failure as a service outage instead
                # of multiplying the timeout across every candidate field.
                return None, available
            available = True
            for hit in payload.get("response", {}).get("docs", []):
                hgnc_id = hit.get("hgnc_id")
                if not hgnc_id:
                    continue
                fetch_url = f"https://rest.genenames.org/fetch/hgnc_id/{quote(hgnc_id, safe='')}"
                record_payload = self._get_json(fetch_url)
                if record_payload is None:
                    continue
                docs = record_payload.get("response", {}).get("docs", [])
                if not docs:
                    continue
                record = docs[0]
                if record.get("status") != "Approved":
                    continue
                field_value = record.get(field)
                values = field_value if isinstance(field_value, list) else [field_value]
                if any(isinstance(item, str) and item.casefold() == value.casefold() for item in values):
                    candidates.add(record["symbol"])
        if len(candidates) == 1:
            return GeneResolution(candidates.pop(), source="HGNC"), available
        if len(candidates) > 1:
            return GeneResolution(value, "gene_ambiguous", "HGNC"), available
        return None, available

    def _resolve_mygene(self, value: str, *, taxid: int) -> tuple[Optional[GeneResolution], bool]:
        available = False
        candidates: set[str] = set()
        for field in ("symbol", "alias"):
            query = f'{field}:"{value}"'
            params = urlencode({
                "q": query,
                "species": str(taxid),
                "fields": "symbol,alias,taxid,entrezgene",
                "size": "10",
            })
            payload = self._get_json(f"https://mygene.info/v3/query?{params}")
            if payload is None:
                return None, available
            available = True
            for hit in payload.get("hits", []):
                if hit.get("taxid") != taxid or not hit.get("symbol"):
                    continue
                field_value = hit.get(field)
                values = field_value if isinstance(field_value, list) else [field_value]
                if any(isinstance(item, str) and item.casefold() == value.casefold() for item in values):
                    candidates.add(hit["symbol"])
        if len(candidates) == 1:
            return GeneResolution(candidates.pop(), source="MyGene.info"), available
        if len(candidates) > 1:
            return GeneResolution(value, "gene_ambiguous", "MyGene.info"), available
        return None, available


_DEFAULT_GENE_RESOLVER = GeneSymbolResolver()


# ---------------------------------------------------------------------------
# Modification
# ---------------------------------------------------------------------------

class Modification(BaseModel):
    """One chemical modification at one position on one strand."""
    model_config = {"extra": "forbid"}

    type: Literal[
        "2'-OMe", "2'-Fluoro", "2'-MOE", "2'-Amino",
        "phosphorothioate", "GalNAc", "LNA", "dT-overhang",
        "cholesterol", "other"
    ] = Field(
        description=(
            "Decoded modification chemistry. Use abbrev_map first, then "
            "modification_decoding_reference.md defaults. "
            "Use 'other' only when the type cannot be determined."
        )
    )
    position: Optional[str] = Field(
        default=None,
        description=(
            "1-indexed position from the 5' end of this strand, e.g. '7' or '1-2'. "
            "Omit if the patent does not state a position."
        )
    )
    strand: Literal["sense", "antisense"] = Field(
        description="Which strand this modification is on."
    )
    raw_symbol: Optional[str] = Field(
        default=None,
        description=(
            "The original shorthand exactly as it appears in the source "
            "(e.g. 'Am', 'fC', 's'). Populate when the symbol was ambiguous "
            "or could not be fully decoded; leave empty when decoding was certain."
        )
    )
    flagged: bool = Field(
        default=False,
        description=(
            "True if this modification symbol was ambiguous or could not be "
            "confidently decoded. When True, also set needs_review=True on the "
            "parent SiRNARecord."
        )
    )

    def to_str(self) -> str:
        """Serialise to one semicolon-list element: e.g. \"2'-OMe@7(Am)\" or \"phosphorothioate@20!\"."""
        pos = f"@{self.position}" if self.position else ""
        sym = f"({self.raw_symbol})" if self.raw_symbol else ""
        flag = "!" if self.flagged else ""
        return f"{self.type}{pos}{sym}{flag}"


# ---------------------------------------------------------------------------
# Activity
# ---------------------------------------------------------------------------

class Activity(BaseModel):
    """One activity measurement for a siRNA strand."""
    model_config = {"extra": "forbid"}
    _normalization_warnings: set[str] = PrivateAttr(default_factory=set)

    metric: Optional[Literal[
        "IC50", "EC50", "knockdown_percent", "mRNA_remaining_percent",
        "ED50", "cell_viability_percent", "KSP_GAPDH_ratio", "other"
    ]] = Field(
        default=None,
        description=(
            "What was measured. "
            "'mRNA_remaining_percent': % of target mRNA remaining after treatment (lower = more knockdown). "
            "'knockdown_percent': direct % silenced. "
            "'IC50'/'EC50': concentration for 50% effect, with units in value. "
            "'cell_viability_percent': cell survival assay, not knockdown. "
            "'KSP_GAPDH_ratio': in-vivo ratio of target to housekeeping gene in tissue. "
            "'other': a target-silencing magnitude that none of the above fits -- a "
            "fold-change, a log reduction, a ratio to control. It still needs a number "
            "in value. It is NOT a free-text field: a duration ('18 weeks'), a serum "
            "stability half-life ('t1/2 = 32 days'), or a conclusion with no number "
            "('effective at inhibiting angiogenesis') is not an activity measurement. "
            "Put those in notes and leave metric empty -- the row is still a valid "
            "sequence record. Validation demotes them for you and sets needs_review. "
            "Leave metric empty rather than naming one whose number you do not have: "
            "a metric with an empty value reads as a measurement that was made and "
            "lost, and is flagged as an unfinished transcription."
        )
    )
    value: Optional[str] = Field(
        default=None,
        description=(
            "Canonical numeric result. Percentage metrics require an explicit '%' "
            "and preserve decimal precision, e.g. '85%' or '0.029%'; "
            "percentage ranges, inequalities, approximations, and uncertainty are invalid. "
            "IC50/EC50 concentrations are converted to nM and may retain a reported "
            "range, inequality, approximation, or uncertainty."
        )
    )
    value_numeric: Optional[float] = Field(
        default=None,
        description=(
            "Machine-readable number parsed from value: the magnitude without the unit, "
            "the percent sign, or the inequality. 85% -> 85.0; <0.05 nM -> 0.05; "
            "1.2E-09 -> 1.2e-09. None when value is absent or has no parseable number."
        )
    )
    unit: Optional[str] = Field(
        default=None,
        description=(
            "Unit of value_numeric: '%' for the percentage metrics, a concentration unit "
            "(nM, uM, pM, ...) for IC50/EC50, or empty when the source states none."
        )
    )
    assay_type: Optional[str] = Field(
        default=None,
        description=(
            "Experimental method used to produce the measurement. "
            "Examples: 'qPCR', 'branched-DNA', 'luciferase reporter', "
            "'single_dose_screen', 'dose_response_IC20_IC50_IC80'."
        )
    )
    concentration: Optional[str] = Field(
        default=None,
        description=(
            "Canonical siRNA dose used in the in-vitro experiment. Supported molar "
            "units are converted to nM; ranges and inequalities are retained. A plain "
            "number is accepted but causes needs_review unless a default unit is supplied "
            "through validation context. Leave empty for in-vivo rows."
        )
    )
    cell_line: Optional[str] = Field(
        default=None,
        description=(
            "Cell line where in-vitro activity was measured, e.g. 'HeLa', 'Hep3B', 'A549'. "
            "Leave empty for in-vivo rows."
        )
    )
    in_vivo: bool = Field(
        default=False,
        description=(
            "True if this measurement came from an in-vivo experiment. "
            "When True, populate species, route, dose_mg_per_kg, tissue, timepoint."
        )
    )
    species: Optional[str] = Field(
        default=None,
        description="Species/strain for in-vivo experiments, e.g. 'mouse', 'rat', 'NHP', 'cynomolgus'."
    )
    route: Optional[str] = Field(
        default=None,
        description=(
            "Route of administration for in-vivo experiments. "
            "Examples: 'IV', 'SC', 'GalNAc-SC', 'IT', 'intravitreal'. "
            "'GalNAc-SC' distinguishes subcutaneous GalNAc-conjugate dosing from plain SC."
        )
    )
    dose_mg_per_kg: Optional[str] = Field(
        default=None,
        description=(
            "Dose in mg/kg for in-vivo experiments, e.g. '10', '0.3'. "
            "Extract from group labels or table headers if not in a dedicated column."
        )
    )
    tissue: Optional[str] = Field(
        default=None,
        description=(
            "Tissue or biofluid sampled for in-vivo measurement. "
            "Examples: 'liver', 'serum', 'lung', 'kidney'."
        )
    )
    timepoint: Optional[str] = Field(
        default=None,
        description=(
            "Time post-dose when the measurement was taken. "
            "Examples: 'day 27', '72h', '4 weeks'. "
            "Extract from group labels or study design text if not in a dedicated column."
        )
    )
    notes: Optional[str] = Field(
        default=None,
        description=(
            "Overflow field for anything that does not fit the structured columns above. "
            "Use semicolons to separate multiple items. "
            "Examples: secondary IC50 values ('IC20=2.3nM;IC80=31nM'), "
            "cross-cell-line results ('SKOV3_IC50=142nM;HCT116_IC50=31.6nM'), "
            "comparison to control ('vs AD-1955 luciferase control'), "
            "or a note that data lives in a figure ('see Fig.5, page 85')."
        )
    )

    @field_validator("value", "concentration", mode="before")
    @classmethod
    def coerce_measurement_to_text(cls, value: Any) -> Optional[str]:
        return _clean_text(value)

    @model_validator(mode="after")
    def normalize_measurements(self, info: ValidationInfo) -> Activity:
        """Store canonical strings in the existing fields; do not add CSV columns."""
        warnings: set[str] = set()
        context = info.context or {}
        default_unit = context.get("concentration_default_unit")

        if self.concentration is not None:
            normalized = _normalise_quantity(
                self.concentration,
                quantity_type="concentration",
                default_unit=default_unit,
            )
            object.__setattr__(self, "concentration", normalized.value or None)
            if normalized.warning:
                warnings.add(normalized.warning)

        # `other` has no unit and no shape of its own, so it is the only metric
        # a non-measurement can hide in. Demote instead of rejecting: the row is
        # still a real strand record, the text is kept verbatim in notes, and
        # the metric column stops asserting a measurement that was never made.
        # Raising here would abort the whole CSV write and lose the patent's
        # good rows along with the bad one -- 1,127 rows in WO2005007855's case.
        if self.metric == "other":
            kind = _non_measurement_kind(self.value)
            if kind is not None:
                if self.value is not None:
                    note = f"non_activity_readout[{kind}]: {self.value}"
                    object.__setattr__(
                        self, "notes", f"{self.notes}; {note}" if self.notes else note)
                object.__setattr__(self, "value", None)
                object.__setattr__(self, "metric", None)
                warnings.add(f"activity_not_a_measurement_{kind}")
        elif self.metric is not None and self.value is None:
            # A named metric is a specific claim -- this duplex has an IC50,
            # this one was measured as percent knockdown -- so an empty value is
            # a transcription that did not finish, not a row to reclassify.
            # Keep the metric and make the gap visible instead of hiding it.
            warnings.add("activity_metric_without_value")

        if self.value is not None:
            if self.metric in _PERCENT_METRICS:
                normalized = _normalise_percentage(self.value, remaining=self.metric == "mRNA_remaining_percent")
                object.__setattr__(self, "value", normalized.value or None)
                if normalized.warning:
                    warnings.add(normalized.warning)
            elif self.metric in _CONCENTRATION_METRICS:
                normalized = _normalise_quantity(
                    self.value,
                    quantity_type="concentration",
                    default_unit=default_unit,
                )
                object.__setattr__(self, "value", normalized.value or None)
                if normalized.warning:
                    warnings.add(normalized.warning)
            else:
                object.__setattr__(self, "value", _clean_text(self.value))

        # Union, not assignment. Pydantic runs this validator again when the
        # Activity is nested into a SiRNARecord, and the second pass sees the
        # already-normalised object: a demoted row now has metric=None and
        # value=None, so nothing re-derives the warning and a plain assignment
        # would silently drop it -- leaving needs_review False on exactly the
        # rows that need review. Normalisations that leave their evidence in
        # place (a unitless concentration) re-warn either way; this one cannot.
        self._normalization_warnings = warnings | self._normalization_warnings
        return self

    @model_validator(mode="after")
    def derive_numeric_value(self) -> Activity:
        """Split value into a number and a unit; never raise, never guess a unit."""
        import re as _re

        value = self.value
        if not value:
            object.__setattr__(self, "value_numeric", None)
            object.__setattr__(self, "unit", None)
            return self
        text = str(value).strip()
        match = _re.search(r"[<>~≈≤≥]?\s*(-?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?)\s*%?", text)
        numeric = None
        if match:
            try:
                numeric = float(match.group(1).replace(",", "."))
            except ValueError:
                numeric = None
        unit = None
        if self.metric in _PERCENT_METRICS or "%" in text:
            unit = "%"
        else:
            unit_match = _re.search(r"(nM|µM|μM|uM|mM|pM|ng/mL|mg/kg|mg\s*/\s*kg)", text, _re.I)
            unit = unit_match.group(1) if unit_match else None
        object.__setattr__(self, "value_numeric", numeric)
        object.__setattr__(self, "unit", unit)
        return self


# ---------------------------------------------------------------------------
# SiRNARecord
# ---------------------------------------------------------------------------

class SiRNARecord(BaseModel):
    """
    One extracted strand record — sense OR antisense.
    The agent builds two per duplex row (one per strand).

    Construction validates all fields. extra="forbid" rejects any key not
    declared below, so the agent cannot accidentally invent columns.
    Call write_activity_csv(records, path) to flatten and write the CSV.
    """
    model_config = {"extra": "forbid"}
    _normalization_warnings: set[str] = PrivateAttr(default_factory=set)

    # --- Identity ---
    sirna_id: str = Field(
        description=(
            "Duplex identifier exactly as it appears in the patent, "
            "e.g. 'AD-12115', 'AL-DP-4440'. "
            "Preserve hyphens and case as printed — do not normalise. "
            "If the table has no duplex ID column, synthesise one from table "
            "label and row index, e.g. 'Table6A_row_001', and set needs_review=True."
        )
    )
    strand: Literal["sense", "antisense"] = Field(
        description="Which strand this record represents."
    )
    paired_strand_id: Optional[str] = Field(
        default=None,
        description=(
            "sirna_id of the complementary strand record. "
            "Allows the pair to be rejoined after CSV export. "
            "Should be identical to sirna_id since both strands share the same duplex ID."
        )
    )
    target_gene: Optional[str] = Field(
        default=None,
        description=(
            "Canonical approved symbol for the gene or mRNA being silenced. During "
            "validation, exact symbols, aliases, and previous symbols are resolved using "
            "HGNC for human genes and MyGene.info as a species-scoped fallback. Ambiguous "
            "terms remain unchanged and set needs_review; use a patent-specific resolver "
            "override only when the source explicitly defines the abbreviation."
        )
    )

    @field_validator("target_gene", mode="before")
    @classmethod
    def clean_target_gene(cls, value: Any) -> Optional[str]:
        return _clean_text(value)

    # --- Provenance ---
    source_table: Optional[str] = Field(
        default=None,
        description=(
            "Table or figure label as it appears in the patent, "
            "e.g. 'Table 3', 'Example 13'. "
            "Preserve verbatim including any OCR artefacts in the title "
            "(e.g. 'EGNL1' instead of 'EGLN1') — note the discrepancy in the summary."
        )
    )
    page_index: Optional[int] = Field(
        default=None,
        description="page_index from the page cache for the page this record was extracted from."
    )
    source_quote: Optional[str] = Field(
        default=None,
        description=(
            "Verbatim table row or sentence from which this record was drawn. "
            "Primary audit trail — always populate. Not written to CSV; kept in records.json."
        )
    )

    # --- Sequence ---
    seq_id: Optional[str] = Field(
        default=None,
        description=(
            "SEQ ID NO as cited in the patent table for this strand, e.g. '1240'. "
            "Used to join against seq_cache. Leave empty if no SEQ ID is given."
        )
    )
    sequence: Optional[str] = Field(
        default=None,
        description=(
            "Strand sequence 5'→3', shorthand preserved exactly as in source. "
            "Populated from seq_cache (sequence.txt ground truth) when the SEQ ID resolves; "
            "otherwise from the OCR'd table verbatim. "
            "Never auto-correct or infer missing bases."
        )
    )
    sequence_source: Literal["sequence_txt", "fulltext", "ocr"] = Field(
        description=(
            "How the bases were read. Two of the three are machine-readable, one is not. "
            "'sequence_txt': joined from the ST.25/ST.26 ground-truth listing. "
            "'fulltext': read from machine-readable publication full text, e.g. an EPO OPS "
            "description, so the characters are exact even though they were not cross-checked "
            "against a sequence listing (for many such patents there is none). "
            "'ocr': came from an OCR'd page, so the characters themselves may be wrong; also "
            "used when a listing exists but the SEQ ID was not found or is outside the cached "
            "range. Only this value implies the bases may be misread."
        )
    )

    # --- Modifications ---
    modifications: list[ModificationGroup] = Field(
        default_factory=list,
        description=(
            "Chemical modifications on this strand, grouped by canonical chemistry: "
            "one group per canonical_name with the 1-based positions it covers and the "
            "printed symbol at each position. canonical_name must be a name from "
            "references/modification_vocabulary.json (resolve with abbrev_map first, then "
            "modification_decoding_reference.md). Flattened to a ' * '-joined string in "
            "the CSV by to_csv_row()."
        )
    )

    @field_validator("modifications", mode="before")
    @classmethod
    def _group_legacy_modifications(cls, value):
        """Accept the old per-position entries and fold them into groups.

        Records written before the grouped format carry one dict per
        type/position (type, position, strand, raw_symbol, flagged). They must
        keep loading, so the grouping happens at validation time instead of a
        separate migration step.
        """
        if value is None:
            return []
        if not isinstance(value, list):
            return value
        legacy = [e for e in value if isinstance(e, dict) and "canonical_name" not in e]
        if not legacy:
            return value
        grouped: dict[str, dict] = {}
        for entry in legacy:
            name = str(entry.get("type") or "other")
            position = entry.get("position")
            positions = []
            text = str(position or "").strip()
            if text.isdigit():
                positions = [int(text)]
            else:
                import re as _re
                m = _re.fullmatch(r"(\d+)\s*-\s*(\d+)", text)
                if m:
                    a, b = int(m.group(1)), int(m.group(2))
                    positions = list(range(min(a, b), max(a, b) + 1))
            group = grouped.setdefault(name, {"canonical_name": name, "positions": [], "raw_names": [], "flagged": False})
            for pos in positions:
                if pos < 1:  # 0 means "no position stated"; the grouped form needs 1-based
                    continue
                group["positions"].append(pos)
                group["raw_names"].append(str(entry.get("raw_symbol") or ""))
            if entry.get("flagged"):
                group["flagged"] = True
        merged = [g for g in value if isinstance(g, dict) and "canonical_name" in g]
        merged += [g for g in grouped.values() if g["positions"]]
        return merged

    # --- Activity ---
    activity: Optional[Activity] = Field(
        default=None,
        description=(
            "Primary reported activity measurement for this siRNA. "
            "If the patent reports multiple measurements (e.g. multiple doses, "
            "multiple cell lines), capture the most prominent one here and pack "
            "the rest into activity.notes as semicolon-delimited items."
        )
    )

    # --- Flags ---
    needs_review: bool = Field(
        default=False,
        description=(
            "True if any field is uncertain, ambiguous, or absent but expected. "
            "Auto-set to True when sequence_source='ocr', since only OCR can get the bases "
            "themselves wrong. 'fulltext' is exact text and is not auto-flagged, though it is "
            "not listing-verified either. "
            "Also set True when: a modification is flagged, the sequence length "
            "is outside 18-23 nt, the duplex ID was synthesised, or any value "
            "was inferred rather than read directly from the source."
        )
    )

    @model_validator(mode="after")
    def normalize_and_validate_record(self, info: ValidationInfo) -> SiRNARecord:
        warnings: set[str] = set()

        if self.target_gene:
            context = info.context or {}
            resolver = context.get("gene_resolver", _DEFAULT_GENE_RESOLVER)
            taxid = int(context.get("gene_taxid", 9606))
            if resolver is not None:
                resolution = resolver.resolve(self.target_gene, taxid=taxid)
                object.__setattr__(self, "target_gene", resolution.symbol)
                if resolution.warning:
                    warnings.add(resolution.warning)

        if self.sequence_source == "ocr":
            self.needs_review = True
        if any(g.flagged for g in self.modifications):
            self.needs_review = True
        if self.activity and self.activity._normalization_warnings:
            warnings.update(self.activity._normalization_warnings)
        if self.sequence and re.fullmatch(r"[ACGTUNacgtun]+", self.sequence):
            if not 18 <= len(self.sequence) <= 23:
                warnings.add("sequence_length_outside_18_23")
        if warnings:
            self.needs_review = True
        self._normalization_warnings = warnings
        return self

    @property
    def normalization_warnings(self) -> tuple[str, ...]:
        """Machine-readable review reasons without adding them to JSON or CSV."""
        return tuple(sorted(self._normalization_warnings))

    # ------------------------------------------------------------------
    # Flattening
    # ------------------------------------------------------------------

    def modifications_str(self) -> str:
        """Canonical names joined by ' * ', in group order (dataset standard).

        A strand with no chemistry renders '0' in both the names and the
        positions column, matching the reference dataset's convention.
        """
        if not self.modifications:
            return "0"
        return " * ".join(g.canonical_name for g in self.modifications)

    def modification_positions_str(self) -> str:
        """One comma-separated position list per group, joined by ' * '.

        Groups stay in the same order as modifications_str so the two columns
        read as aligned pairs: 'Phosphorothioate * 2-O-Methyl' against
        '1,2 * 1,2,3,4'. A strand with no chemistry renders '0'.
        """
        if not self.modifications:
            return "0"
        return " * ".join(",".join(str(p) for p in g.positions) for g in self.modifications)

    def to_csv_row(self) -> dict:
        """
        Flatten to the 23-column CSV dict. This is the only flattening path.
        The column order here defines CSV_COLUMNS and the output column order.
        """
        act = self.activity or Activity()
        return {
            "sirna_id":               self.sirna_id,
            "strand":                 self.strand,
            "paired_strand_id":       self.paired_strand_id or "",
            "target_gene":            self.target_gene or "",
            "source_table":           self.source_table or "",
            "page_index":             self.page_index if self.page_index is not None else "",
            "seq_id":                 self.seq_id or "",
            "sequence":               self.sequence or "",
            "sequence_source":        self.sequence_source,
            "modifications":          self.modifications_str(),
            "modification_positions": self.modification_positions_str(),
            "activity_metric":        act.metric or "",
            "activity_value":         act.value or "",
            "activity_value_numeric": "" if act.value_numeric is None else repr(act.value_numeric),
            "activity_unit":          act.unit or "",
            "activity_assay_type":    act.assay_type or "",
            "activity_concentration": act.concentration or "",
            "activity_cell_line":     act.cell_line or "",
            "activity_in_vivo":       str(act.in_vivo),
            "activity_species":       act.species or "",
            "activity_route":         act.route or "",
            "activity_dose_mg_per_kg": act.dose_mg_per_kg or "",
            "activity_tissue":        act.tissue or "",
            "activity_timepoint":     act.timepoint or "",
            "activity_notes":         act.notes or "",
            "needs_review":           str(self.needs_review),
        }


# ---------------------------------------------------------------------------
# CSV column contract — derived from to_csv_row(), never declared separately
# ---------------------------------------------------------------------------

# The column contract is built after ModificationGroup is defined, further down.


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def write_activity_csv(
    records: list[SiRNARecord | dict],
    path: str | Path,
    *,
    gene_resolver: Optional[GeneSymbolResolver] = _DEFAULT_GENE_RESOLVER,
    gene_taxid: int = 9606,
    concentration_default_unit: Optional[str] = None,
) -> None:
    """
    Validate and write the activity CSV.

    Accepts SiRNARecord instances or raw dicts. Raw dictionaries are validated
    with normalisation context before flattening. ``gene_taxid`` defaults to
    human; pass a different NCBI taxonomy ID for non-human target genes.

    Raises ValidationError listing all failures before writing any output.
    The file is created only if every record passes.
    """
    from pydantic import ValidationError

    errors: list[str] = []
    rows: list[dict] = []
    context = {
        "gene_resolver": gene_resolver,
        "gene_taxid": gene_taxid,
        "concentration_default_unit": concentration_default_unit,
    }

    for i, rec in enumerate(records):
        try:
            if not isinstance(rec, SiRNARecord):
                rec = SiRNARecord.model_validate(rec, context=context)
            rows.append(rec.to_csv_row())
        except (ValidationError, TypeError) as exc:
            sid = rec.get("sirna_id", "?") if isinstance(rec, dict) else "?"
            errors.append(f"  Row {i} ({sid}): {exc}")

    if errors:
        raise ValueError(
            f"Validation failed — {len(errors)} error(s):\n" + "\n".join(errors)
        )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows × {len(CSV_COLUMNS)} columns → {path}")


# ---------------------------------------------------------------------------
# Sequence–Tm output contract
# ---------------------------------------------------------------------------

TM_CSV_COLUMNS = [
    "patent_id_origin",
    "sense_sequence",
    "antisense_sequence",
    "modification_sense",
    "sense_position",
    "modification_antisense",
    "antisense_position",
    "tm_value",
]

_SCALAR_TM_RE = re.compile(
    r"^[+-]?(?:\d+(?:[.,]\d+)?|\.\d+)(?:[eE][+-]?\d+)?$"
)
_POSITION_GROUP_RE = re.compile(r"^[1-9]\d*(?:,[1-9]\d*)*$")


class ModificationGroup(BaseModel):
    """One canonical chemistry applied at one or more nucleotide positions."""

    model_config = {"extra": "forbid"}

    canonical_name: str
    positions: tuple[int, ...]
    raw_names: tuple[str, ...] = ()
    flagged: bool = False
    raw_symbols: tuple[str, ...] = ()

    @field_validator("canonical_name", mode="before")
    @classmethod
    def clean_canonical_name(cls, value: Any) -> str:
        text = _clean_text(value)
        if text is None:
            raise ValueError("canonical modification name is required")
        return text

    @field_validator("positions", mode="before")
    @classmethod
    def normalize_positions(cls, value: Any) -> tuple[int, ...]:
        if isinstance(value, str):
            value = [item.strip() for item in value.split(",") if item.strip()]
        positions = tuple(sorted({int(item) for item in value}))
        if not positions or positions[0] < 1:
            raise ValueError("modification positions must be positive one-indexed integers")
        return positions


class UnresolvedChemistry(BaseModel):
    """Raw chemistry evidence retained until the agent resolves it."""

    model_config = {"extra": "forbid"}

    raw_name: str
    raw_symbol: Optional[str] = None
    positions: tuple[int, ...]
    patent_definition: Optional[str] = None
    sequence_ocr: str
    source_quote: Optional[str] = None
    page_index: Optional[int] = None
    table_id: Optional[str] = None
    occurrences: int = 1


class ChemistryResolution(BaseModel):
    """Auditable agent decision for one previously unresolved chemistry."""

    model_config = {"extra": "forbid"}

    raw_name: str
    canonical_name: str
    resolution: Literal["existing_canonical", "canonical_alias", "new_canonical"]
    evidence: str
    confidence: Literal["high", "medium", "low"]
    needs_review: bool = False


# ModificationGroup is declared after SiRNARecord, so the forward reference is
# resolved here, before the CSV contract instantiates the model.
SiRNARecord.model_rebuild()

CSV_COLUMNS: list[str] = list(
    SiRNARecord(sirna_id="_", strand="sense", sequence_source="ocr").to_csv_row().keys()
)
assert len(CSV_COLUMNS) == 26, f"Column count drifted: {len(CSV_COLUMNS)}"


def serialize_modification_groups(groups: list[ModificationGroup]) -> tuple[str, str]:
    """Serialize groups into the training dataset's aligned pair of columns."""
    if not groups:
        return "0", "0"
    ordered = sorted(
        groups,
        key=lambda group: (min(group.positions), group.canonical_name.casefold()),
    )
    modification = " * ".join(group.canonical_name for group in ordered)
    position = " * ".join(
        ",".join(str(item) for item in group.positions) for group in ordered
    )
    return modification, position


def _normalise_tm_scalar(value: Any) -> Decimal:
    """Return a finite scalar Celsius value quantized to two decimals."""
    if isinstance(value, bool):
        raise ValueError("Tm must be a numeric scalar, not a boolean")
    if isinstance(value, Decimal):
        result = value
    else:
        text = _clean_text(value)
        if text is None:
            raise ValueError("Tm is required")
        text = re.sub(r"(?i)\s*(?:°\s*C|deg(?:ree)?s?\s*C|Celsius)\s*$", "", text).strip()
        if not _SCALAR_TM_RE.fullmatch(text):
            raise ValueError(
                "Tm must be one scalar number; inequalities, ranges, approximations, "
                "uncertainty, and missing-value labels are not allowed"
            )
        result = _decimal(text)
    if not result.is_finite():
        raise ValueError("Tm must be finite")
    return result.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)


class SequenceTmOutputRecord(BaseModel):
    """One merge-safe duplex/Tm row validated before CSV serialization."""

    model_config = {"extra": "forbid"}

    patent_id_origin: str
    sense_sequence: str
    antisense_sequence: str
    modification_sense: str = "0"
    sense_position: str = "0"
    modification_antisense: str = "0"
    antisense_position: str = "0"
    tm_value: Decimal

    @field_validator("patent_id_origin", mode="before")
    @classmethod
    def normalize_patent_id(cls, value: Any) -> str:
        text = _clean_text(value)
        if text is None:
            raise ValueError("patent_id_origin is required")
        normalized = re.sub(r"[^A-Za-z0-9]", "", text).upper()
        if not normalized:
            raise ValueError("patent_id_origin contains no letters or digits")
        return normalized

    @field_validator("sense_sequence", "antisense_sequence", mode="before")
    @classmethod
    def normalize_sequence(cls, value: Any) -> str:
        text = _clean_text(value)
        if text is None:
            raise ValueError("sequence is required")
        # The reference dataset sometimes appends a terminal inverted-abasic
        # label to the nucleotide string while representing the chemistry in
        # the parallel modification/position fields. Keep only the bases here.
        text = re.sub(r"(?i)\((?:invAb|iAb|inverted\s+abasic)\)\s*$", "", text)
        normalized = re.sub(r"\s+", "", text).upper()
        if not re.fullmatch(r"[ACGTUN]+", normalized):
            raise ValueError("sequence must contain only A/C/G/T/U/N after whitespace removal")
        return normalized

    @field_validator(
        "modification_sense",
        "sense_position",
        "modification_antisense",
        "antisense_position",
        mode="before",
    )
    @classmethod
    def normalize_group_text(cls, value: Any) -> str:
        text = _clean_text(value)
        if text is None:
            raise ValueError("modification and position are required; use '0' when unmodified")
        parts = [part.strip() for part in text.split("*")]
        if any(not part for part in parts):
            raise ValueError("modification and position groups cannot be empty")
        return " * ".join(parts)

    @field_validator("tm_value", mode="before")
    @classmethod
    def normalize_tm(cls, value: Any) -> Decimal:
        return _normalise_tm_scalar(value)

    @field_serializer("tm_value")
    def serialize_tm(self, value: Decimal) -> str:
        return format(value, ".2f")

    @model_validator(mode="after")
    def validate_modification_positions(self) -> SequenceTmOutputRecord:
        for strand, sequence, modification, position in (
            ("sense", self.sense_sequence, self.modification_sense, self.sense_position),
            ("antisense", self.antisense_sequence, self.modification_antisense, self.antisense_position),
        ):
            if modification == "0" or position == "0":
                if modification != "0" or position != "0":
                    raise ValueError(f"unmodified {strand} must use modification='0' and position='0'")
                continue
            modifications = [part.strip() for part in modification.split(" * ")]
            position_groups = [part.strip() for part in position.split(" * ")]
            if len(modifications) != len(position_groups):
                raise ValueError(f"{strand} modification and position group counts must match")
            normalized_groups: list[str] = []
            for group in position_groups:
                if not _POSITION_GROUP_RE.fullmatch(group):
                    raise ValueError(f"invalid {strand} modification position group: {group!r}")
                positions = sorted({int(item) for item in group.split(",")})
                if positions[-1] > len(sequence) + 2:
                    raise ValueError(
                        f"{strand} modification position {positions[-1]} exceeds the allowed "
                        f"two-position terminal extension for sequence length {len(sequence)}"
                    )
                normalized_groups.append(",".join(str(item) for item in positions))
            setattr(self, f"modification_{strand}", " * ".join(modifications))
            setattr(self, f"{strand}_position", " * ".join(normalized_groups))
        return self

    def to_csv_row(self) -> dict[str, str]:
        return {
            "patent_id_origin": self.patent_id_origin,
            "sense_sequence": self.sense_sequence,
            "antisense_sequence": self.antisense_sequence,
            "modification_sense": self.modification_sense,
            "sense_position": self.sense_position,
            "modification_antisense": self.modification_antisense,
            "antisense_position": self.antisense_position,
            "tm_value": format(self.tm_value, ".2f"),
        }


def _validate_tm_records(
    records: list[SequenceTmOutputRecord | Mapping[str, Any]],
) -> list[SequenceTmOutputRecord]:
    from pydantic import ValidationError

    validated: list[SequenceTmOutputRecord] = []
    errors: list[str] = []
    for index, record in enumerate(records):
        try:
            if isinstance(record, SequenceTmOutputRecord):
                # Revalidate instances so mutated or constructed data cannot bypass checks.
                record = SequenceTmOutputRecord.model_validate(record.model_dump())
            else:
                record = SequenceTmOutputRecord.model_validate(record)
            validated.append(record)
        except (ValidationError, TypeError, ValueError) as exc:
            errors.append(f"  Row {index}: {exc}")
    if errors:
        raise ValueError(
            f"Sequence–Tm validation failed — {len(errors)} error(s):\n" + "\n".join(errors)
        )
    return validated


def write_sequence_tm_csv(
    records: list[SequenceTmOutputRecord | Mapping[str, Any]],
    path: str | Path,
) -> None:
    """Validate a complete batch, then atomically write the five-column CSV."""
    validated = _validate_tm_records(records)
    rows = [record.to_csv_row() for record in validated]

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            newline="",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=TM_CSV_COLUMNS, extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, destination)
    except Exception:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
        raise


def merge_sequence_tm_csv(inputs: list[str | Path], output: str | Path) -> None:
    """Validate, exactly deduplicate, sort, and atomically merge Tm CSV files."""
    raw_records: list[dict[str, str]] = []
    for input_path in inputs:
        with open(input_path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != TM_CSV_COLUMNS:
                raise ValueError(
                    f"{input_path}: expected header {TM_CSV_COLUMNS}, got {reader.fieldnames}"
                )
            raw_records.extend(dict(row) for row in reader)

    validated = _validate_tm_records(raw_records)
    unique = {tuple(record.to_csv_row()[column] for column in TM_CSV_COLUMNS): record for record in validated}
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
    write_sequence_tm_csv(ordered, output)

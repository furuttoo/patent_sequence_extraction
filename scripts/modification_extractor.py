"""Tokenize modified oligonucleotide notation into canonical grouped chemistry.

Known chemistry is normalized deterministically with the dataset vocabulary.
Unknown chemistry is retained as structured evidence and must be resolved by an
agent-authored chemistry map before a final output row is emitted.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from pathlib import Path
from typing import Any, Mapping, Optional
import unicodedata

from schemas_pydantic import (
    ChemistryResolution,
    ModificationGroup,
    UnresolvedChemistry,
)


VOCABULARY_PATH = Path(__file__).resolve().parents[1] / "references" / "modification_vocabulary.json"
_BASES = "ACGUTN"


def _key(value: str) -> str:
    text = unicodedata.normalize("NFKC", value).replace("′", "'").replace("’", "'")
    return re.sub(r"[\s_-]+", "", text).casefold()


class ChemistryVocabulary:
    def __init__(self, path: str | Path = VOCABULARY_PATH) -> None:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
        self.canonical: set[str] = set()
        self.aliases: dict[str, str] = {}
        for entry in payload["chemistries"]:
            canonical = entry["canonical_name"]
            self.canonical.add(canonical)
            self.aliases[_key(canonical)] = canonical
            for alias in entry.get("aliases", []):
                self.aliases[_key(alias)] = canonical

    def resolve(self, value: str) -> Optional[str]:
        return self.aliases.get(_key(value))


@dataclass
class ModificationExtraction:
    base_sequence: str
    groups: list[ModificationGroup]
    unresolved: list[UnresolvedChemistry]


def _add(
    grouped: dict[str, dict[str, Any]],
    canonical: str,
    position: int,
    raw_name: str,
    raw_symbol: str,
) -> None:
    entry = grouped.setdefault(
        canonical,
        {"positions": set(), "raw_names": set(), "raw_symbols": set()},
    )
    entry["positions"].add(position)
    entry["raw_names"].add(raw_name)
    entry["raw_symbols"].add(raw_symbol)


def _resolve_or_unresolved(
    *,
    raw_name: str,
    raw_symbol: str,
    position: int,
    sequence_ocr: str,
    vocabulary: ChemistryVocabulary,
    resolutions: Mapping[str, str],
    grouped: dict[str, dict[str, Any]],
    unresolved: list[UnresolvedChemistry],
    patent_definition: Optional[str] = None,
) -> None:
    canonical = vocabulary.resolve(raw_name) or resolutions.get(raw_name) or resolutions.get(raw_symbol)
    if canonical:
        _add(grouped, canonical, position, raw_name, raw_symbol)
    else:
        unresolved.append(
            UnresolvedChemistry(
                raw_name=raw_name,
                raw_symbol=raw_symbol,
                positions=(position,),
                patent_definition=patent_definition,
                sequence_ocr=sequence_ocr,
            )
        )


def extract_modifications(
    notation: str,
    *,
    abbrev_map: Optional[Mapping[str, str]] = None,
    resolutions: Optional[Mapping[str, str]] = None,
    lowercase_means_ome: bool = True,
    vocabulary: Optional[ChemistryVocabulary] = None,
) -> ModificationExtraction:
    """Decode common Alnylam-style notation without repairing source text."""
    vocabulary = vocabulary or ChemistryVocabulary()
    abbrev_map = abbrev_map or {}
    resolutions = resolutions or {}
    grouped: dict[str, dict[str, Any]] = {}
    unresolved: list[UnresolvedChemistry] = []
    bases: list[str] = []
    index = 0

    def record(raw_name: str, raw_symbol: str, position: int) -> None:
        definition = abbrev_map.get(raw_symbol)
        name = definition or raw_name
        _resolve_or_unresolved(
            raw_name=name,
            raw_symbol=raw_symbol,
            position=position,
            sequence_ocr=notation,
            vocabulary=vocabulary,
            resolutions=resolutions,
            grouped=grouped,
            unresolved=unresolved,
            patent_definition=definition,
        )

    while index < len(notation):
        char = notation[index]

        latex_modifier = re.match(
            r"\s*\$\s*(?P<base>[ACGUTNacgutn])?\s*_\s*\{+\s*"
            r"(?P<name>[A-Za-z][A-Za-z0-9'-]*)\s*\}+\s*\$",
            notation[index:],
            re.I,
        )
        if latex_modifier:
            annotated_base = latex_modifier.group("base")
            if annotated_base:
                bases.append(annotated_base.upper())
            position = max(1, len(bases))
            record(latex_modifier.group("name"), latex_modifier.group(0).strip(), position)
            index += len(latex_modifier.group(0))
            continue

        if char.isspace() or char in "5'3'′–-":
            index += 1
            continue

        parenthetical_terminal = re.match(
            r"\((?P<name>Chol(?:-TEG)?|invdT|iAb|invAb)\)",
            notation[index:],
            re.I,
        )
        if parenthetical_terminal:
            symbol = parenthetical_terminal.group(0)
            name = parenthetical_terminal.group("name")
            position = max(1, len(bases) + (1 if bases else 0))
            record(name, symbol, position)
            index += len(symbol)
            continue

        terminal = re.match(r"(?:L96|GalNAc|Chol(?:-TEG)?|VP|P)", notation[index:], re.I)
        if terminal:
            symbol = terminal.group(0)
            position = max(1, len(bases))
            record(symbol, symbol, position)
            index += len(symbol)
            continue

        parenthetical = re.match(r"\(([ACGUTNacgutn])([^)]+)\)", notation[index:])
        if parenthetical:
            base, suffix = parenthetical.groups()
            bases.append(base.upper())
            position = len(bases)
            symbol = parenthetical.group(0)
            record(suffix, symbol, position)
            index += len(symbol)
            continue

        deoxy = re.match(r"d([ACGTUacgtu])", notation[index:])
        if deoxy:
            base = deoxy.group(1).upper()
            bases.append(base)
            position = len(bases)
            symbol = deoxy.group(0)
            record("2-Deoxythymidine" if base == "T" else "2-Deoxy", symbol, position)
            index += len(symbol)
            continue

        prefixed = re.match(r"([mf])([ACGUTacgut])", notation[index:])
        if prefixed:
            prefix, base = prefixed.groups()
            bases.append(base.upper())
            position = len(bases)
            record("2-O-Methyl" if prefix.lower() == "m" else "2-F", prefixed.group(0), position)
            index += len(prefixed.group(0))
            continue

        if char.upper() in _BASES:
            bases.append(char.upper())
            position = len(bases)
            raw_symbol = char
            index += 1
            if char.islower() and lowercase_means_ome:
                record("2-O-Methyl", char, position)

            suffixes: list[str] = []
            while index < len(notation) and notation[index] in "mfbs":
                suffixes.append(notation[index])
                raw_symbol += notation[index]
                index += 1
            for suffix in suffixes:
                if suffix == "m":
                    record("2-O-Methyl", raw_symbol, position)
                elif suffix == "f":
                    record("2-F", raw_symbol, position)
                elif suffix == "s":
                    record("Phosphorothioate", raw_symbol, position)
                else:
                    record(suffix, raw_symbol, position)
            continue

        # A phosphorothioate linkage may follow a deoxy token (``dTsdT``),
        # where it is not consumed by the ordinary base-suffix loop.
        if char == "s" and bases:
            record("Phosphorothioate", char, len(bases))
            index += 1
            continue

        unknown = re.match(r"[^ACGUTNacgutn\s]+", notation[index:])
        symbol = unknown.group(0) if unknown else char
        _resolve_or_unresolved(
            raw_name=symbol,
            raw_symbol=symbol,
            position=max(1, len(bases)),
            sequence_ocr=notation,
            vocabulary=vocabulary,
            resolutions=resolutions,
            grouped=grouped,
            unresolved=unresolved,
        )
        index += len(symbol)

    groups = [
        ModificationGroup(
            canonical_name=canonical,
            positions=tuple(details["positions"]),
            raw_names=tuple(sorted(details["raw_names"])),
            raw_symbols=tuple(sorted(details["raw_symbols"])),
        )
        for canonical, details in grouped.items()
    ]
    return ModificationExtraction("".join(bases), groups, unresolved)


def load_resolution_map(path: str | Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    entries = payload.get("resolutions", payload) if isinstance(payload, dict) else payload
    resolutions: dict[str, str] = {}
    if isinstance(entries, dict):
        return {str(key): str(value) for key, value in entries.items()}
    for entry in entries:
        resolution = ChemistryResolution.model_validate(entry)
        resolutions[resolution.raw_name] = resolution.canonical_name
    return resolutions

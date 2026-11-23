from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from modification_extractor import extract_modifications
from schemas_pydantic import serialize_modification_groups


class ModificationExtractorTests(unittest.TestCase):
    def test_suffix_notation(self):
        result = extract_modifications("AmCfGmU")
        self.assertEqual(result.base_sequence, "ACGU")
        self.assertEqual(result.unresolved, [])
        self.assertEqual(
            serialize_modification_groups(result.groups),
            ("2-O-Methyl * 2-F", "1,3 * 2"),
        )

    def test_lowercase_and_phosphorothioate(self):
        result = extract_modifications("csasug")
        self.assertEqual(result.base_sequence, "CAUG")
        self.assertEqual(result.unresolved, [])
        self.assertEqual(
            serialize_modification_groups(result.groups),
            ("2-O-Methyl * Phosphorothioate", "1,2,3,4 * 1,2"),
        )

    def test_deoxy_overhang_and_galnac(self):
        result = extract_modifications("ACGUdTdTL96")
        self.assertEqual(result.base_sequence, "ACGUTT")
        self.assertEqual(result.unresolved, [])
        self.assertEqual(
            serialize_modification_groups(result.groups),
            ("2-Deoxythymidine * N-acetyl-galactosamine", "5,6 * 6"),
        )

    def test_unknown_is_retained(self):
        result = extract_modifications("ACGUZ99")
        self.assertEqual(result.base_sequence, "ACGU")
        self.assertTrue(result.unresolved)
        self.assertEqual(result.unresolved[0].raw_symbol, "Z99")

    def test_resolution_map_applies_to_unknown(self):
        result = extract_modifications(
            "ACGUZ99",
            resolutions={"Z99": "New lipid conjugate"},
        )
        self.assertEqual(result.unresolved, [])
        self.assertEqual(
            serialize_modification_groups(result.groups),
            ("New lipid conjugate", "4"),
        )

    def test_parenthetical_terminal_modifications_are_not_split(self):
        result = extract_modifications(
            "(Chol-TEG)AfuCfaAf(invdT)",
            resolutions={"invdT": "Inverted DeoxyThymidine"},
        )
        self.assertEqual(result.unresolved, [])
        names = {group.canonical_name for group in result.groups}
        self.assertIn("Chol-TEG", names)
        self.assertIn("Inverted DeoxyThymidine", names)

    def test_phosphorothioate_between_deoxy_tokens(self):
        result = extract_modifications("dTsdT")
        self.assertEqual(result.base_sequence, "TT")
        self.assertEqual(result.unresolved, [])
        self.assertEqual(
            serialize_modification_groups(result.groups),
            ("2-Deoxythymidine * Phosphorothioate", "1,2 * 1"),
        )

    def test_latex_una_suffix_and_annotated_base(self):
        suffix = extract_modifications(r"C $ _{UNA} $u")
        annotated = extract_modifications(r"C $ U_{una} $U")
        self.assertEqual(suffix.unresolved, [])
        self.assertEqual(annotated.unresolved, [])
        self.assertEqual(suffix.base_sequence, "CU")
        self.assertEqual(annotated.base_sequence, "CUU")
        self.assertIn("UnLocked nucleic acid", {g.canonical_name for g in suffix.groups})
        self.assertIn("UnLocked nucleic acid", {g.canonical_name for g in annotated.groups})


if __name__ == "__main__":
    unittest.main()

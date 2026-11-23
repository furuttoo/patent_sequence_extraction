import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

from pydantic import ValidationError


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from extract_sequence_tm import aggregate_unresolved_chemistry, extract_from_cache, extract_patent
from schemas_pydantic import (
    SequenceTmOutputRecord,
    merge_sequence_tm_csv,
    write_sequence_tm_csv,
)
from tm_page_filter import build_cache


class SequenceTmSchemaTests(unittest.TestCase):
    def test_scalar_tm_is_two_decimal(self):
        record = SequenceTmOutputRecord.model_validate(
            {
                "patent_id_origin": "WO 2023/014765",
                "sense_sequence": "ACGUACGUACGUACGUACG",
                "antisense_sequence": "UGCAUGCAUGCAUGCAUGC",
                "tm_value": "68.405 °C",
            }
        )
        self.assertEqual(record.patent_id_origin, "WO2023014765")
        self.assertEqual(record.to_csv_row()["tm_value"], "68.41")

    def test_non_scalar_tm_is_rejected(self):
        for value in (">70", "62–64", "65 ± 0.5", "~68", "ND"):
            with self.subTest(value=value), self.assertRaises(ValidationError):
                SequenceTmOutputRecord.model_validate(
                    {
                        "patent_id_origin": "WO1",
                        "sense_sequence": "ACGUACGUACGUACGUACG",
                        "antisense_sequence": "UGCAUGCAUGCAUGCAUGC",
                        "tm_value": value,
                    }
                )

    def test_modification_position_alignment(self):
        record = SequenceTmOutputRecord.model_validate(
            {
                "patent_id_origin": "WO1",
                "sense_sequence": "ACGUACGUACGUACGUACG",
                "antisense_sequence": "UGCAUGCAUGCAUGCAUGC",
                "modification_sense": "Phosphorothioate * 2-O-Methyl",
                "sense_position": "2,1,2 * 4,3",
                "tm_value": 70,
            }
        )
        self.assertEqual(record.sense_position, "1,2 * 3,4")

    def test_empty_modification_group_is_rejected(self):
        with self.assertRaises(ValidationError):
            SequenceTmOutputRecord.model_validate(
                {
                    "patent_id_origin": "WO1",
                    "sense_sequence": "ACGUACGUACGUACGUACG",
                    "antisense_sequence": "UGCAUGCAUGCAUGCAUGC",
                    "modification_sense": "2-O-Methyl ** 2-F",
                    "sense_position": "1,2 * 3,4",
                    "tm_value": "68.4",
                }
            )

    def test_atomic_writer_does_not_replace_on_failure(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "out.csv"
            output.write_text("sentinel", encoding="utf-8")
            with self.assertRaises(ValueError):
                write_sequence_tm_csv(
                    [
                        {
                            "patent_id_origin": "WO1",
                            "sense_sequence": "ACGUACGUACGUACGUACG",
                            "antisense_sequence": "UGCAUGCAUGCAUGCAUGC",
                            "tm_value": ">70",
                        }
                    ],
                    output,
                )
            self.assertEqual(output.read_text(encoding="utf-8"), "sentinel")


class TmFilterAndExtractionTests(unittest.TestCase):
    def _write_json(self, directory: str, name: str, data: object) -> Path:
        path = Path(directory) / name
        path.write_text(json.dumps(data), encoding="utf-8")
        return path

    def test_isolated_modification_tm_does_not_hit(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            patent = self._write_json(
                temp_dir,
                "WO1.json",
                {
                    "patent_id": "WO1",
                    "pages": [
                        {
                            "page_index": 0,
                            "blocks": [
                                {"block_label": "text", "block_content": "Modified sequence AmCmTmGm"}
                            ],
                        }
                    ],
                },
            )
            self.assertTrue(build_cache(patent)["skipped_no_tm"])

    def test_spaced_tm_modification_in_sequence_table_does_not_hit(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            patent = self._write_json(
                temp_dir,
                "WO1.json",
                {
                    "patent_id": "WO1",
                    "pages": [
                        {
                            "page_index": 0,
                            "blocks": [
                                {
                                    "block_label": "table",
                                    "block_content": (
                                        "<table><tr><th>Sequence</th></tr>"
                                        "<tr><td>Am Cm Tm Gm</td></tr></table>"
                                    ),
                                }
                            ],
                        }
                    ],
                },
            )
            self.assertTrue(build_cache(patent)["skipped_no_tm"])

    def test_unresolved_chemistry_is_aggregated(self):
        aggregated = aggregate_unresolved_chemistry(
            [
                {
                    "page_index": 2,
                    "source_row": ["ACGUZ99", "70"],
                    "unresolved": [
                        {
                            "raw_name": "Z99",
                            "raw_symbol": "Z99",
                            "positions": [4],
                            "sequence_ocr": "ACGUZ99",
                        }
                    ],
                },
                {
                    "page_index": 3,
                    "source_row": ["UGCAZ99", "71"],
                    "unresolved": [
                        {
                            "raw_name": "Z99",
                            "raw_symbol": "Z99",
                            "positions": [4],
                            "sequence_ocr": "UGCAZ99",
                        }
                    ],
                },
            ]
        )
        self.assertEqual(len(aggregated), 1)
        self.assertEqual(aggregated[0]["occurrences"], 2)

    def test_table_extracts_scalar_and_rejects_range(self):
        table = """
        <table>
          <tr><th>Antisense sequence</th><th>Sense sequence</th><th>SEQ ID pair</th><th>Tm (°C)</th></tr>
          <tr><td>ACGUACGUACGUACGUACG</td><td>UGCAUGCAUGCAUGCAUGC</td><td>1/2</td><td>68.4</td></tr>
          <tr><td>AAAAACGUACGUACGUACG</td><td>UUUUUGCAUGCAUGCAUGC</td><td>3/4</td><td>62-64</td></tr>
        </table>
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            patent = self._write_json(
                temp_dir,
                "patent.json",
                {
                    "patent_id": "WO 2023/014765",
                    "pages": [
                        {
                            "page_index": 7,
                            "blocks": [{"block_label": "table", "block_content": table}],
                        }
                    ],
                },
            )
            output = Path(temp_dir) / "out.csv"
            rejected = Path(temp_dir) / "rejected.json"
            result = extract_patent(patent, output, rejected_json=rejected)
            self.assertEqual(result["rows"], 1)
            self.assertEqual(result["rejected"], 1)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["patent_id_origin"], "WO2023014765")
            self.assertEqual(rows[0]["tm_value"], "68.40")

    def test_modified_ocr_sequence_is_decoded_before_tm_output(self):
        table = """
        <table>
          <tr><th>Antisense sequence</th><th>Sense sequence</th><th>SEQ ID pair</th><th>Tm (°C)</th></tr>
          <tr><td>AmCfGmU</td><td>UmGfCmA</td><td>12/13</td><td>71</td></tr>
        </table>
        """
        cache = {
            "patent_id": "WO1",
            "pages": [
                {
                    "page_index": 4,
                    "blocks": [{"block_label": "table", "block_content": table}],
                }
            ],
        }
        records, rejected = extract_from_cache(
            cache,
            seq_cache={
                12: {"sequence": "ACGUACGUACGUACGUACG"},
                13: {"sequence": "UGCAUGCAUGCAUGCAUGC"},
            },
        )
        self.assertEqual(rejected, [])
        row = records[0].to_csv_row()
        self.assertEqual(row["antisense_sequence"], "ACGUACGUACGUACGUACG")
        self.assertEqual(row["sense_sequence"], "UGCAUGCAUGCAUGCAUGC")
        self.assertEqual(row["modification_antisense"], "2-O-Methyl * 2-F")
        self.assertEqual(row["antisense_position"], "1,3 * 2")

    def test_sequence_listing_value_wins_for_plain_ocr_sequence(self):
        table = """
        <table>
          <tr><th>Antisense sequence</th><th>Sense sequence</th><th>SEQ ID pair</th><th>Tm (°C)</th></tr>
          <tr><td>ACGUACGUACGUACGUACA</td><td>UGCAUGCAUGCAUGCAUGA</td><td>12/13</td><td>71</td></tr>
        </table>
        """
        cache = {
            "patent_id": "WO1",
            "pages": [
                {
                    "page_index": 4,
                    "blocks": [{"block_label": "table", "block_content": table}],
                }
            ],
        }
        records, rejected = extract_from_cache(
            cache,
            seq_cache={
                12: {"sequence": "ACGUACGUACGUACGUACG"},
                13: {"sequence": "UGCAUGCAUGCAUGCAUGC"},
            },
        )
        self.assertEqual(rejected, [])
        self.assertEqual(records[0].antisense_sequence, "ACGUACGUACGUACGUACG")
        self.assertEqual(records[0].sense_sequence, "UGCAUGCAUGCAUGCAUGC")
        self.assertEqual(records[0].to_csv_row()["tm_value"], "71.00")

    def test_mixed_case_sequence_uses_modification_decoder(self):
        table = """
        <table>
          <tr><th>Antisense sequence</th><th>Sense sequence</th><th>SEQ ID pair</th><th>Tm (°C)</th></tr>
          <tr><td>AcGu</td><td>UgCa</td><td>1/2</td><td>66</td></tr>
        </table>
        """
        cache = {
            "patent_id": "WO1",
            "pages": [
                {
                    "page_index": 1,
                    "blocks": [{"block_label": "table", "block_content": table}],
                }
            ],
        }
        records, rejected = extract_from_cache(cache)
        self.assertEqual(rejected, [])
        self.assertEqual(records[0].antisense_sequence, "ACGU")
        self.assertEqual(records[0].modification_antisense, "2-O-Methyl")
        self.assertEqual(records[0].antisense_position, "2,4")

    def test_latex_tm_header_and_seq_id_pair_side_table_are_joined(self):
        sequence_table = """
        <table>
          <tr><th>SEQ ID No.</th><th>Antisense strand sequence</th>
              <th>SEQ ID No.</th><th>Sense strand sequence</th><th>SEQ ID pair</th></tr>
          <tr><td>50</td><td>ACGU</td><td>62</td><td>UGCA</td><td>50/62</td></tr>
          <tr><td>53</td><td>AAAA</td><td>67</td><td>UUUU</td><td>53/67</td></tr>
        </table>
        """
        tm_table = """
        <table>
          <tr><th>SEQ ID pair</th><th>$ T_{{m}} $</th></tr>
          <tr><td>50/62</td><td>77.1</td></tr>
          <tr><td>53/67</td><td>76.8</td></tr>
        </table>
        """
        cache = {
            "patent_id": "WO1",
            "pages": [{"page_index": 1, "blocks": [
                {"block_label": "table", "block_content": sequence_table},
                {"block_label": "table", "block_content": tm_table},
            ]}],
        }
        records, rejected = extract_from_cache(cache)
        self.assertEqual(rejected, [])
        self.assertEqual(len(records), 2)
        self.assertEqual(
            {(record.antisense_sequence, record.sense_sequence, record.to_csv_row()["tm_value"]) for record in records},
            {("ACGU", "UGCA", "77.10"), ("AAAA", "UUUU", "76.80")},
        )

    def test_merged_tm_and_pair_cell_recovers_mero_components(self):
        table = r"""
        <table>
          <tr><th>Antisense strand sequence</th>
              <th>Sense strand sequence $ T_{m} $</th><th>SEQ ID mero</th></tr>
          <tr><td>ACGU</td><td>UGCA\nAAAA</td><td>46.3 50/64/81</td></tr>
        </table>
        """
        cache = {
            "patent_id": "WO1",
            "pages": [{"page_index": 2, "blocks": [
                {"block_label": "table", "block_content": table},
            ]}],
        }
        records, rejected = extract_from_cache(cache)
        self.assertEqual(rejected, [])
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].antisense_sequence, "ACGU")
        self.assertEqual(records[0].sense_sequence, "UGCAAAAA")
        self.assertEqual(records[0].to_csv_row()["tm_value"], "46.30")

    def test_latex_tm_table_passes_page_gate(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            patent = self._write_json(
                temp_dir,
                "WO1.json",
                {
                    "patent_id": "WO1",
                    "pages": [{
                        "page_index": 0,
                        "blocks": [{
                            "block_label": "table",
                            "block_content": (
                                "<table><tr><th>SEQ ID pair</th><th>$ T_{{m}} $</th></tr>"
                                "<tr><td>50/62</td><td>77.1</td></tr></table>"
                            ),
                        }],
                    }],
                },
            )
            cache = build_cache(patent)
            self.assertFalse(cache["skipped_no_tm"])
            self.assertEqual(cache["tm_hit_count"], 1)

    def test_no_tm_patent_writes_no_csv(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            patent = self._write_json(
                temp_dir,
                "WO2.json",
                {"patent_id": "WO2", "pages": [{"page_index": 0, "blocks": []}]},
            )
            output = Path(temp_dir) / "out.csv"
            result = extract_patent(patent, output)
            self.assertEqual(result["status"], "skipped_no_tm")
            self.assertFalse(output.exists())

    def test_merge_validates_and_deduplicates(self):
        record = {
            "patent_id_origin": "WO1",
            "sense_sequence": "ACGUACGUACGUACGUACG",
            "antisense_sequence": "UGCAUGCAUGCAUGCAUGC",
            "modification_sense": "0",
            "sense_position": "0",
            "modification_antisense": "0",
            "antisense_position": "0",
            "tm_value": "68.40",
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            first = Path(temp_dir) / "first.csv"
            second = Path(temp_dir) / "second.csv"
            merged = Path(temp_dir) / "merged.csv"
            write_sequence_tm_csv([record], first)
            write_sequence_tm_csv([record], second)
            merge_sequence_tm_csv([first, second], merged)
            with merged.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)


if __name__ == "__main__":
    unittest.main()

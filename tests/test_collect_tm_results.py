import json
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from collect_tm_results import (  # noqa: E402
    audit_skipped_patent,
    candidates_from_job_list,
    validate_skip_marker,
)


class CollectTmResultsTests(unittest.TestCase):
    def test_skip_marker_is_terminal_not_pending(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            patent = root / "WO2024000001"
            output = patent / "output"
            output.mkdir(parents=True)
            marker = output / "WO2024000001_tm_skipped.json"
            marker.write_text(
                json.dumps(
                    {
                        "patent_id_origin": "WO2024000001",
                        "status": "skipped_no_tm",
                        "source": "tm_page_filter_preflight",
                    }
                ),
                encoding="utf-8",
            )
            jobs = root / "jobs.txt"
            jobs.write_text(f"{patent}\n", encoding="utf-8")

            candidates, skips, conflicts, listed, pending = candidates_from_job_list(
                jobs, root / "merged.csv"
            )

            self.assertEqual(candidates, [])
            self.assertEqual(skips, [marker.resolve()])
            self.assertEqual(conflicts, [])
            self.assertEqual(listed, 1)
            self.assertEqual(pending, 0)
            self.assertEqual(
                validate_skip_marker(marker)["audit_status"],
                "not_audited_for_false_negative",
            )

    def test_csv_and_skip_marker_are_reported_as_conflict(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            patent = root / "WO2024000002"
            output = patent / "output"
            output.mkdir(parents=True)
            (output / "WO2024000002_sequence_tm.csv").write_text("header\n", encoding="utf-8")
            (output / "WO2024000002_tm_skipped.json").write_text("{}\n", encoding="utf-8")
            jobs = root / "jobs.txt"
            jobs.write_text(f"{patent}\n", encoding="utf-8")

            candidates, skips, conflicts, listed, pending = candidates_from_job_list(
                jobs, root / "merged.csv"
            )

            self.assertEqual(candidates, [])
            self.assertEqual(skips, [])
            self.assertEqual(len(conflicts), 1)
            self.assertEqual(listed, 1)
            self.assertEqual(pending, 0)

    def test_relaxed_audit_flags_broad_tm_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            patent = root / "WO2024000003"
            output = patent / "output"
            output.mkdir(parents=True)
            patent_json = patent / "WO2024000003.json"
            patent_json.write_text('{"text": "OCR table has $T_{{m}}$ data"}', encoding="utf-8")
            marker = output / "WO2024000003_tm_skipped.json"
            marker.write_text("{}\n", encoding="utf-8")
            item = {"patent_id": patent.name, "path": str(marker)}

            audited = audit_skipped_patent(
                item,
                lambda *_args, **_kwargs: {"tm_hit_count": 0},
            )

            self.assertEqual(audited["audit_status"], "needs_agent_review_broad_tm_evidence")
            self.assertTrue(audited["broad_lexical_hit"])


if __name__ == "__main__":
    unittest.main()

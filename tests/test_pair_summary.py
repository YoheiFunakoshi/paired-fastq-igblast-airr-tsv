from __future__ import annotations

from pathlib import Path
import csv
import shutil
import unittest
import uuid

from airr_igblast_paired.pair_summary import (
    default_derived_tsv_paths,
    pair_id_and_read_label,
    split_and_integrate_airr_tsv,
)


class PairSummaryTests(unittest.TestCase):
    def test_default_derived_paths_keep_sample_name(self) -> None:
        paths = default_derived_tsv_paths(Path("results") / "sample.airr.tsv")

        self.assertEqual(paths.r1_tsv, Path("results") / "sample.R1.airr.tsv")
        self.assertEqual(paths.r2_tsv, Path("results") / "sample.R2.airr.tsv")
        self.assertEqual(paths.integrated_tsv, Path("results") / "sample.integrated.tsv")

    def test_pair_id_and_read_label_uses_last_pipe_component(self) -> None:
        self.assertEqual(pair_id_and_read_label("read-1|R1"), ("read-1", "R1"))
        self.assertEqual(pair_id_and_read_label("read-1|extra|R2"), ("read-1|extra", "R2"))
        self.assertEqual(pair_id_and_read_label("read-1"), ("read-1", None))

    def test_split_and_integrate_keeps_conflicts_with_simple_rule(self) -> None:
        root = Path(f"test_pair_summary_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            input_tsv = root / "sample.airr.tsv"
            input_tsv.write_text(
                "\n".join(
                    [
                        "sequence_id\tv_call\td_call\tj_call\tjunction\tjunction_aa\tproductive",
                        "readA|R1\tIGHV1\tIGHD1\tIGHJ4\tAAATTT\tCAR\tT",
                        "readA|R2\tIGHV1\tIGHD1\tIGHJ4\tAAATTT\tCAR\tT",
                        "readB|R1\tIGHV3\tIGHD2\tIGHJ5\tCCCC\tCA\tT",
                        "readB|R2\tIGHV4\tIGHD2\tIGHJ6\tCCCCGG\tCARG\tT",
                        "readC|R2\tIGHV7\t\tIGHJ2\tGGGG\tCSS\tF",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            paths, stats = split_and_integrate_airr_tsv(input_tsv)

            self.assertEqual(stats.total_pairs, 3)
            self.assertEqual(stats.r1_rows, 2)
            self.assertEqual(stats.r2_rows, 3)
            self.assertEqual(stats.junction_aa_conflicts, 1)
            self.assertIn("readA|R1", paths.r1_tsv.read_text(encoding="utf-8"))
            self.assertIn("readA|R2", paths.r2_tsv.read_text(encoding="utf-8"))

            with paths.integrated_tsv.open("rt", encoding="utf-8", newline="") as handle:
                rows = {row["pair_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

            self.assertEqual(rows["readA"]["junction_aa_status"], "match")
            self.assertEqual(rows["readB"]["junction_aa_status"], "conflict")
            self.assertEqual(rows["readB"]["final_junction_aa"], "CARG")
            self.assertEqual(rows["readB"]["preferred_read"], "R2")
            self.assertEqual(rows["readB"]["final_v_call"], "IGHV4")
            self.assertEqual(rows["readC"]["junction_aa_status"], "r2_only")
            self.assertEqual(rows["readC"]["usable_for_qasas"], "true")
        finally:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

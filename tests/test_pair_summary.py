from __future__ import annotations

from pathlib import Path
import csv
import shutil
import unittest
import uuid
import zipfile

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
        self.assertEqual(paths.counts_tsv, Path("results") / "sample.integrated_counts.tsv")
        self.assertEqual(paths.counts_xlsx, Path("results") / "sample.integrated_counts.xlsx")

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
                        "readA|R1\tIGHV1-1*01,IGHV1-2*01\tIGHD1\tIGHJ4*01\tAAATTT\tCARYW\tT",
                        "readA|R2\tIGHV1-2*02,IGHV1-1*02\tIGHD1\tIGHJ4*02\tAAATTT\tCARYW\tT",
                        "readA_dup|R1\tIGHV1-1*01,IGHV1-2*01\tIGHD1\tIGHJ4*01\tAAATTT\tCARYW\tT",
                        "readA_dup|R2\tIGHV1-2*02,IGHV1-1*02\tIGHD1\tIGHJ4*02\tAAATTT\tCARYW\tT",
                        "readB|R1\tIGHV3\tIGHD2\tIGHJ5\tCCCC\tCA\tT",
                        "readB|R2\tIGHV4\tIGHD2\tIGHJ6\tCCCCGG\tCARGW\tT",
                        "readC|R2\tIGHV7\t\tIGHJ2\tGGGG\tCSS\tF",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            paths, stats = split_and_integrate_airr_tsv(input_tsv)

            self.assertEqual(stats.total_pairs, 4)
            self.assertEqual(stats.r1_rows, 3)
            self.assertEqual(stats.r2_rows, 4)
            self.assertEqual(stats.junction_aa_conflicts, 1)
            self.assertEqual(stats.included_in_counts, 3)
            self.assertEqual(stats.unique_final_clonotypes, 2)
            self.assertIn("readA|R1", paths.r1_tsv.read_text(encoding="utf-8"))
            self.assertIn("readA|R2", paths.r2_tsv.read_text(encoding="utf-8"))

            with paths.integrated_tsv.open("rt", encoding="utf-8", newline="") as handle:
                rows = {row["pair_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

            self.assertEqual(rows["readA"]["junction_aa_status"], "match")
            self.assertEqual(rows["readB"]["junction_aa_status"], "conflict")
            self.assertEqual(rows["readB"]["final_junction_aa"], "CARGW")
            self.assertEqual(rows["readB"]["preferred_read"], "R2")
            self.assertEqual(rows["readB"]["final_v_call"], "IGHV4")
            self.assertEqual(rows["readC"]["junction_aa_status"], "r2_only")
            self.assertEqual(rows["readC"]["include_in_counts"], "false")
            self.assertIn("not_productive", rows["readC"]["exclude_reason"])

            with paths.counts_tsv.open("rt", encoding="utf-8", newline="") as handle:
                count_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(count_rows[0]["unique_v_gene_set"], "IGHV1-1,IGHV1-2")
            self.assertEqual(count_rows[0]["unique_j_gene_set"], "IGHJ4")
            self.assertEqual(count_rows[0]["final_junction_aa"], "CARYW")
            self.assertEqual(count_rows[0]["read_pair_count"], "2")
            self.assertEqual(count_rows[0]["match_count"], "2")
            self.assertEqual(count_rows[0]["productive_true_count"], "2")
            self.assertEqual(count_rows[0]["canonical_junction_aa_count"], "2")
            self.assertTrue(paths.counts_xlsx.exists())
            with zipfile.ZipFile(paths.counts_xlsx) as archive:
                self.assertIn("xl/worksheets/sheet1.xml", archive.namelist())
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_counts_use_exact_gene_sets_and_canonical_productive_filter(self) -> None:
        root = Path(f"test_pair_summary_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            input_tsv = root / "sample.airr.tsv"
            input_tsv.write_text(
                "\n".join(
                    [
                        "sequence_id\tv_call\td_call\tj_call\tjunction\tjunction_aa\tproductive",
                        "setAB_1|R2\tIGHV4-61*01,IGHV4-59*01\tIGHD1\tIGHJ4*01\tAAATTT\tCVQGFDYW\tT",
                        "setA_1|R2\tIGHV4-61*02\tIGHD1\tIGHJ4*02\tAAATTT\tCVQGFDYW\tT",
                        "setAB_2|R2\tIGHV4-59*02,IGHV4-61*03\tIGHD2\tIGHJ4*03\tAAATTT\tCVQGFDYW\tT",
                        "bad_start|R2\tIGHV4-61*01\tIGHD1\tIGHJ4*01\tAAATTT\tAVQGFDYW\tT",
                        "bad_productive|R2\tIGHV4-61*01\tIGHD1\tIGHJ4*01\tAAATTT\tCVQGFDYW\tF",
                        "missing_j|R2\tIGHV4-61*01\tIGHD1\t\tAAATTT\tCVQGFDYW\tT",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            paths, stats = split_and_integrate_airr_tsv(input_tsv)

            self.assertEqual(stats.total_pairs, 6)
            self.assertEqual(stats.included_in_counts, 3)
            self.assertEqual(stats.unique_final_clonotypes, 2)

            with paths.counts_tsv.open("rt", encoding="utf-8", newline="") as handle:
                count_rows = list(csv.DictReader(handle, delimiter="\t"))

            count_lookup = {
                (row["unique_v_gene_set"], row["unique_j_gene_set"], row["final_junction_aa"]): row
                for row in count_rows
            }
            self.assertEqual(
                count_lookup[("IGHV4-59,IGHV4-61", "IGHJ4", "CVQGFDYW")]["read_pair_count"],
                "2",
            )
            self.assertEqual(
                count_lookup[("IGHV4-61", "IGHJ4", "CVQGFDYW")]["read_pair_count"],
                "1",
            )

            with paths.integrated_tsv.open("rt", encoding="utf-8", newline="") as handle:
                integrated = {row["pair_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

            self.assertEqual(integrated["bad_start"]["include_in_counts"], "false")
            self.assertIn("junction_aa_not_c_start", integrated["bad_start"]["exclude_reason"])
            self.assertEqual(integrated["bad_productive"]["include_in_counts"], "false")
            self.assertIn("not_productive", integrated["bad_productive"]["exclude_reason"])
            self.assertEqual(integrated["missing_j"]["include_in_counts"], "false")
            self.assertIn("missing_j_call", integrated["missing_j"]["exclude_reason"])
        finally:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

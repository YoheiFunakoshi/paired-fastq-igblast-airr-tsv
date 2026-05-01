from __future__ import annotations

from pathlib import Path
import unittest

from airr_igblast_paired.naming import (
    RESULTS_FOLDER_NAME,
    default_output_tsv_path,
    fastq_stem,
    sample_name_from_fastqs,
)


class NamingTests(unittest.TestCase):
    def test_fastq_stem_removes_gzip_fastq_suffix(self) -> None:
        self.assertEqual(fastq_stem("sample_R1.fastq.gz"), "sample_R1")

    def test_sample_name_from_illumina_pair(self) -> None:
        sample = sample_name_from_fastqs(
            "KKF103hG_S57_L001_R1_001.fastq",
            "KKF103hG_S57_L001_R2_001.fastq",
        )

        self.assertEqual(sample, "KKF103hG_S57_L001")

    def test_sample_name_from_simple_pair(self) -> None:
        sample = sample_name_from_fastqs("patient-A_R1.fq.gz", "patient-A_R2.fq.gz")

        self.assertEqual(sample, "patient-A")

    def test_default_output_uses_results_folder_and_sample_name(self) -> None:
        output = default_output_tsv_path(
            "KKF103hG_S57_L001_R1_001.fastq",
            "KKF103hG_S57_L001_R2_001.fastq",
            Path("work"),
        )

        self.assertEqual(
            output,
            Path("work") / RESULTS_FOLDER_NAME / "KKF103hG_S57_L001.airr.tsv",
        )


if __name__ == "__main__":
    unittest.main()

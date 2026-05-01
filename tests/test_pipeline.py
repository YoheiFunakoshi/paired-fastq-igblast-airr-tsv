from __future__ import annotations

from pathlib import Path
import shutil
import unittest
import uuid
from unittest.mock import patch

from airr_igblast_paired.igblast import IgBlastConfig
from airr_igblast_paired.pipeline import run_paired_igblast


class PipelineTests(unittest.TestCase):
    def test_work_dir_copies_final_tsv_and_query_fasta(self) -> None:
        root = Path(f"test_pipeline_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            r1 = root / "sample_R1.fastq"
            r2 = root / "sample_R2.fastq"
            out = root / "results" / "sample.airr.tsv"
            query = root / "results" / "sample.queries.fasta"
            work = root / "work"
            r1.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            def fake_run_igblast(query_fasta: Path, output_tsv: Path, _: IgBlastConfig) -> list[str]:
                self.assertTrue(str(query_fasta).startswith(str(work)))
                self.assertTrue(str(output_tsv).startswith(str(work)))
                output_tsv.write_text("sequence_id\tsequence\nread1|R1\tAACCGG\n", encoding="utf-8")
                return ["igblastn", "-query", str(query_fasta), "-out", str(output_tsv)]

            with patch("airr_igblast_paired.pipeline.run_igblast", side_effect=fake_run_igblast):
                result = run_paired_igblast(
                    r1_path=r1,
                    r2_path=r2,
                    output_tsv=out,
                    query_fasta=query,
                    igblast_config=IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                    work_dir=work,
                )

            self.assertEqual(result.output_tsv, out)
            self.assertEqual(result.query_fasta, query)
            self.assertTrue(out.exists())
            self.assertTrue(query.exists())
            self.assertIn(">read1|R2\nAACCTT\n", query.read_text(encoding="utf-8"))
        finally:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from pathlib import Path
import unittest

from airr_igblast_paired.fastq import FastqRecord
from airr_igblast_paired.prepare import (
    ReadTransform,
    prepare_paired_fastq_to_fasta,
    reverse_complement,
    transform_sequence,
)


def record(read_id: str, sequence: str, quality: str | None = None) -> FastqRecord:
    if quality is None:
        quality = "I" * len(sequence)
    return FastqRecord(read_id=read_id, header=f"@{read_id}", sequence=sequence, quality=quality)


class PrepareTests(unittest.TestCase):
    def test_reverse_complement(self) -> None:
        self.assertEqual(reverse_complement("ACGTNry"), "ryNACGT")

    def test_transform_trims_before_reverse_complement(self) -> None:
        transformed = transform_sequence(
            record("read1", "AACCGGTT"),
            ReadTransform("reverse-complement", trim_left=2, trim_right=2),
        )

        self.assertEqual(transformed, "CCGG")

    def test_prepare_writes_r1_and_r2_as_separate_queries(self) -> None:
        r1_path = Path("test_prepare_R1.fastq")
        r2_path = Path("test_prepare_R2.fastq")
        out_path = Path("test_prepare.fasta")
        try:
            r1_path.write_text("@read1/1\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
            r2_path.write_text("@read1/2\nAAGGTT\n+\nIIIIII\n", encoding="utf-8")

            stats = prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path)
            output = out_path.read_text(encoding="utf-8")
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

        self.assertEqual(stats.total_pairs, 1)
        self.assertEqual(stats.records_written, 2)
        self.assertIn(">read1|R1\nAACCGG\n", output)
        self.assertIn(">read1|R2\nAACCTT\n", output)

    def test_prepare_can_filter_by_read_length(self) -> None:
        r1_path = Path("test_filter_R1.fastq")
        r2_path = Path("test_filter_R2.fastq")
        out_path = Path("test_filter.fasta")
        try:
            r1_path.write_text("@read1/1\nAAAA\n+\nIIII\n", encoding="utf-8")
            r2_path.write_text("@read1/2\nCCCC\n+\nIIII\n", encoding="utf-8")

            stats = prepare_paired_fastq_to_fasta(r1_path, r2_path, out_path, min_length=5)
            output = out_path.read_text(encoding="utf-8")
        finally:
            r1_path.unlink(missing_ok=True)
            r2_path.unlink(missing_ok=True)
            out_path.unlink(missing_ok=True)

        self.assertEqual(stats.records_written, 0)
        self.assertEqual(stats.skipped_too_short, 2)
        self.assertEqual(output, "")


if __name__ == "__main__":
    unittest.main()

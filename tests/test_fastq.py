from __future__ import annotations

from pathlib import Path
import unittest

from airr_igblast_paired.fastq import normalize_read_id, read_fastq


class FastqTests(unittest.TestCase):
    def test_normalize_read_id(self) -> None:
        self.assertEqual(normalize_read_id("@sample/1"), "sample")
        self.assertEqual(normalize_read_id("@instrument:run 2:N:0:ACGT"), "instrument:run")

    def test_read_fastq(self) -> None:
        path = Path("test_read_fastq.fastq")
        try:
            path.write_text("@read1/1\nACGT\n+\nIIII\n", encoding="utf-8")
            records = list(read_fastq(path))
        finally:
            path.unlink(missing_ok=True)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].read_id, "read1")
        self.assertEqual(records[0].sequence, "ACGT")
        self.assertEqual(records[0].quality, "IIII")


if __name__ == "__main__":
    unittest.main()

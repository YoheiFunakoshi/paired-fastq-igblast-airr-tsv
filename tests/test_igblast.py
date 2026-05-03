from __future__ import annotations

from pathlib import Path
import shutil
import unittest
import uuid
from unittest.mock import patch

from airr_igblast_paired.igblast import IgBlastConfig, build_igblast_command, run_igblast_batched


class IgBlastTests(unittest.TestCase):
    def test_default_threads_is_four(self) -> None:
        command = build_igblast_command(
            Path("queries.fasta"),
            Path("out.tsv"),
            IgBlastConfig(
                igblastn="igblastn",
                germline_db_v="human_gl_V",
                germline_db_j="human_gl_J",
            ),
        )

        self.assertEqual(command[command.index("-num_threads") + 1], "4")

    def test_build_igblast_command_uses_airr_outfmt_19(self) -> None:
        command = build_igblast_command(
            Path("queries.fasta"),
            Path("out.tsv"),
            IgBlastConfig(
                igblastn="igblastn",
                germline_db_v="human_gl_V",
                germline_db_d="human_gl_D",
                germline_db_j="human_gl_J",
                auxiliary_data="human_gl.aux",
                num_threads=4,
            ),
        )

        self.assertIn("-outfmt", command)
        self.assertEqual(command[command.index("-outfmt") + 1], "19")
        self.assertIn("-germline_db_D", command)
        self.assertIn("-auxiliary_data", command)
        self.assertEqual(command[command.index("-num_threads") + 1], "4")

    def test_run_igblast_batched_appends_one_header(self) -> None:
        root = Path(f"test_igblast_tmp_{uuid.uuid4().hex[:8]}")
        shutil.rmtree(root, ignore_errors=True)
        try:
            root.mkdir()
            query = root / "queries.fasta"
            output = root / "out.airr.tsv"
            query.write_text(
                ">read1|R1\nAAA\n"
                ">read1|R2\nTTT\n"
                ">read2|R1\nCCC\n",
                encoding="utf-8",
            )

            def fake_run_igblast(batch_query: Path, batch_output: Path, _: IgBlastConfig) -> list[str]:
                names = [
                    line[1:].strip()
                    for line in batch_query.read_text(encoding="utf-8").splitlines()
                    if line.startswith(">")
                ]
                batch_output.write_text(
                    "sequence_id\tproductive\n"
                    + "".join(f"{name}\tT\n" for name in names),
                    encoding="utf-8",
                )
                return ["igblastn", "-query", str(batch_query), "-out", str(batch_output)]

            with patch("airr_igblast_paired.igblast.run_igblast", side_effect=fake_run_igblast):
                command = run_igblast_batched(
                    query,
                    output,
                    IgBlastConfig(germline_db_v="v", germline_db_j="j"),
                    batch_size=2,
                )

            lines = output.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines.count("sequence_id\tproductive"), 1)
            self.assertEqual(lines[1:], ["read1|R1\tT", "read1|R2\tT", "read2|R1\tT"])
            self.assertIn("# batches", command)
            self.assertIn("2", command)
        finally:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

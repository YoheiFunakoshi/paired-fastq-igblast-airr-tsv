from __future__ import annotations

from pathlib import Path
import unittest

from airr_igblast_paired.igblast import IgBlastConfig, build_igblast_command


class IgBlastTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()

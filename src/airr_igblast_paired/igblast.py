from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import subprocess


@dataclass(frozen=True)
class IgBlastConfig:
    germline_db_v: str
    germline_db_j: str
    germline_db_d: str | None = None
    igblastn: str = "igblastn"
    organism: str = "human"
    domain_system: str = "imgt"
    ig_seqtype: str = "Ig"
    auxiliary_data: str | None = None
    num_threads: int = 1
    extra_args: list[str] = field(default_factory=list)


def build_igblast_command(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
) -> list[str]:
    command = [
        config.igblastn,
        "-query",
        str(query_fasta),
        "-out",
        str(output_tsv),
        "-outfmt",
        "19",
        "-germline_db_V",
        config.germline_db_v,
        "-germline_db_J",
        config.germline_db_j,
        "-organism",
        config.organism,
        "-domain_system",
        config.domain_system,
        "-ig_seqtype",
        config.ig_seqtype,
        "-num_threads",
        str(config.num_threads),
    ]

    if config.germline_db_d:
        command.extend(["-germline_db_D", config.germline_db_d])
    if config.auxiliary_data:
        command.extend(["-auxiliary_data", config.auxiliary_data])
    command.extend(config.extra_args)
    return command


def run_igblast(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
) -> list[str]:
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = build_igblast_command(query_fasta, output_tsv, config)
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        command_text = " ".join(command)
        message = result.stderr.strip() or result.stdout.strip() or "IgBLAST failed without output"
        raise RuntimeError(f"IgBLAST failed with exit code {result.returncode}\n{command_text}\n{message}")
    return command

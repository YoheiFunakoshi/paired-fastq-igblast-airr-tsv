from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import tempfile

from .igblast import IgBlastConfig, run_igblast
from .prepare import (
    PrepareStats,
    ReadSelection,
    ReadTransform,
    prepare_paired_fastq_to_fasta,
)


@dataclass(frozen=True)
class PipelineResult:
    stats: PrepareStats
    command: list[str]
    query_fasta: Path | None
    output_tsv: Path


def run_paired_igblast(
    *,
    r1_path: str | Path,
    r2_path: str | Path,
    output_tsv: str | Path,
    igblast_config: IgBlastConfig,
    query_fasta: str | Path | None = None,
    read_selection: ReadSelection = "both",
    r1_transform: ReadTransform = ReadTransform("forward"),
    r2_transform: ReadTransform = ReadTransform("reverse-complement"),
    query_name_template: str = "{read_id}|{read}",
    min_length: int = 0,
    max_n_rate: float = 1.0,
    strict_ids: bool = True,
) -> PipelineResult:
    if not str(r1_path).strip():
        raise ValueError("R1 FASTQ is required")
    if not str(r2_path).strip():
        raise ValueError("R2 FASTQ is required")
    if not str(output_tsv).strip():
        raise ValueError("Output TSV is required")
    if not igblast_config.germline_db_v.strip():
        raise ValueError("V DB prefix is required")
    if not igblast_config.germline_db_j.strip():
        raise ValueError("J DB prefix is required")

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    if query_fasta:
        query_fasta_path = Path(query_fasta)
        stats = prepare_paired_fastq_to_fasta(
            r1_path,
            r2_path,
            query_fasta_path,
            read_selection=read_selection,
            r1_transform=r1_transform,
            r2_transform=r2_transform,
            query_name_template=query_name_template,
            min_length=min_length,
            max_n_rate=max_n_rate,
            strict_ids=strict_ids,
        )
        command = run_igblast(query_fasta_path, output_tsv, igblast_config)
        return PipelineResult(stats, command, query_fasta_path, output_tsv)

    fd, temp_name = tempfile.mkstemp(
        prefix=f"{output_tsv.stem}.",
        suffix=".queries.fasta",
        dir=output_tsv.parent,
    )
    os.close(fd)
    temp_fasta = Path(temp_name)
    try:
        stats = prepare_paired_fastq_to_fasta(
            r1_path,
            r2_path,
            temp_fasta,
            read_selection=read_selection,
            r1_transform=r1_transform,
            r2_transform=r2_transform,
            query_name_template=query_name_template,
            min_length=min_length,
            max_n_rate=max_n_rate,
            strict_ids=strict_ids,
        )
        command = run_igblast(temp_fasta, output_tsv, igblast_config)
        return PipelineResult(stats, command, None, output_tsv)
    finally:
        temp_fasta.unlink(missing_ok=True)

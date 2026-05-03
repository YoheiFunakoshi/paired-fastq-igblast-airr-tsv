from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from contextlib import contextmanager
from collections.abc import Callable
import ctypes
import os
import shutil
import tempfile
import uuid

from .igblast import IgBlastConfig, run_igblast, run_igblast_batched
from .pair_summary import PairSummaryStats, split_and_integrate_airr_tsv
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
    r1_tsv: Path | None = None
    r2_tsv: Path | None = None
    integrated_tsv: Path | None = None
    counts_tsv: Path | None = None
    counts_xlsx: Path | None = None
    pair_summary_stats: PairSummaryStats | None = None


def _build_derived_outputs(output_tsv: Path) -> tuple[Path, Path, Path, Path, Path, PairSummaryStats]:
    derived_paths, pair_stats = split_and_integrate_airr_tsv(output_tsv)
    return (
        derived_paths.r1_tsv,
        derived_paths.r2_tsv,
        derived_paths.integrated_tsv,
        derived_paths.counts_tsv,
        derived_paths.counts_xlsx,
        pair_stats,
    )


def default_work_dir() -> Path:
    base = os.environ.get("LOCALAPPDATA")
    if base:
        return Path(base) / "PairedFastqIgblastAirrTsv" / "work"
    return Path(tempfile.gettempdir()) / "PairedFastqIgblastAirrTsv" / "work"


@contextmanager
def keep_windows_awake():
    if os.name != "nt":
        yield
        return

    es_continuous = 0x80000000
    es_system_required = 0x00000001
    es_awaymode_required = 0x00000040
    ctypes.windll.kernel32.SetThreadExecutionState(es_continuous | es_system_required | es_awaymode_required)
    try:
        yield
    finally:
        ctypes.windll.kernel32.SetThreadExecutionState(es_continuous)


def _run_igblast_maybe_batched(
    query_fasta: Path,
    output_tsv: Path,
    igblast_config: IgBlastConfig,
    igblast_batch_size: int | None,
    progress_callback: Callable[[str], None] | None = None,
) -> list[str]:
    if igblast_batch_size and igblast_batch_size > 0:
        return run_igblast_batched(
            query_fasta,
            output_tsv,
            igblast_config,
            batch_size=igblast_batch_size,
            progress_callback=progress_callback,
        )
    if progress_callback:
        progress_callback("Starting IgBLAST without batching...")
    return run_igblast(query_fasta, output_tsv, igblast_config)


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
    igblast_batch_size: int | None = None,
    progress_callback: Callable[[str], None] | None = None,
    work_dir: str | Path | None = None,
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
    final_query_fasta = Path(query_fasta) if query_fasta else None
    if final_query_fasta:
        final_query_fasta.parent.mkdir(parents=True, exist_ok=True)

    if work_dir:
        work_root = Path(work_dir)
        work_root.mkdir(parents=True, exist_ok=True)
        scratch_dir = work_root / f"{output_tsv.stem}.{uuid.uuid4().hex[:8]}"
        scratch_dir.mkdir()
        scratch_query = scratch_dir / (final_query_fasta.name if final_query_fasta else f"{output_tsv.stem}.queries.fasta")
        scratch_output = scratch_dir / output_tsv.name
        success = False
        try:
            with keep_windows_awake():
                if progress_callback:
                    progress_callback("Preparing IgBLAST query FASTA...")
                stats = prepare_paired_fastq_to_fasta(
                    r1_path,
                    r2_path,
                    scratch_query,
                    read_selection=read_selection,
                    r1_transform=r1_transform,
                    r2_transform=r2_transform,
                    query_name_template=query_name_template,
                    min_length=min_length,
                    max_n_rate=max_n_rate,
                    strict_ids=strict_ids,
                )
                command = _run_igblast_maybe_batched(
                    scratch_query,
                    scratch_output,
                    igblast_config,
                    igblast_batch_size,
                    progress_callback,
                )
                if progress_callback:
                    progress_callback("Copying final AIRR TSV to Results folder...")
                shutil.copy2(scratch_output, output_tsv)
                if progress_callback:
                    progress_callback("Creating R1/R2, integrated, and count outputs...")
                r1_tsv, r2_tsv, integrated_tsv, counts_tsv, counts_xlsx, pair_stats = _build_derived_outputs(output_tsv)
                if final_query_fasta:
                    shutil.copy2(scratch_query, final_query_fasta)
            success = True
            return PipelineResult(
                stats,
                command,
                final_query_fasta,
                output_tsv,
                r1_tsv,
                r2_tsv,
                integrated_tsv,
                counts_tsv,
                counts_xlsx,
                pair_stats,
            )
        finally:
            if success:
                shutil.rmtree(scratch_dir, ignore_errors=True)

    if final_query_fasta:
        with keep_windows_awake():
            if progress_callback:
                progress_callback("Preparing IgBLAST query FASTA...")
            stats = prepare_paired_fastq_to_fasta(
                r1_path,
                r2_path,
                final_query_fasta,
                read_selection=read_selection,
                r1_transform=r1_transform,
                r2_transform=r2_transform,
                query_name_template=query_name_template,
                min_length=min_length,
                max_n_rate=max_n_rate,
                strict_ids=strict_ids,
            )
            command = _run_igblast_maybe_batched(
                final_query_fasta,
                output_tsv,
                igblast_config,
                igblast_batch_size,
                progress_callback,
            )
            if progress_callback:
                progress_callback("Creating R1/R2, integrated, and count outputs...")
            r1_tsv, r2_tsv, integrated_tsv, counts_tsv, counts_xlsx, pair_stats = _build_derived_outputs(output_tsv)
        return PipelineResult(
            stats,
            command,
            final_query_fasta,
            output_tsv,
            r1_tsv,
            r2_tsv,
            integrated_tsv,
            counts_tsv,
            counts_xlsx,
            pair_stats,
        )

    fd, temp_name = tempfile.mkstemp(
        prefix=f"{output_tsv.stem}.",
        suffix=".queries.fasta",
        dir=output_tsv.parent,
    )
    os.close(fd)
    temp_fasta = Path(temp_name)
    try:
        with keep_windows_awake():
            if progress_callback:
                progress_callback("Preparing IgBLAST query FASTA...")
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
            command = _run_igblast_maybe_batched(
                temp_fasta,
                output_tsv,
                igblast_config,
                igblast_batch_size,
                progress_callback,
            )
            if progress_callback:
                progress_callback("Creating R1/R2, integrated, and count outputs...")
            r1_tsv, r2_tsv, integrated_tsv, counts_tsv, counts_xlsx, pair_stats = _build_derived_outputs(output_tsv)
        return PipelineResult(
            stats,
            command,
            None,
            output_tsv,
            r1_tsv,
            r2_tsv,
            integrated_tsv,
            counts_tsv,
            counts_xlsx,
            pair_stats,
        )
    finally:
        temp_fasta.unlink(missing_ok=True)

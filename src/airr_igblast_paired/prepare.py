from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from .fastq import FastqRecord, read_fastq


Orientation = Literal["forward", "reverse-complement"]
ReadSelection = Literal["both", "r1", "r2"]

_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


@dataclass(frozen=True)
class ReadTransform:
    orientation: Orientation = "forward"
    trim_left: int = 0
    trim_right: int = 0


@dataclass
class PrepareStats:
    total_pairs: int = 0
    records_written: int = 0
    r1_written: int = 0
    r2_written: int = 0
    skipped_too_short: int = 0
    skipped_n_rate: int = 0


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def transform_sequence(record: FastqRecord, transform: ReadTransform) -> str:
    if transform.trim_left < 0 or transform.trim_right < 0:
        raise ValueError("trim values must be 0 or greater")

    sequence = record.sequence.upper()
    if transform.trim_left:
        sequence = sequence[transform.trim_left :]
    if transform.trim_right:
        sequence = sequence[: -transform.trim_right]

    if transform.orientation == "forward":
        return sequence
    if transform.orientation == "reverse-complement":
        return reverse_complement(sequence)
    raise ValueError(f"unsupported orientation: {transform.orientation}")


def n_rate(sequence: str) -> float:
    if not sequence:
        return 1.0
    return sequence.count("N") / len(sequence)


def should_write_sequence(sequence: str, *, min_length: int, max_n_rate: float) -> tuple[bool, str | None]:
    if min_length < 0:
        raise ValueError("min_length must be 0 or greater")
    if not 0 <= max_n_rate <= 1:
        raise ValueError("max_n_rate must be between 0 and 1")
    if len(sequence) < min_length:
        return False, "too_short"
    if n_rate(sequence) > max_n_rate:
        return False, "n_rate"
    return True, None


def write_fasta_record(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def make_query_name(read_id: str, read_label: Literal["R1", "R2"], template: str) -> str:
    try:
        query_name = template.format(read_id=read_id, read=read_label)
    except KeyError as exc:
        raise ValueError("query name template may only use {read_id} and {read}") from exc
    if not query_name or any(char.isspace() for char in query_name):
        raise ValueError("query names must be non-empty and contain no whitespace")
    return query_name


def _write_one_read(
    *,
    fasta,
    stats: PrepareStats,
    read_id: str,
    read_label: Literal["R1", "R2"],
    record: FastqRecord,
    transform: ReadTransform,
    query_name_template: str,
    min_length: int,
    max_n_rate: float,
) -> None:
    sequence = transform_sequence(record, transform)
    should_write, reason = should_write_sequence(sequence, min_length=min_length, max_n_rate=max_n_rate)
    if not should_write:
        if reason == "too_short":
            stats.skipped_too_short += 1
        elif reason == "n_rate":
            stats.skipped_n_rate += 1
        return

    query_name = make_query_name(read_id, read_label, query_name_template)
    write_fasta_record(fasta, query_name, sequence)
    stats.records_written += 1
    if read_label == "R1":
        stats.r1_written += 1
    else:
        stats.r2_written += 1


def prepare_paired_fastq_to_fasta(
    r1_path: str | Path,
    r2_path: str | Path,
    fasta_path: str | Path,
    *,
    read_selection: ReadSelection = "both",
    r1_transform: ReadTransform = ReadTransform("forward"),
    r2_transform: ReadTransform = ReadTransform("reverse-complement"),
    query_name_template: str = "{read_id}|{read}",
    min_length: int = 0,
    max_n_rate: float = 1.0,
    strict_ids: bool = True,
) -> PrepareStats:
    fasta_path = Path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    stats = PrepareStats()

    if read_selection not in ("both", "r1", "r2"):
        raise ValueError("read_selection must be one of: both, r1, r2")

    r1_iter = read_fastq(r1_path)
    r2_iter = read_fastq(r2_path)

    with fasta_path.open("wt", encoding="utf-8", newline="\n") as fasta:
        while True:
            try:
                r1 = next(r1_iter)
            except StopIteration:
                try:
                    next(r2_iter)
                except StopIteration:
                    return stats
                raise ValueError("R2 has more records than R1")

            try:
                r2 = next(r2_iter)
            except StopIteration as exc:
                raise ValueError("R1 has more records than R2") from exc

            if strict_ids and r1.read_id != r2.read_id:
                raise ValueError(f"read ID mismatch: R1={r1.read_id!r}, R2={r2.read_id!r}")

            stats.total_pairs += 1
            if read_selection in ("both", "r1"):
                _write_one_read(
                    fasta=fasta,
                    stats=stats,
                    read_id=r1.read_id,
                    read_label="R1",
                    record=r1,
                    transform=r1_transform,
                    query_name_template=query_name_template,
                    min_length=min_length,
                    max_n_rate=max_n_rate,
                )
            if read_selection in ("both", "r2"):
                _write_one_read(
                    fasta=fasta,
                    stats=stats,
                    read_id=r2.read_id,
                    read_label="R2",
                    record=r2,
                    transform=r2_transform,
                    query_name_template=query_name_template,
                    min_length=min_length,
                    max_n_rate=max_n_rate,
                )

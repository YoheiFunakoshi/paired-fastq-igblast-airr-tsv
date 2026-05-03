from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import os
import shutil
import subprocess
import ctypes
from collections.abc import Callable, Iterator


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
    num_threads: int = 4
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
        "-organism",
        config.organism,
        "-domain_system",
        config.domain_system,
        "-ig_seqtype",
        config.ig_seqtype,
        "-num_threads",
        str(config.num_threads),
    ]

    if config.germline_db_v:
        command.extend(["-germline_db_V", config.germline_db_v])
    if config.germline_db_d:
        command.extend(["-germline_db_D", config.germline_db_d])
    if config.germline_db_j:
        command.extend(["-germline_db_J", config.germline_db_j])
    if config.auxiliary_data:
        command.extend(["-auxiliary_data", config.auxiliary_data])
    command.extend(config.extra_args)
    return command


def _windows_short_path(path: str | Path) -> str:
    text = str(path)
    if os.name != "nt" or not text:
        return text

    buffer = ctypes.create_unicode_buffer(32768)
    result = ctypes.windll.kernel32.GetShortPathNameW(str(text), buffer, len(buffer))
    if result:
        return buffer.value
    return text


def _db_prefix_to_windows_short_path(prefix: str) -> str:
    if os.name != "nt" or not prefix:
        return prefix

    path = Path(prefix)
    if path.exists():
        return _windows_short_path(path)

    db_suffixes = (
        ".ndb",
        ".nhr",
        ".nin",
        ".nog",
        ".nos",
        ".not",
        ".nsq",
        ".ntf",
        ".nto",
        ".phr",
        ".pin",
        ".pog",
        ".psd",
        ".psi",
        ".psq",
    )
    if any(Path(str(path) + suffix).exists() for suffix in db_suffixes):
        return str(Path(_windows_short_path(path.parent)) / path.name)
    return prefix


def _file_to_windows_short_path(path_text: str) -> str:
    if os.name != "nt" or not path_text:
        return path_text

    path = Path(path_text)
    if path.exists():
        return _windows_short_path(path)
    if path.parent.exists():
        return str(Path(_windows_short_path(path.parent)) / path.name)
    return path_text


def _normalize_command_for_windows(command: list[str]) -> list[str]:
    if os.name != "nt":
        return command

    normalized = list(command)
    path_flags = {
        "-query": "file",
        "-out": "file",
        "-germline_db_V": "db",
        "-germline_db_D": "db",
        "-germline_db_J": "db",
        "-auxiliary_data": "file",
    }
    for index, value in enumerate(normalized):
        if index == 0:
            resolved = shutil.which(value) or value
            normalized[index] = _windows_short_path(resolved) if Path(resolved).exists() else value
            continue
        previous = normalized[index - 1] if index > 0 else ""
        if previous in path_flags:
            if path_flags[previous] == "db":
                normalized[index] = _db_prefix_to_windows_short_path(value)
            else:
                normalized[index] = _file_to_windows_short_path(value)
    return normalized


def _command_value(command: list[str], flag: str) -> str | None:
    try:
        index = command.index(flag)
    except ValueError:
        return None
    if index + 1 >= len(command):
        return None
    return command[index + 1]


def _refdata_root_from_command(command: list[str]) -> Path | None:
    for flag in ("-germline_db_V", "-germline_db_D", "-germline_db_J"):
        value = _command_value(command, flag)
        if not value:
            continue
        prefix = Path(value)
        parent = prefix.parent
        if parent.name.lower() != "db":
            continue
        root = parent.parent
        if (root / "internal_data").exists():
            return root
    return None


def _igblast_runtime_context(command: list[str]) -> tuple[Path | None, dict[str, str]]:
    env = os.environ.copy()
    if os.name != "nt":
        return None, env

    refdata_root = _refdata_root_from_command(command)
    if refdata_root is not None:
        env["IGDATA"] = _windows_short_path(refdata_root)
        return None, env

    resolved = shutil.which(command[0]) or command[0]
    exe = Path(resolved)
    install_root = exe.parent.parent if exe.exists() else None
    internal_data = install_root / "internal_data" if install_root else None
    if internal_data and internal_data.exists():
        # IgBLAST looks up annotation files relative to IGDATA.  On Windows,
        # using the short cwd plus IGDATA='.' avoids BLAST database path issues
        # when the install path contains spaces such as "Program Files".
        env["IGDATA"] = "."
        return Path(_windows_short_path(internal_data)), env
    return None, env


def run_igblast(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
) -> list[str]:
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = _normalize_command_for_windows(build_igblast_command(query_fasta, output_tsv, config))
    cwd, env = _igblast_runtime_context(command)
    result = subprocess.run(command, capture_output=True, text=True, check=False, cwd=cwd, env=env)
    if result.returncode != 0:
        command_text = " ".join(command)
        message = result.stderr.strip() or result.stdout.strip() or "IgBLAST failed without output"
        raise RuntimeError(f"IgBLAST failed with exit code {result.returncode}\n{command_text}\n{message}")
    return command


def _read_fasta_records(path: str | Path) -> Iterator[list[str]]:
    record: list[str] = []
    with Path(path).open("rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith(">") and record:
                yield record
                record = []
            record.append(line)
    if record:
        yield record


def _write_fasta_batch(path: Path, records: list[list[str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.writelines(record)


def _append_airr_tsv_batch(final_tsv: Path, batch_tsv: Path, *, wrote_header: bool) -> bool:
    with batch_tsv.open("rt", encoding="utf-8", newline="") as source, final_tsv.open(
        "at" if wrote_header else "wt",
        encoding="utf-8",
        newline="",
    ) as target:
        for line_number, line in enumerate(source):
            if wrote_header and line_number == 0 and line.startswith("sequence_id\t"):
                continue
            target.write(line)
            if line_number == 0 and line.startswith("sequence_id\t"):
                wrote_header = True
    return wrote_header


def run_igblast_batched(
    query_fasta: str | Path,
    output_tsv: str | Path,
    config: IgBlastConfig,
    *,
    batch_size: int,
    progress_callback: Callable[[str], None] | None = None,
) -> list[str]:
    if batch_size <= 0:
        raise ValueError("batch_size must be greater than 0")

    query_fasta = Path(query_fasta)
    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.unlink(missing_ok=True)

    commands: list[list[str]] = []
    records: list[list[str]] = []
    batch_index = 0
    wrote_header = False

    def flush_batch() -> None:
        nonlocal batch_index, wrote_header, records
        if not records:
            return
        batch_index += 1
        batch_query = output_tsv.parent / f"{output_tsv.stem}.batch{batch_index:04d}.queries.fasta"
        batch_output = output_tsv.parent / f"{output_tsv.stem}.batch{batch_index:04d}.airr.tsv"
        query_count = len(records)
        try:
            if progress_callback:
                progress_callback(f"Starting IgBLAST batch {batch_index} ({query_count} queries)...")
            _write_fasta_batch(batch_query, records)
            commands.append(run_igblast(batch_query, batch_output, config))
            wrote_header = _append_airr_tsv_batch(output_tsv, batch_output, wrote_header=wrote_header)
            if progress_callback:
                progress_callback(f"Finished IgBLAST batch {batch_index}.")
        finally:
            batch_query.unlink(missing_ok=True)
            batch_output.unlink(missing_ok=True)
            records = []

    for record in _read_fasta_records(query_fasta):
        records.append(record)
        if len(records) >= batch_size:
            flush_batch()
    flush_batch()

    if not commands:
        output_tsv.write_text("", encoding="utf-8")
        return []

    first_command = list(commands[0])
    first_command.extend(["# batches", str(len(commands)), "# batch_size", str(batch_size)])
    return first_command

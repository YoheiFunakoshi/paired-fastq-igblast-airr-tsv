from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import os
import shutil
import subprocess
import ctypes


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
    if any(path.with_suffix(suffix).exists() for suffix in db_suffixes):
        return str(Path(_windows_short_path(path.parent)) / path.name)
    return prefix


def _normalize_command_for_windows(command: list[str]) -> list[str]:
    if os.name != "nt":
        return command

    normalized = list(command)
    path_flags = {
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
            elif Path(value).exists():
                normalized[index] = _windows_short_path(value)
    return normalized


def _igblast_runtime_context(igblastn: str) -> tuple[Path | None, dict[str, str]]:
    env = os.environ.copy()
    if os.name != "nt":
        return None, env

    resolved = shutil.which(igblastn) or igblastn
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
    cwd, env = _igblast_runtime_context(command[0])
    result = subprocess.run(command, capture_output=True, text=True, check=False, cwd=cwd, env=env)
    if result.returncode != 0:
        command_text = " ".join(command)
        message = result.stderr.strip() or result.stdout.strip() or "IgBLAST failed without output"
        raise RuntimeError(f"IgBLAST failed with exit code {result.returncode}\n{command_text}\n{message}")
    return command

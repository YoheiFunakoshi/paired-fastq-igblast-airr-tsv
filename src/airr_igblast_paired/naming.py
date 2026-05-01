from __future__ import annotations

from pathlib import Path
import os
import re


PROJECT_FOLDER_NAME = "Paired Fastq IgBLAST AIRR tsv"
RESULTS_FOLDER_NAME = "Results of Paired Fastq IgBLAST AIRR tsv"

_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_INVALID_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*]+')
_READ_MARKER_PATTERNS = (
    re.compile(r"(?i)([._-])R[12]([._-]\d{3})?$"),
    re.compile(r"(?i)([._-])read[12]$"),
    re.compile(r"(?i)([._-])[12]$"),
)


def default_data_folder() -> Path:
    return Path.home() / "Desktop" / PROJECT_FOLDER_NAME


def default_results_folder(data_folder: str | Path | None = None) -> Path:
    root = Path(data_folder) if data_folder else default_data_folder()
    return root / RESULTS_FOLDER_NAME


def fastq_stem(path: str | Path) -> str:
    name = Path(path).name
    lower_name = name.lower()
    for suffix in _FASTQ_SUFFIXES:
        if lower_name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def sample_name_from_fastqs(r1_path: str | Path, r2_path: str | Path | None = None) -> str:
    r1_stem = fastq_stem(r1_path)
    if not r2_path:
        return _clean_sample_name(_strip_read_marker(r1_stem))

    r2_stem = fastq_stem(r2_path)
    r1_sample = _strip_read_marker(r1_stem)
    r2_sample = _strip_read_marker(r2_stem)

    if r1_sample.lower() == r2_sample.lower():
        return _clean_sample_name(r1_sample)

    common = os.path.commonprefix([r1_stem, r2_stem]).rstrip("._- ")
    if common:
        return _clean_sample_name(_strip_read_marker(common))

    return _clean_sample_name(r1_sample)


def default_output_tsv_path(
    r1_path: str | Path,
    r2_path: str | Path | None = None,
    data_folder: str | Path | None = None,
) -> Path:
    sample_name = sample_name_from_fastqs(r1_path, r2_path)
    return default_results_folder(data_folder) / f"{sample_name}.airr.tsv"


def default_query_fasta_path(
    r1_path: str | Path,
    r2_path: str | Path | None = None,
    data_folder: str | Path | None = None,
) -> Path:
    sample_name = sample_name_from_fastqs(r1_path, r2_path)
    return default_results_folder(data_folder) / f"{sample_name}.queries.fasta"


def _strip_read_marker(stem: str) -> str:
    cleaned = stem
    for pattern in _READ_MARKER_PATTERNS:
        cleaned = pattern.sub("", cleaned)
        if cleaned != stem:
            break
    return cleaned.rstrip("._- ")


def _clean_sample_name(name: str) -> str:
    cleaned = _INVALID_FILENAME_CHARS.sub("_", name).strip(" ._")
    return cleaned or "sample"

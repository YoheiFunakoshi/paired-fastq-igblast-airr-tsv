from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO
import gzip


@dataclass(frozen=True)
class FastqRecord:
    read_id: str
    header: str
    sequence: str
    quality: str


def normalize_read_id(header: str) -> str:
    text = header.strip()
    if text.startswith("@"):
        text = text[1:]
    read_id = text.split()[0]
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]
    return read_id


def open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def read_fastq(path: str | Path) -> Iterator[FastqRecord]:
    path = Path(path)
    with open_text(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if header == "":
                return

            record_number += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError(f"{path}: incomplete FASTQ record at record {record_number}")

            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")

            if not header.startswith("@"):
                raise ValueError(f"{path}: FASTQ header does not start with @ at record {record_number}")
            if not plus.startswith("+"):
                raise ValueError(f"{path}: FASTQ plus line does not start with + at record {record_number}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: sequence and quality lengths differ at record {record_number}"
                )

            yield FastqRecord(
                read_id=normalize_read_id(header),
                header=header,
                sequence=sequence,
                quality=quality,
            )

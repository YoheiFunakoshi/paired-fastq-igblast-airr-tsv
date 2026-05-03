from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Literal


ReadLabel = Literal["R1", "R2"]


@dataclass(frozen=True)
class DerivedTsvPaths:
    r1_tsv: Path
    r2_tsv: Path
    integrated_tsv: Path
    counts_tsv: Path


@dataclass(frozen=True)
class PairSummaryStats:
    total_pairs: int
    r1_rows: int
    r2_rows: int
    junction_aa_conflicts: int
    unique_final_clonotypes: int


def default_derived_tsv_paths(output_tsv: str | Path) -> DerivedTsvPaths:
    output = Path(output_tsv)
    name = output.name
    lower_name = name.lower()
    if lower_name.endswith(".airr.tsv"):
        sample = name[: -len(".airr.tsv")]
    elif lower_name.endswith(".tsv"):
        sample = name[: -len(".tsv")]
    else:
        sample = output.stem

    return DerivedTsvPaths(
        r1_tsv=output.with_name(f"{sample}.R1.airr.tsv"),
        r2_tsv=output.with_name(f"{sample}.R2.airr.tsv"),
        integrated_tsv=output.with_name(f"{sample}.integrated.tsv"),
        counts_tsv=output.with_name(f"{sample}.integrated_counts.tsv"),
    )


def split_and_integrate_airr_tsv(
    input_tsv: str | Path,
    paths: DerivedTsvPaths | None = None,
) -> tuple[DerivedTsvPaths, PairSummaryStats]:
    input_path = Path(input_tsv)
    derived = paths or default_derived_tsv_paths(input_path)
    for path in (derived.r1_tsv, derived.r2_tsv, derived.integrated_tsv, derived.counts_tsv):
        path.parent.mkdir(parents=True, exist_ok=True)

    pairs: dict[str, dict[ReadLabel, dict[str, str]]] = {}
    r1_rows = 0
    r2_rows = 0

    with input_path.open("rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames:
            _write_empty_tsv(derived.r1_tsv)
            _write_empty_tsv(derived.r2_tsv)
            _write_integrated_tsv(derived.integrated_tsv, [])
            _write_counts_tsv(derived.counts_tsv, [])
            return derived, PairSummaryStats(0, 0, 0, 0, 0)

        with (
            derived.r1_tsv.open("wt", encoding="utf-8", newline="") as r1_handle,
            derived.r2_tsv.open("wt", encoding="utf-8", newline="") as r2_handle,
        ):
            r1_writer = csv.DictWriter(r1_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            r2_writer = csv.DictWriter(r2_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            r1_writer.writeheader()
            r2_writer.writeheader()

            for row in reader:
                pair_id, read_label = pair_id_and_read_label(row.get("sequence_id", ""))
                if not pair_id or read_label is None:
                    continue
                pairs.setdefault(pair_id, {})[read_label] = row
                if read_label == "R1":
                    r1_writer.writerow(row)
                    r1_rows += 1
                else:
                    r2_writer.writerow(row)
                    r2_rows += 1

    integrated_rows = [_integrated_row(pair_id, pair_rows) for pair_id, pair_rows in sorted(pairs.items())]
    _write_integrated_tsv(derived.integrated_tsv, integrated_rows)
    counts_rows = _counts_rows(integrated_rows)
    _write_counts_tsv(derived.counts_tsv, counts_rows)
    conflicts = sum(1 for row in integrated_rows if row["junction_aa_status"] == "conflict")
    return derived, PairSummaryStats(len(integrated_rows), r1_rows, r2_rows, conflicts, len(counts_rows))


def pair_id_and_read_label(sequence_id: str) -> tuple[str, ReadLabel | None]:
    text = sequence_id.strip()
    if not text or "|" not in text:
        return text, None
    pair_id, read_label = text.rsplit("|", 1)
    if read_label in ("R1", "R2"):
        return pair_id, read_label
    return text, None


def _write_empty_tsv(path: Path) -> None:
    path.write_text("", encoding="utf-8", newline="")


def _is_value(value: str | None) -> bool:
    if value is None:
        return False
    text = value.strip()
    return bool(text) and text.lower() not in {"na", "n/a", "none", "null"}


def _get(row: dict[str, str] | None, field: str) -> str:
    if row is None:
        return ""
    return row.get(field, "") or ""


def _choose_junction_aa(r1: str, r2: str) -> tuple[str, str, str, str]:
    r1_has = _is_value(r1)
    r2_has = _is_value(r2)
    if r1_has and r2_has:
        if r1 == r2:
            return r1, "match", "both", "same_junction_aa"
        if len(r1) > len(r2):
            return r1, "conflict", "R1", "conflict_longer_r1"
        if len(r2) > len(r1):
            return r2, "conflict", "R2", "conflict_longer_r2"
        return r2, "conflict", "R2", "conflict_same_length_r2_priority"
    if r1_has:
        return r1, "r1_only", "R1", "only_r1_has_junction_aa"
    if r2_has:
        return r2, "r2_only", "R2", "only_r2_has_junction_aa"
    return "", "none", "", "no_junction_aa"


def _choose_prefer_r2(r1: str, r2: str) -> tuple[str, str, str]:
    r1_has = _is_value(r1)
    r2_has = _is_value(r2)
    if r1_has and r2_has and r1 == r2:
        return r1, "both", "match"
    if r2_has:
        return r2, "R2", "r2_priority"
    if r1_has:
        return r1, "R1", "r1_only"
    return "", "", "none"


def _choose_by_preferred_read(
    r1: str,
    r2: str,
    preferred_read: str,
    *,
    fallback_r2: bool = False,
) -> tuple[str, str, str]:
    r1_has = _is_value(r1)
    r2_has = _is_value(r2)
    if r1_has and r2_has and r1 == r2:
        return r1, "both", "match"
    if preferred_read == "R1" and r1_has:
        return r1, "R1", "preferred_read"
    if preferred_read == "R2" and r2_has:
        return r2, "R2", "preferred_read"
    if fallback_r2 and r2_has:
        return r2, "R2", "fallback_r2"
    if r1_has:
        return r1, "R1", "fallback_r1"
    if r2_has:
        return r2, "R2", "fallback_r2"
    return "", "", "none"


def _integrated_row(pair_id: str, pair_rows: dict[ReadLabel, dict[str, str]]) -> dict[str, str]:
    r1 = pair_rows.get("R1")
    r2 = pair_rows.get("R2")
    r1_junction_aa = _get(r1, "junction_aa")
    r2_junction_aa = _get(r2, "junction_aa")
    final_junction_aa, junction_status, preferred_read, decision_reason = _choose_junction_aa(
        r1_junction_aa,
        r2_junction_aa,
    )

    final_v, v_source, v_reason = _choose_prefer_r2(_get(r1, "v_call"), _get(r2, "v_call"))
    final_j, j_source, j_reason = _choose_by_preferred_read(
        _get(r1, "j_call"),
        _get(r2, "j_call"),
        preferred_read,
        fallback_r2=True,
    )
    final_d, d_source, d_reason = _choose_by_preferred_read(
        _get(r1, "d_call"),
        _get(r2, "d_call"),
        preferred_read,
        fallback_r2=True,
    )
    final_productive, productive_source, productive_reason = _choose_by_preferred_read(
        _get(r1, "productive"),
        _get(r2, "productive"),
        preferred_read,
        fallback_r2=True,
    )
    final_junction, junction_source, junction_reason = _choose_by_preferred_read(
        _get(r1, "junction"),
        _get(r2, "junction"),
        preferred_read,
        fallback_r2=True,
    )

    usable_for_qasas = all(_is_value(value) for value in (final_v, final_j, final_junction_aa))

    return {
        "pair_id": pair_id,
        "r1_sequence_id": _get(r1, "sequence_id"),
        "r2_sequence_id": _get(r2, "sequence_id"),
        "final_junction_aa": final_junction_aa,
        "junction_aa_status": junction_status,
        "preferred_read": preferred_read,
        "junction_aa_decision_reason": decision_reason,
        "r1_junction_aa": r1_junction_aa,
        "r2_junction_aa": r2_junction_aa,
        "final_junction": final_junction,
        "junction_source": junction_source,
        "junction_decision_reason": junction_reason,
        "r1_junction": _get(r1, "junction"),
        "r2_junction": _get(r2, "junction"),
        "final_v_call": final_v,
        "v_call_source": v_source,
        "v_call_decision_reason": v_reason,
        "r1_v_call": _get(r1, "v_call"),
        "r2_v_call": _get(r2, "v_call"),
        "final_d_call": final_d,
        "d_call_source": d_source,
        "d_call_decision_reason": d_reason,
        "r1_d_call": _get(r1, "d_call"),
        "r2_d_call": _get(r2, "d_call"),
        "final_j_call": final_j,
        "j_call_source": j_source,
        "j_call_decision_reason": j_reason,
        "r1_j_call": _get(r1, "j_call"),
        "r2_j_call": _get(r2, "j_call"),
        "final_productive": final_productive,
        "productive_source": productive_source,
        "productive_decision_reason": productive_reason,
        "r1_productive": _get(r1, "productive"),
        "r2_productive": _get(r2, "productive"),
        "usable_for_qasas": "true" if usable_for_qasas else "false",
    }


INTEGRATED_FIELDNAMES = [
    "pair_id",
    "r1_sequence_id",
    "r2_sequence_id",
    "final_junction_aa",
    "junction_aa_status",
    "preferred_read",
    "junction_aa_decision_reason",
    "r1_junction_aa",
    "r2_junction_aa",
    "final_junction",
    "junction_source",
    "junction_decision_reason",
    "r1_junction",
    "r2_junction",
    "final_v_call",
    "v_call_source",
    "v_call_decision_reason",
    "r1_v_call",
    "r2_v_call",
    "final_d_call",
    "d_call_source",
    "d_call_decision_reason",
    "r1_d_call",
    "r2_d_call",
    "final_j_call",
    "j_call_source",
    "j_call_decision_reason",
    "r1_j_call",
    "r2_j_call",
    "final_productive",
    "productive_source",
    "productive_decision_reason",
    "r1_productive",
    "r2_productive",
    "usable_for_qasas",
]


def _write_integrated_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INTEGRATED_FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


COUNTS_FIELDNAMES = [
    "final_v_call",
    "final_d_call",
    "final_j_call",
    "final_junction_aa",
    "read_pair_count",
    "match_count",
    "conflict_count",
    "r1_only_count",
    "r2_only_count",
    "none_count",
    "productive_true_count",
    "productive_false_count",
    "usable_for_qasas_count",
]


def _counts_rows(integrated_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts: dict[tuple[str, str, str], dict[str, int | str]] = {}
    for row in integrated_rows:
        key = (
            row.get("final_v_call", ""),
            row.get("final_d_call", ""),
            row.get("final_j_call", ""),
            row.get("final_junction_aa", ""),
        )
        bucket = counts.setdefault(
            key,
            {
                "final_v_call": key[0],
                "final_d_call": key[1],
                "final_j_call": key[2],
                "final_junction_aa": key[3],
                "read_pair_count": 0,
                "match_count": 0,
                "conflict_count": 0,
                "r1_only_count": 0,
                "r2_only_count": 0,
                "none_count": 0,
                "productive_true_count": 0,
                "productive_false_count": 0,
                "usable_for_qasas_count": 0,
            },
        )
        bucket["read_pair_count"] = int(bucket["read_pair_count"]) + 1
        status_key = row.get("junction_aa_status", "none")
        if status_key not in {"match", "conflict", "r1_only", "r2_only", "none"}:
            status_key = "none"
        bucket[f"{status_key}_count"] = int(bucket[f"{status_key}_count"]) + 1

        productive = row.get("final_productive", "").strip().lower()
        if productive in {"t", "true", "yes", "1"}:
            bucket["productive_true_count"] = int(bucket["productive_true_count"]) + 1
        elif productive in {"f", "false", "no", "0"}:
            bucket["productive_false_count"] = int(bucket["productive_false_count"]) + 1

        if row.get("usable_for_qasas", "").strip().lower() == "true":
            bucket["usable_for_qasas_count"] = int(bucket["usable_for_qasas_count"]) + 1

    sorted_rows = sorted(
        counts.values(),
        key=lambda item: (
            -int(item["read_pair_count"]),
            str(item["final_v_call"]),
            str(item["final_d_call"]),
            str(item["final_j_call"]),
            str(item["final_junction_aa"]),
        ),
    )
    return [{field: str(row[field]) for field in COUNTS_FIELDNAMES} for row in sorted_rows]


def _write_counts_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COUNTS_FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

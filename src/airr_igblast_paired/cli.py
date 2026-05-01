from __future__ import annotations

import argparse
import sys

from . import __version__
from .igblast import IgBlastConfig
from .pipeline import run_paired_igblast
from .prepare import PrepareStats, ReadTransform, prepare_paired_fastq_to_fasta


def _add_fastq_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--r1", required=True, help="R1 FASTQ path. .gz is supported.")
    parser.add_argument("--r2", required=True, help="R2 FASTQ path. .gz is supported.")
    parser.add_argument(
        "--read-selection",
        choices=("both", "r1", "r2"),
        default="both",
        help="Which reads to send to IgBLAST. R1/R2 are never merged.",
    )
    parser.add_argument(
        "--r1-orientation",
        choices=("forward", "reverse-complement"),
        default="forward",
        help="Orientation used for R1 query sequences.",
    )
    parser.add_argument(
        "--r2-orientation",
        choices=("forward", "reverse-complement"),
        default="reverse-complement",
        help="Orientation used for R2 query sequences.",
    )
    parser.add_argument("--trim-left-r1", type=int, default=0, help="Bases trimmed from the left side of R1.")
    parser.add_argument("--trim-right-r1", type=int, default=0, help="Bases trimmed from the right side of R1.")
    parser.add_argument("--trim-left-r2", type=int, default=0, help="Bases trimmed from the left side of R2.")
    parser.add_argument("--trim-right-r2", type=int, default=0, help="Bases trimmed from the right side of R2.")
    parser.add_argument("--min-length", type=int, default=0, help="Skip query sequences shorter than this length.")
    parser.add_argument(
        "--max-n-rate",
        type=float,
        default=1.0,
        help="Skip query sequences with an N fraction greater than this value.",
    )
    parser.add_argument(
        "--query-name-template",
        default="{read_id}|{read}",
        help="FASTA query name template. Allowed fields: {read_id}, {read}.",
    )
    parser.add_argument(
        "--allow-id-mismatch",
        action="store_true",
        help="Do not fail when normalized R1/R2 read IDs differ.",
    )


def _add_igblast_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--igblastn", default="igblastn", help="Path to igblastn executable.")
    parser.add_argument("--germline-db-v", required=True, help="IgBLAST V germline database prefix.")
    parser.add_argument("--germline-db-d", help="IgBLAST D germline database prefix.")
    parser.add_argument("--germline-db-j", required=True, help="IgBLAST J germline database prefix.")
    parser.add_argument("--organism", default="human", help="IgBLAST organism value.")
    parser.add_argument("--domain-system", default="imgt", help="IgBLAST domain system.")
    parser.add_argument("--ig-seqtype", default="Ig", help="IgBLAST sequence type.")
    parser.add_argument("--auxiliary-data", help="IgBLAST auxiliary data file, for example human_gl.aux.")
    parser.add_argument("--num-threads", type=int, default=1, help="Number of IgBLAST threads.")
    parser.add_argument(
        "--extra-igblast-arg",
        action="append",
        help="Additional igblastn argument. Repeat once per token.",
    )


def _r1_transform(args: argparse.Namespace) -> ReadTransform:
    return ReadTransform(args.r1_orientation, args.trim_left_r1, args.trim_right_r1)


def _r2_transform(args: argparse.Namespace) -> ReadTransform:
    return ReadTransform(args.r2_orientation, args.trim_left_r2, args.trim_right_r2)


def _igblast_config(args: argparse.Namespace) -> IgBlastConfig:
    return IgBlastConfig(
        igblastn=args.igblastn,
        germline_db_v=args.germline_db_v,
        germline_db_d=args.germline_db_d,
        germline_db_j=args.germline_db_j,
        organism=args.organism,
        domain_system=args.domain_system,
        ig_seqtype=args.ig_seqtype,
        auxiliary_data=args.auxiliary_data,
        num_threads=args.num_threads,
        extra_args=args.extra_igblast_arg or [],
    )


def _print_stats(stats: PrepareStats) -> None:
    print(
        "prepare stats: "
        f"total_pairs={stats.total_pairs}, "
        f"records_written={stats.records_written}, "
        f"r1_written={stats.r1_written}, "
        f"r2_written={stats.r2_written}, "
        f"skipped_too_short={stats.skipped_too_short}, "
        f"skipped_n_rate={stats.skipped_n_rate}",
        file=sys.stderr,
    )


def _prepare(args: argparse.Namespace) -> int:
    stats = prepare_paired_fastq_to_fasta(
        args.r1,
        args.r2,
        args.out_fasta,
        read_selection=args.read_selection,
        r1_transform=_r1_transform(args),
        r2_transform=_r2_transform(args),
        query_name_template=args.query_name_template,
        min_length=args.min_length,
        max_n_rate=args.max_n_rate,
        strict_ids=not args.allow_id_mismatch,
    )
    _print_stats(stats)
    print(f"wrote {args.out_fasta}", file=sys.stderr)
    return 0


def _run(args: argparse.Namespace) -> int:
    result = run_paired_igblast(
        r1_path=args.r1,
        r2_path=args.r2,
        output_tsv=args.out,
        igblast_config=_igblast_config(args),
        query_fasta=args.query_fasta,
        read_selection=args.read_selection,
        r1_transform=_r1_transform(args),
        r2_transform=_r2_transform(args),
        query_name_template=args.query_name_template,
        min_length=args.min_length,
        max_n_rate=args.max_n_rate,
        strict_ids=not args.allow_id_mismatch,
    )
    _print_stats(result.stats)
    print("ran: " + " ".join(result.command), file=sys.stderr)
    print(f"wrote {result.output_tsv}", file=sys.stderr)
    return 0


def _gui(_: argparse.Namespace) -> int:
    from .gui import main as gui_main

    gui_main()
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paired-fastq-igblast-airr-tsv",
        description="Create AIRR TSV output from paired FASTQ files using IgBLAST outfmt 19.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare",
        help="Build an IgBLAST query FASTA from R1/R2 without merging reads.",
    )
    _add_fastq_args(prepare_parser)
    prepare_parser.add_argument("--out-fasta", required=True, help="Output FASTA path for IgBLAST input.")
    prepare_parser.set_defaults(func=_prepare)

    run_parser = subparsers.add_parser("run", help="Prepare R1/R2 queries and run igblastn -outfmt 19.")
    _add_fastq_args(run_parser)
    _add_igblast_args(run_parser)
    run_parser.add_argument("--out", required=True, help="Output AIRR TSV path.")
    run_parser.add_argument("--query-fasta", help="Optional path to keep the IgBLAST query FASTA.")
    run_parser.set_defaults(func=_run)

    gui_parser = subparsers.add_parser("gui", help="Open the local GUI.")
    gui_parser.set_defaults(func=_gui)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

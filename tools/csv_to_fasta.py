#!/usr/bin/env python3
"""
Convert CSV file(s) with target_id and sequence columns into a FASTA file.

Usage examples:
  # Basic: one CSV to one FASTA
  python csv_to_fasta.py inputs.csv --output outputs.fasta

  # Multiple CSVs concatenated
  python csv_to_fasta.py a.csv b.csv c.csv -o merged.fasta

  # Custom column names and wrap sequences at 80 chars
  python csv_to_fasta.py data.csv -o out.fasta --id-column id --seq-column seq --wrap 80

  # Dedupe by target_id (keep first occurrence)
  python csv_to_fasta.py data.csv -o out.fasta --dedupe

Contract:
- Inputs: one or more CSV files; must include id and sequence columns (configurable)
- Output: FASTA where each record is two lines: ">id" then the raw sequence (optionally wrapped)
- Error modes: missing files/columns -> non-zero exit; empty or NA sequence rows are skipped with a warning
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import time
from typing import Iterable, List, Tuple

import pandas as pd


def read_sequences(
    csv_paths: Iterable[Path], id_col: str, seq_col: str, group_by: str = None
) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for p in csv_paths:
        if not p.exists():
            raise FileNotFoundError(f"CSV file not found: {p}")
        try:
            df = pd.read_csv(p, keep_default_na=False)
        except Exception as e:
            raise RuntimeError(f"Failed to read CSV: {p}: {e}") from e

        missing = [c for c in (id_col, seq_col) if c not in df.columns]
        if missing:
            raise KeyError(
                f"Missing required column(s) {missing} in {p}. Available: {list(df.columns)}"
            )

        # Keep only the needed columns, normalize types
        df = df[[id_col, seq_col, group_by]].copy()
        # Cast to string while preserving NaNs to filter later
        df[id_col] = df[id_col].astype(str)
        df[seq_col] = df[seq_col].astype(str)
        frames.append(df)

    if not frames:
        return pd.DataFrame(columns=[id_col, seq_col])
    return pd.concat(frames, axis=0, ignore_index=True)


def wrap_sequence(seq: str, line_len: int) -> str:
    if line_len <= 0:
        return seq
    return "\n".join(seq[i : i + line_len] for i in range(0, len(seq), line_len))


def write_fasta(
    rows: Iterable[Tuple[str, str]],
    out_path: Path,
    wrap: int = 0,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for tid, seq in rows:
            f.write(f">{tid}\n")
            if wrap and wrap > 0:
                f.write(wrap_sequence(seq, wrap) + "\n")
            else:
                f.write(seq + "\n")


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Convert CSV(s) with target_id and sequence columns into a FASTA file."
        )
    )
    parser.add_argument(
        "csv_files",
        nargs="+",
        type=Path,
        help="Path(s) to input CSV file(s)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output FASTA file path",
    )
    parser.add_argument(
        "--id-column",
        default="target_id",
        help="Name of the ID column (default: target_id)",
    )
    parser.add_argument(
        "--seq-column",
        default="sequence",
        help="Name of the sequence column (default: sequence)",
    )
    parser.add_argument(
        "--dedupe",
        action="store_true",
        help="If set, keep the first record per ID and drop duplicates",
    )
    parser.add_argument(
        "--wrap",
        type=int,
        default=0,
        help="Wrap sequence lines to this length (0 = no wrap)",
    )
    parser.add_argument(
        "--uppercase",
        action="store_true",
        help="Convert sequences to uppercase",
    )
    parser.add_argument(
        "--group-by",
        type=str,
        default=None,
        help="Optional column name to group by, selects longest then first occurrence per group",
    )

    args = parser.parse_args(argv)
    group_by = args.group_by
    try:
        df = read_sequences(
            args.csv_files, args.id_column, args.seq_column, group_by=group_by
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2

    # Clean up rows: drop NA/empty sequences or IDs
    def _norm(s: str) -> str:
        return (s or "").strip()

    df[args.id_column] = df[args.id_column].map(_norm)
    df[args.seq_column] = df[args.seq_column].map(_norm)
    before = len(df)
    df = df[(df[args.id_column] != "") & (df[args.seq_column] != "")]
    dropped = before - len(df)
    if dropped > 0:
        print(
            f"Warning: skipped {dropped} row(s) with empty id/sequence", file=sys.stderr
        )
    # If grouping, select longest sequence per group, then first occurrence in original order
    if group_by is not None:
        df["seq_length"] = df[args.seq_column].map(len)
        df = df.sort_values(by="seq_length", ascending=False, kind="stable")
        df = df.drop_duplicates(subset=[group_by], keep="first")
        # Reassign ids to the grouped key
        df[args.id_column] = df[group_by]
        df = df.drop(columns=["seq_length", group_by])

    if args.dedupe:
        df = df.drop_duplicates(subset=[args.id_column], keep="first")

    if args.uppercase:
        df[args.seq_column] = df[args.seq_column].str.upper()

    rows = list(zip(df[args.id_column].tolist(), df[args.seq_column].tolist()))
    try:
        write_fasta(rows, args.output, wrap=args.wrap)
    except Exception as e:
        print(f"Error writing FASTA: {e}", file=sys.stderr)
        return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

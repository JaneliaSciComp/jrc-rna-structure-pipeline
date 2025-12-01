#!/usr/bin/env python3
"""
Script to merge two CSV files based on specified columns.

Usage:
    python merge_csv_files.py <left_file> <right_file> <output_file> \
        --left-on <column> --right-on <column> [--how <merge_type>]

Example:
    python merge_csv_files.py file1.csv file2.csv merged.csv \
        --left-on id --right-on id --how inner
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def merge_csv_files(
    left_file: str,
    right_file: str,
    output_file: str,
    left_on: str = None,
    right_on: str = None,
    on: str = None,
    how: str = "inner",
    suffixes: tuple = ("_left", "_right"),
) -> None:
    """
    Merge two CSV files based on specified columns.

    Parameters:
    -----------
    left_file : str
        Path to the left CSV file
    right_file : str
        Path to the right CSV file
    output_file : str
        Path to save the merged CSV file
    left_on : str, optional
        Column name in left DataFrame to merge on
    right_on : str, optional
        Column name in right DataFrame to merge on
    on : str, optional
        Column name to merge on (if same in both DataFrames)
    how : str, default 'inner'
        Type of merge: 'inner', 'outer', 'left', 'right'
    suffixes : tuple, default ('_left', '_right')
        Suffixes to apply to overlapping column names
    """
    # Check if files exist
    if not Path(left_file).exists():
        raise FileNotFoundError(f"Left file not found: {left_file}")
    if not Path(right_file).exists():
        raise FileNotFoundError(f"Right file not found: {right_file}")

    # Read CSV files
    print(f"Reading {left_file}...")
    df_left = pd.read_csv(left_file, keep_default_na=False)
    print(f"  Shape: {df_left.shape}")
    print(f"  Columns: {list(df_left.columns)}")

    print(f"\nReading {right_file}...")
    df_right = pd.read_csv(right_file, keep_default_na=False)
    print(f"  Shape: {df_right.shape}")
    print(f"  Columns: {list(df_right.columns)}")

    # Perform merge
    print(
        f"\nMerging on {'column: ' + str(on) if on else 'left: ' + str(left_on) + ', right: ' + str(right_on)}..."
    )
    print(f"Merge type: {how}")

    if on:
        df_merged = pd.merge(df_left, df_right, on=on, how=how, suffixes=suffixes)
    else:
        df_merged = pd.merge(
            df_left,
            df_right,
            left_on=left_on,
            right_on=right_on,
            how=how,
            suffixes=suffixes,
        )

    print(f"\nMerged DataFrame shape: {df_merged.shape}")
    print(f"Columns: {list(df_merged.columns)}")

    # Save merged file
    print(f"\nSaving to {output_file}...")
    df_merged.to_csv(output_file, index=False)
    print("Done!")

    # Print summary statistics
    print("\n--- Merge Summary ---")
    print(f"Left file rows: {len(df_left)}")
    print(f"Right file rows: {len(df_right)}")
    print(f"Merged rows: {len(df_merged)}")

    if how == "inner":
        print(f"Rows matched: {len(df_merged)}")
    elif how == "left":
        print(f"Left rows retained: {len(df_merged)}")
    elif how == "right":
        print(f"Right rows retained: {len(df_merged)}")
    elif how == "outer":
        print(f"Total unique rows: {len(df_merged)}")

    # Identify missing rows
    merge_key_left = left_on if left_on else on
    merge_key_right = right_on if right_on else on

    if merge_key_left and merge_key_right:
        # Get unique keys from each file
        left_keys = set(df_left[merge_key_left].dropna().drop_duplicates())
        right_keys = set(df_right[merge_key_right].dropna().drop_duplicates())
        print("Unique keys in left file:", len(left_keys))

        print("Unique keys in right file:", len(right_keys))
        # List those keys that are not unique
        non_unique_left = df_left[merge_key_left][
            df_left[merge_key_left].duplicated()
        ].drop_duplicates()
        non_unique_right = df_right[merge_key_right][
            df_right[merge_key_right].duplicated()
        ].drop_duplicates()
        print("Non-unique keys in left file:", len(non_unique_left))
        print("Non-unique keys in right file:", len(non_unique_right))

        # Find missing keys
        missing_in_right = left_keys - right_keys
        missing_in_left = right_keys - left_keys

        if missing_in_right:
            print(
                f"\n--- Keys in left file but NOT in right file: {len(missing_in_right)} ---"
            )
            if len(missing_in_right) <= 20:
                for key in sorted(missing_in_right):
                    print(f"  {key}")
            else:
                for key in sorted(list(missing_in_right)[:10]):
                    print(f"  {key}")
                print(f"  ... and {len(missing_in_right) - 10} more")

        if missing_in_left:
            print(
                f"\n--- Keys in right file but NOT in left file: {len(missing_in_left)} ---"
            )
            if len(missing_in_left) <= 20:
                for key in sorted(missing_in_left):
                    print(f"  {key}")
            else:
                for key in sorted(list(missing_in_left)[:10]):
                    print(f"  {key}")
                print(f"  ... and {len(missing_in_left) - 10} more")
        if len(non_unique_left):
            print(f"\n--- Non-unique keys in left file: {len(non_unique_left)} ---")
            if len(non_unique_left) <= 20:
                for key in sorted(non_unique_left):
                    print(f"  {key}")
            else:
                for key in sorted(list(non_unique_left)[:10]):
                    print(f"  {key}")
                print(f"  ... and {len(non_unique_left) - 10} more")
        if len(non_unique_right):
            print(f"\n--- Non-unique keys in right file: {len(non_unique_right)} ---")
            if len(non_unique_right) <= 20:
                for key in sorted(non_unique_right):
                    print(f"  {key}")
            else:
                for key in sorted(list(non_unique_right)[:10]):
                    print(f"  {key}")
                print(f"  ... and {len(non_unique_right) - 10} more")


def main():
    parser = argparse.ArgumentParser(
        description="Merge two CSV files based on specified columns.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Merge types (--how):
  inner    : Use intersection of keys from both frames (default)
  outer    : Use union of keys from both frames
  left     : Use only keys from left frame
  right    : Use only keys from right frame

Examples:
  # Merge on the same column name in both files
  python merge_csv_files.py file1.csv file2.csv output.csv --on id
  
  # Merge on different column names
  python merge_csv_files.py file1.csv file2.csv output.csv \\
      --left-on pdb_id --right-on structure_id
  
  # Left join (keep all rows from left file)
  python merge_csv_files.py file1.csv file2.csv output.csv \\
      --on id --how left
  
  # Outer join with custom suffixes for duplicate columns
  python merge_csv_files.py file1.csv file2.csv output.csv \\
      --on id --how outer --suffix-left _x --suffix-right _y
        """,
    )

    parser.add_argument("left_file", help="Path to the left CSV file")
    parser.add_argument("right_file", help="Path to the right CSV file")
    parser.add_argument("output_file", help="Path to save the merged CSV file")

    # Merge column arguments
    merge_group = parser.add_argument_group("merge columns")
    merge_group.add_argument(
        "--on",
        help="Column names to merge on (if same in both files)",
        nargs="+",
    )
    merge_group.add_argument(
        "--left-on", dest="left_on", help="Column name in left file to merge on"
    )
    merge_group.add_argument(
        "--right-on", dest="right_on", help="Column name in right file to merge on"
    )

    # Merge options
    parser.add_argument(
        "--how",
        choices=["inner", "outer", "left", "right"],
        default="inner",
        help="Type of merge to perform (default: inner)",
    )
    parser.add_argument(
        "--suffix-left",
        dest="suffix_left",
        default="_left",
        help="Suffix for overlapping columns from left file (default: _left)",
    )
    parser.add_argument(
        "--suffix-right",
        dest="suffix_right",
        default="_right",
        help="Suffix for overlapping columns from right file (default: _right)",
    )

    args = parser.parse_args()

    # Validate merge column arguments
    if args.on and (args.left_on or args.right_on):
        parser.error("Cannot specify both --on and --left-on/--right-on")

    if not args.on and not (args.left_on and args.right_on):
        parser.error("Must specify either --on or both --left-on and --right-on")

    try:
        merge_csv_files(
            left_file=args.left_file,
            right_file=args.right_file,
            output_file=args.output_file,
            left_on=args.left_on,
            right_on=args.right_on,
            on=args.on,
            how=args.how,
            suffixes=(args.suffix_left, args.suffix_right),
        )
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

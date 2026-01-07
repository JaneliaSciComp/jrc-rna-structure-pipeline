#!/usr/bin/env python3
"""
Add groups from one file to another (e.g., from metadata to bioassembly metadata).
Combines all groups for chains in the bioassembly.

For example, if a bioassembly contains chains A, B, C with group_ids 1, 2, 3,
the combined group_id would be "1,2,3".

This script handles:
- group_id
- seq_group_id
- mmseqs_* columns (any column matching the pattern)
"""

import argparse
import pandas as pd
import re
from pathlib import Path


def merge_groupings(
    right_file: str,
    left_file: str,
    output_file: str,
    right_on: str = "pdb_id",
    left_on: str = "pdb_id",
    group_columns: list[str] = None,
    group_pattern: str = r"^(group_id|seq_group_id|mmseqs_.+)$",
) -> None:
    """
    Merge grouping columns from right file to left file.

    For each entry in the left file, finds all matching chains in the right file
    and combines their group values.

    Parameters:
    -----------
    right_file : str
        Path to the source CSV file with grouping information
    left_file : str
        Path to the target CSV file to add groupings to (e.g., bioassembly metadata)
    output_file : str
        Path to save the merged CSV file
    right_on : str
        Column name in right DataFrame to match on (default: 'pdb_id')
    left_on : str
        Column name in left DataFrame to match on (default: 'pdb_id')
    group_columns : list[str], optional
        Explicit list of group columns to merge. If None, uses group_pattern.
    group_pattern : str
        Regex pattern to match group columns (default matches group_id, seq_group_id, mmseqs_*)
    """
    # Check if files exist
    if not Path(right_file).exists():
        raise FileNotFoundError(f"Right file not found: {right_file}")
    if not Path(left_file).exists():
        raise FileNotFoundError(f"Left file not found: {left_file}")

    # Read CSV files
    print(f"Reading {right_file}...")
    df_right = pd.read_csv(right_file, keep_default_na=False)
    print(f"  Shape: {df_right.shape}")
    print(f"  Columns: {list(df_right.columns)}")

    print(f"\nReading {left_file}...")
    df_left = pd.read_csv(left_file, keep_default_na=False)
    print(f"  Shape: {df_left.shape}")
    print(f"  Columns: {list(df_left.columns)}")

    # Identify group columns to merge
    if group_columns is None:
        pattern = re.compile(group_pattern)
        group_columns = [col for col in df_right.columns if pattern.match(col)]

    if not group_columns:
        raise ValueError(f"No group columns found matching pattern: {group_pattern}")

    print(f"\nGroup columns to merge: {group_columns}")

    # Check that required columns exist
    if right_on not in df_right.columns:
        raise ValueError(f"Column '{right_on}' not found in right file")
    if left_on not in df_left.columns:
        raise ValueError(f"Column '{left_on}' not found in left file")

    # Verify all group columns exist in right file
    missing_cols = [col for col in group_columns if col not in df_right.columns]
    if missing_cols:
        raise ValueError(f"Group columns not found in right file: {missing_cols}")

    # Preprocessing: Explode concatenated strings to lists
    df_right_processed = df_right.copy()

    # Convert all group columns to lists
    for col in group_columns:
        df_right_processed[col] = df_right_processed[col].str.split(";")

    # Prepare data: Select only needed columns from right file
    right_cols = [right_on] + group_columns
    df_right_subset = df_right_processed[right_cols].copy()

    # Merge and aggregate
    print("\nMerging left to right...")

    # Step 1: Merge left to right
    merged_df = pd.merge(
        df_left,
        df_right_subset,
        left_on=left_on,
        right_on=right_on,
        how="left",
    )

    print(f"Merged DataFrame shape: {merged_df.shape}")

    # Step 2: Group by left_on
    print(f"Grouping by '{left_on}'...")

    # Step 3: Aggregate the groups by extending lists
    def aggregate_lists(series):
        """Aggregate list values by extending them into a single list."""
        all_values = []
        for val in series:
            if isinstance(val, list):
                all_values.extend(val)
            elif pd.notna(val) and val != "":
                all_values.append(str(val))

        # Remove duplicates while preserving some order (sorted)
        unique_values = sorted(set(v for v in all_values if v and v != "nan"))
        return ";".join(unique_values) if unique_values else ""

    # Get all columns from left file
    left_columns = list(df_left.columns)

    # Create aggregation dictionary
    agg_dict = {}
    for col in merged_df.columns:
        if col == left_on:
            continue  # Skip the grouping column
        elif col in group_columns:
            # Aggregate group columns by extending lists
            agg_dict[col] = aggregate_lists
        elif col in left_columns:
            # For original left columns, take first value
            agg_dict[col] = "first"

    # Perform aggregation
    df_result = merged_df.groupby(left_on, as_index=False).agg(agg_dict)

    # Reorder columns to match original left file order, then add new group columns
    result_columns = [col for col in left_columns if col in df_result.columns]
    new_group_cols = [col for col in group_columns if col not in left_columns]
    final_column_order = result_columns + new_group_cols
    df_result = df_result[final_column_order]

    print(f"\nResult DataFrame shape: {df_result.shape}")
    print(f"Columns: {list(df_result.columns)}")

    # Save merged file
    print(f"\nSaving to {output_file}...")
    df_result.to_csv(output_file, index=False)
    print("Done!")

    # Print summary statistics
    print("\n--- Merge Summary ---")
    print(f"Right file rows: {len(df_right)}")
    print(f"Left file rows: {len(df_left)}")
    print(f"Output rows: {len(df_result)}")
    print(f"Group columns merged: {len(group_columns)}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge grouping columns from one metadata file to another"
    )

    parser.add_argument(
        "left_file", type=str, help="Path to target CSV file to add groupings to"
    )

    parser.add_argument(
        "right_file", type=str, help="Path to source CSV file with grouping information"
    )

    parser.add_argument(
        "output_file", type=str, help="Path to save the merged CSV file"
    )

    parser.add_argument(
        "--right-on",
        type=str,
        default="pdb_id",
        help="Column name in right DataFrame to match on (default: pdb_id)",
    )

    parser.add_argument(
        "--left-on",
        type=str,
        default="pdb_id",
        help="Column name in left DataFrame to match on (default: pdb_id)",
    )

    parser.add_argument(
        "--group-columns",
        nargs="+",
        type=str,
        default=None,
        help="Explicit list of group columns to merge (default: auto-detect)",
    )

    parser.add_argument(
        "--group-pattern",
        type=str,
        default=r"^(group_id|seq_group_id|mmseqs_.+)$",
        help="Regex pattern to match group columns (default: matches group_id, seq_group_id, mmseqs_*)",
    )

    args = parser.parse_args()

    merge_groupings(
        right_file=args.right_file,
        left_file=args.left_file,
        output_file=args.output_file,
        right_on=args.right_on,
        left_on=args.left_on,
        group_columns=args.group_columns,
        group_pattern=args.group_pattern,
    )


if __name__ == "__main__":
    main()

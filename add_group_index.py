import pandas as pd


def add_group_index(
    df: pd.DataFrame,
    group_key: str,
    sort_by: list[str] | None = None,
    group_index_name: str = "group_index",
) -> pd.DataFrame:
    """
    Add group index to dataframe based on grouping key and optional sorting.

    Args:
        df: Input dataframe
        group_key: Column name to group by
        sort_by: List of columns to sort by within each group before assigning indices
        group_index_name: Name of the column to add with group indices

    Returns:
        Dataframe with added group index column
    """
    df_copy = df.copy()

    # Sort if requested
    if sort_by:
        df_copy = df_copy.sort_values(by=sort_by)

    # Add group index
    df_copy[group_index_name] = df_copy.groupby(group_key).cumcount() + 1

    return df_copy


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Add group index to entries based on grouping from reference file or inline grouping key"
    )
    parser.add_argument(
        "input_file", type=str, help="Path to the CSV file to add group indices to"
    )
    parser.add_argument(
        "output_file", type=str, help="Output CSV file with added group indices"
    )

    parser.add_argument(
        "--group_key",
        type=str,
        default="group_id",
        help="Column defining groups (default: group_id)",
    )

    parser.add_argument(
        "--group_index_name",
        type=str,
        default="group_index",
        help="Name of the group index column to add (default: group_index)",
    )

    parser.add_argument(
        "--sort_by",
        nargs="+",
        type=str,
        default=None,
        help="Columns to sort by within each group before assigning indices",
    )

    parser.add_argument(
        "--reference_file",
        type=str,
        default=None,
        help="Optional reference CSV file to get grouping from. If specified, will use group_key from this file and join with input_file on join_key",
    )

    parser.add_argument(
        "--join_key",
        type=str,
        default="target_id",
        help="Column to join on when using reference_file (default: target_id)",
    )

    parser.add_argument(
        "--drop",
        nargs="+",
        type=str,
        default=None,
        help="Columns to drop before saving",
    )

    parser.add_argument(
        "--sort_by_after",
        nargs="+",
        type=str,
        default=None,
        help="Columns to sort by after adding group indices",
    )

    parser.add_argument(
        "--rename",
        nargs=2,
        action="append",
        metavar=("old_name", "new_name"),
        help="Rename columns after processing, applied in order",
    )

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)

    # Read input file
    df = pd.read_csv(input_file)

    # If reference file is provided, get grouping from it
    if args.reference_file:
        reference_file = Path(args.reference_file)
        ref_df = pd.read_csv(reference_file)

        # Extract only the join key and group key from reference
        ref_cols = set(
            [args.join_key, args.group_key] + (args.sort_by if args.sort_by else [])
        )
        ref_cols = [col for col in ref_cols if col in ref_df.columns]
        ref_df = ref_df[ref_cols].drop_duplicates()
    else:
        # Usse input df as ref df
        # Check if all key columns are present
        missing_cols = [
            col for col in [args.join_key, args.group_key] if col not in df.columns
        ]
        if missing_cols:
            raise ValueError(
                f"Missing columns in input file for grouping: {missing_cols}"
            )
        ref_df = df[[args.join_key, args.group_key]].drop_duplicates()

    # Add group index
    ref_df = add_group_index(
        ref_df,
        group_key=args.group_key,
        sort_by=args.sort_by,
        group_index_name=args.group_index_name,
    )
    # Merge back to original dataframe to get group indices
    df = df.merge(
        ref_df[[args.join_key, args.group_key, args.group_index_name]],
        on=[args.join_key],
    )

    if args.rename:
        # Do the sequential rename instead of one shot to avoid conflicts, renames are applied in order
        for old_name, new_name in args.rename:
            print(f"Renaming column: {old_name} -> {new_name}")
            # drop exising columns to avoid conflicts
            df = df.drop(columns=[new_name], errors="ignore")
            df = df.rename(columns={old_name: new_name})

    if args.sort_by_after:
        df = df.sort_values(by=args.sort_by_after)

    to_drop = args.drop if args.drop else []
    to_drop.extend([col for col in df.columns if col.endswith("_")])
    if to_drop:
        df = df.drop(columns=to_drop)

    # Put common columns at the front if present
    front_cols = ["pdb_id", "chain_id", "target_id", "auth_chain_id", "sequence"]
    front_cols = [col for col in front_cols if col in df.columns]
    other_cols = [col for col in df.columns if col not in front_cols]
    df = df[front_cols + other_cols]

    # Save output
    df.to_csv(output_file, float_format="%.3f", index=False)
    print(f"Added group indices to {len(df)} rows")
    print(f"Output saved to {output_file}")

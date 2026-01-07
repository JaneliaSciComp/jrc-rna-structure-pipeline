if __name__ == "__main__":
    import argparse
    from pathlib import Path
    import pandas as pd

    parser = argparse.ArgumentParser(
        description="Split sequences into two files based on temporal cutoff"
    )
    parser.add_argument(
        "grouped_sequences",
        type=str,
        help="Path to the CSV file with grouped sequences",
    )
    parser.add_argument(
        "split_file", type=str, help="File to split by matching PDB ids"
    )

    parser.add_argument(
        "--group_key",
        type=str,
        default=None,
        help="Column to group ids from GROUPED_SEQUENCES file  if not already grouped",
    )

    parser.add_argument(
        "--reference_key",
        type=str,
        default="target_id",
        help="Column with PDB/chain IDs in the GROUPED_SEQUENCES file",
    )

    parser.add_argument(
        "--id_column",
        type=str,
        default="target_id",
        help="Column with PDB/chain IDs to match",
    )

    parser.add_argument(
        "--before_prefix",
        type=str,
        default="train",
        help="Output prefix for sequences before the cutoff date",
    )

    parser.add_argument(
        "--after_prefix",
        type=str,
        default="test",
        help="Output prefix for sequences after the cutoff date",
    )

    parser.add_argument(
        "--cutoff_date",
        type=str,
        required=True,
        help="Cutoff date (YYYY-MM-DD) to split sequences",
    )

    args = parser.parse_args()
    cutoff_date = pd.to_datetime(args.cutoff_date)

    data = pd.read_csv(
        args.grouped_sequences, parse_dates=["temporal_cutoff"], keep_default_na=False
    )

    if "all_pdb_ids" not in data.columns:
        if not args.group_key:
            print(
                f"Group key needs to be provided if 'all_pdb_ids' are not present in {args.grouped_sequences}"
            )
            exit(1)

        if args.group_key not in data.columns:
            print(
                f"Group key '{args.group_key}' not found in {args.grouped_sequences} columns"
            )
            exit(1)

        # Group by the specified key to get all PDB ids per group and take the oldest cutoff date
        # Get oldest date per group

        cutoff_by_group = data.groupby(args.group_key)["temporal_cutoff"].min()
        # Assign to rows
        data = data.merge(
            cutoff_by_group.rename("group_temporal_cutoff"),
            left_on=args.group_key,
            right_index=True,
        )
        data = data[[args.reference_key, "group_temporal_cutoff"]].copy()
        data = data.rename(columns={"group_temporal_cutoff": "temporal_cutoff"})
        data = data.groupby(args.reference_key)["temporal_cutoff"].min().reset_index()
    else:
        # We already have grouped data
        if args.reference_key is not None:
            print("Warning: reference_key is ignored if  all_pdb_ids column is present")
        data = data[["all_pdb_ids", "temporal_cutoff"]].copy()
        # Split to row per pdb_id
        data = data.assign(all_pdb_ids=data["all_pdb_ids"].str.split(";")).explode(
            "all_pdb_ids"
        )
        data = data.rename(columns={args.reference_key: "all_pdb_ids"})

    ids_before = data[data["temporal_cutoff"] < cutoff_date][
        args.reference_key
    ].unique()
    ids_after = data[data["temporal_cutoff"] >= cutoff_date][
        args.reference_key
    ].unique()

    split_file = Path(args.split_file)

    data_to_split = pd.read_csv(split_file)
    data_before = data_to_split[data_to_split[args.id_column].isin(ids_before)]
    data_after = data_to_split[data_to_split[args.id_column].isin(ids_after)]

    data_before.to_csv(f"{args.before_prefix}_{split_file.name}", index=False)
    data_after.to_csv(f"{args.after_prefix}_{split_file.name}", index=False)

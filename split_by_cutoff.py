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

    data = pd.read_csv(args.grouped_sequences, parse_dates=["temporal_cutoff"])

    ids_before = (
        data[data["temporal_cutoff"] < cutoff_date]["all_pdb_ids"]
        .str.split(";", expand=True)
        .stack()
        .unique()
    )
    ids_after = (
        data[data["temporal_cutoff"] >= cutoff_date]["all_pdb_ids"]
        .str.split(";", expand=True)
        .stack()
        .unique()
    )

    split_file = Path(args.split_file)

    data_to_split = pd.read_csv(split_file)
    data_before = data_to_split[data_to_split[args.id_column].isin(ids_before)]
    data_after = data_to_split[data_to_split[args.id_column].isin(ids_after)]

    data_before.to_csv(f"{args.before_prefix}_{split_file.name}", index=False)
    data_after.to_csv(f"{args.after_prefix}_{split_file.name}", index=False)

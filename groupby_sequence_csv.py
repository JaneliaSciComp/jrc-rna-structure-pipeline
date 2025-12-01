from groupby_sequence import collect_grouped_idx
import pandas as pd


def groupby_sequence_csv(
    input_csv: str,
    output_csv: str,
    min_match_ratio: float = 0.9,
):
    """
    Groups sequences in a CSV file by similarity and saves the grouped data to a new CSV file.

    Args:
        input_csv (str): Path to the input CSV file containing sequences.
        output_csv (str): Path to the output CSV file to save grouped sequences.
        min_match_ratio (float): Minimum match ratio for grouping sequences.
    """
    # Load the CSV file
    data = pd.read_csv(input_csv, keep_default_na=False)

    # Group sequences
    grouped_data, seq_to_id = collect_grouped_idx(
        data["sequence"], min_match_ratio=min_match_ratio
    )
    # Get id to group mapping
    # Match group_idx to group name
    group_names = {
        group_idx: "_".join(
            data.loc[group_idx, ["pdb_id", "chain_id"]].astype(str).tolist()
        )
        for group_idx in grouped_data.keys()
    }

    data["group_id"] = ""
    data["target_id"] = data["pdb_id"] + "_" + data["chain_id"]
    for group_idx, entries in grouped_data.items():
        indexes = list(entries.keys())
        data.loc[indexes, "group_id"] = group_names[group_idx]

    data.to_csv(output_csv, index=False)


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Group sequences in a CSV file by similarity."
    )
    parser.add_argument("input_csv", help="Path to the input CSV file with sequences")
    parser.add_argument(
        "output_csv", help="Path to the output CSV file for grouped data"
    )
    parser.add_argument(
        "--min_match_ratio",
        type=float,
        default=0.9,
        help="Minimum match ratio for grouping sequences",
    )
    args = parser.parse_args()
    input_csv = Path(args.input_csv)
    output_csv = Path(args.output_csv)

    groupby_sequence_csv(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        min_match_ratio=args.min_match_ratio,
    )

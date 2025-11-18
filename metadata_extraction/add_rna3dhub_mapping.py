from pathlib import Path
import pandas as pd


def get_rna3dhub_mapping(rna3dhub_csv: str | Path) -> pd.DataFrame:
    """
    Adds RNA3DHub mapping information to the given DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame containing RNA structure data.

    Returns:
        pd.DataFrame: DataFrame with added RNA3DHub mapping information.
    """
    df = pd.read_csv(rna3dhub_csv)
    # data format
    # "NR_all_81027.2","9MME|1|k+9MME|1|o+9MME|1|s+9MME|1|w","9MME|1|k+9MME|1|o+9MME|1|s+9MME|1|w,9MME|1|U+9MME|1|Y+9MME|1|c+9MME|1|g"
    df.columns = ["rna3dhub_id", "representative", "members"]
    # Create mapping from pdbid_chain to rna3dhub_id
    mapping = {}
    exploded = df.assign(member_list=df["members"].str.split(r",|\+")).explode(
        "member_list"
    )
    parts = exploded["member_list"].str.split("|")
    pdb_chain = parts.str[0] + "_" + parts.str[2]

    mapping = dict(zip(pdb_chain, exploded["rna3dhub_id"]))
    return mapping


def assign_rna3dhub_mapping(df: pd.DataFrame, mapping: dict) -> pd.DataFrame:
    """
    Adds RNA3DHub mapping information to the given DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame containing RNA structure data.
        rna3dhub_csv (str | Path): Path to the RNA3DHub CSV mapping file.

    Returns:
        pd.DataFrame: DataFrame with added RNA3DHub mapping information.
    """
    # Use auth_asym_id, this is what rna3dhub uses
    auth_target_id = df["pdb_id"] + "_" + df["auth_chain_id"]
    df["rna3dhub_id"] = auth_target_id.map(mapping)

    return df


def process_csv_files(
    input_csv: str | Path,
    rna3dhub_csv: str | Path,
    output_csv: str | Path,
):
    """
    Processes the input CSV file to add RNA3DHub mapping information and saves the result.

    Args:
        input_csv (str | Path): Path to the input CSV file.
        rna3dhub_csv (str | Path): Path to the RNA3DHub CSV mapping file.
        output_csv (str | Path): Path to save the output CSV file.
    """
    df = pd.read_csv(input_csv)
    mapping = get_rna3dhub_mapping(rna3dhub_csv)
    df = assign_rna3dhub_mapping(df, mapping)
    df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Add RNA3DHub mapping information to RNA sequences data."
    )
    parser.add_argument(
        "input_csv", type=str, help="Path to the input CSV file with RNA sequences"
    )

    parser.add_argument(
        "rna3dhub_csv", type=str, help="Path to the RNA3DHub CSV mapping file"
    )

    parser.add_argument(
        "--output_csv",
        type=str,
        help="Path to save the output CSV file (default: {input_csv}_with_rna3dhub.csv)",
        default=None,
    )
    args = parser.parse_args()
    rna3dhub_csv = Path(args.rna3dhub_csv)
    input_csv = Path(args.input_csv)
    if args.output_csv is None:
        output_csv = input_csv.with_name(input_csv.stem + "_with_rna3dhub.csv")
    else:
        output_csv = Path(args.output_csv)

    process_csv_files(
        input_csv=input_csv,
        rna3dhub_csv=rna3dhub_csv,
        output_csv=output_csv,
    )

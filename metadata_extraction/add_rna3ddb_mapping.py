from pathlib import Path
import json
import pandas as pd


def get_rna3ddb_mapping(cluster_json: str | Path) -> dict:
    """
    Build a mapping from '<pdb_id>_<auth_chain_id>' to RNA3DDB component id using cluster.json.

    The cluster.json is expected to be a dict where top-level keys are component ids
    (e.g., "component_1") and each value is a dict whose keys are members like
    "8h8e_G", "8qoa_a", etc.

    Mapping is case-insensitive: keys are normalized to lowercase.

    Args:
        cluster_json (str | Path): Path to the cluster.json file.

    Returns:
        dict: Mapping from '<pdb_id>_<auth_chain_id>' (lowercase) to component id (rna3ddb_id).
    """
    cluster_json = Path(cluster_json)
    with cluster_json.open("r") as f:
        data = json.load(f)

    mapping_components: dict[str, str] = {}
    mapping_clusters: dict[str, str] = {}
    for component_id, clusters in data.items():
        if not isinstance(clusters, dict):
            # Defensive: skip unexpected shapes
            continue
        for cluster_id, member_keys in clusters.items():
            if not member_keys:
                continue
            for chain_id in member_keys:
                pdb_id, chain_id = chain_id.split("_")
                pdb_id = pdb_id.upper()
                # Chain id is case-sensitive
                id = f"{pdb_id}_{chain_id}"
                mapping_clusters[id] = cluster_id
                mapping_components[id] = component_id

    return mapping_components, mapping_clusters


def assign_rna3ddb_mapping(
    df: pd.DataFrame, mapping: dict, output_column: str
) -> pd.DataFrame:
    """
    Adds RNA3DDB component mapping information to the given DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame containing RNA structure data.
        mapping (dict): Mapping from '<pdb_id>_<auth_chain_id>' (lowercase) to component id.

    Returns:
        pd.DataFrame: DataFrame with added RNA3DDB mapping information in 'rna3ddb_id'.
    """
    # Build lookup key using auth_asym_id convention and normalize to lowercase.
    auth_target_id = df["pdb_id"] + "_" + df["auth_chain_id"]
    df[output_column] = auth_target_id.map(mapping)

    return df


def process_csv_files(
    input_csv: str | Path,
    cluster_json: str | Path,
    output_csv: str | Path,
):
    """
    Processes the input CSV file to add RNA3DDB component mapping information and saves the result.

    Args:
        input_csv (str | Path): Path to the input CSV file.
        cluster_json (str | Path): Path to the cluster.json mapping file.
        output_csv (str | Path): Path to save the output CSV file.
    """
    df = pd.read_csv(input_csv)
    mapping_components, mapping_clusters = get_rna3ddb_mapping(cluster_json)
    df = assign_rna3ddb_mapping(
        df, mapping_components, output_column="rna3ddb_component_id"
    )
    df = assign_rna3ddb_mapping(
        df, mapping_clusters, output_column="rna3ddb_cluster_id"
    )

    df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Add RNA3DDB component mapping information to RNA sequences data using cluster.json."
    )
    parser.add_argument(
        "input_csv", type=str, help="Path to the input CSV file with RNA sequences"
    )

    parser.add_argument(
        "cluster_json",
        type=str,
        help="Path to the cluster.json mapping file (RNA3DDB components)",
    )

    parser.add_argument(
        "--output_csv",
        type=str,
        help="Path to save the output CSV file (default: {input_csv}_with_rna3ddb.csv)",
        default=None,
    )
    args = parser.parse_args()
    cluster_json = Path(args.cluster_json)
    input_csv = Path(args.input_csv)
    if args.output_csv is None:
        output_csv = input_csv.with_name(input_csv.stem + "_with_rna3ddb.csv")
    else:
        output_csv = Path(args.output_csv)

    process_csv_files(
        input_csv=input_csv,
        cluster_json=cluster_json,
        output_csv=output_csv,
    )

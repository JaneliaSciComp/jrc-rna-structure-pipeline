#!/usr/bin/env python3
"""
Process multiple rounds of clustering results and create an association CSV file.

This script takes clustering results from multiple rounds where each round clusters
the results from the previous round. It creates a CSV file showing the cluster
membership of each initial member across all clustering rounds.

Usage:
    python process_clustering_rounds.py \
        --cluster-files round1.tsv round2.tsv round3.tsv \
        --initial-members members.txt \
        --output associations.csv
"""

import argparse
import csv
import sys
import pandas as pd
from pathlib import Path
from typing import Dict, List


def load_clustering_round(filepath: Path) -> Dict[str, str]:
    """
    Load a clustering round from a TSV file.

    Args:
        filepath: Path to TSV file with format: cluster_id\tmember_id

    Returns:
        Dictionary mapping member_id -> cluster_id
    """
    member_to_cluster = {}

    with open(filepath, "r") as f:
        reader = csv.reader(f, delimiter="\t")

        # Skip header if present
        first_row = next(reader, None)
        if first_row and (
            first_row[0].lower() in ["cluster_id", "cluster"]
            or "cluster" in first_row[0].lower()
        ):
            # This was a header, continue to next row
            pass
        else:
            # This was data, process it
            if first_row and len(first_row) >= 2:
                cluster_id, member_id = first_row[0], first_row[1]
                member_to_cluster[member_id] = cluster_id

        # Process remaining rows
        for row in reader:
            if len(row) >= 2:
                cluster_id, member_id = row[0], row[1]
                member_to_cluster[member_id] = cluster_id

    return member_to_cluster


def load_initial_members(
    filepath: Path, member_column: str, cluster_column: str
) -> Dict[str, str]:
    """
    Load initial member IDs from a file (one per line).

    Args:
        filepath: Path to file containing member IDs
        column_id: Name of the column to use for member IDs
    Returns:
        List of member IDs
    """
    df = pd.read_csv(filepath)
    return (
        df[[member_column, cluster_column]]
        .dropna()
        .set_index(member_column)[cluster_column]
        .to_dict()
    )


def process_clustering_rounds(
    cluster_files: List[Path], initial_members_assignments: Dict[str, str]
) -> List[Dict[str, str]]:
    """
    Process multiple rounds of clustering.

    Args:
        cluster_files: List of paths to clustering TSV files in order
        initial_members: List of initial member IDs to track

    Returns:
        List of dictionaries, each containing member_id and cluster assignments
    """
    # Initialize tracking for each member
    member_assignments = {
        member_id: {"member_id": member_id, "cluster_round_0": cluster_id}
        for member_id, cluster_id in initial_members_assignments.items()
    }

    for round_idx, cluster_file in enumerate(cluster_files, start=1):
        print(f"Processing round {round_idx}: {cluster_file}", file=sys.stderr)

        # Load this round's clustering
        member_to_cluster = load_clustering_round(cluster_file)

        cluster_col_name = f"cluster_round_{round_idx}"

        # For each initial member, propagate their cluster assignment
        for member_id in initial_members_assignments.keys():
            # Determine what to look up in this round's clustering
            prev_col = f"cluster_round_{round_idx - 1}"
            if prev_col in member_assignments[member_id]:
                lookup_id = member_assignments[member_id][prev_col]
            else:
                # No assignment in previous round, skip
                continue

            # Find cluster assignment for this lookup_id
            if lookup_id in member_to_cluster:
                cluster_id = member_to_cluster[lookup_id]
                member_assignments[member_id][cluster_col_name] = cluster_id

        # Count how many unique clusters for next round
        unique_clusters = set()
        for member_id in initial_members_assignments.keys():
            if cluster_col_name in member_assignments[member_id]:
                unique_clusters.add(member_assignments[member_id][cluster_col_name])

        print(
            f"  Assigned {len([m for m in member_assignments.values() if cluster_col_name in m])} members to {len(unique_clusters)} clusters",
            file=sys.stderr,
        )

    return list(member_assignments.values())


def write_associations_csv(associations: List[Dict[str, str]], output_path: Path):
    """
    Write associations to a CSV file.

    Args:
        associations: List of dictionaries with member_id and cluster assignments
        output_path: Path to output CSV file
    """
    if not associations:
        print("Warning: No associations to write", file=sys.stderr)
        return

    # Get all column names
    all_columns = set()
    for assoc in associations:
        all_columns.update(assoc.keys())

    # Sort columns: member_id first, then cluster rounds in order
    columns = []
    columns.extend(sorted(all_columns))

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(associations)

    print(f"Wrote {len(associations)} associations to {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Process multiple rounds of clustering results"
    )
    parser.add_argument(
        "--cluster-files",
        nargs="+",
        required=True,
        help="Clustering TSV files in order (round1.tsv round2.tsv ...)",
    )
    parser.add_argument(
        "--initial-members",
        required=True,
        help="File containing initial member IDs",
    )
    parser.add_argument(
        "--initial-member-id-column",
        default="target_id",
        help="Column name for initial member IDs",
    )
    parser.add_argument(
        "--initial-member-cluster-column",
        default="group_id",
        help="Column name for initial cluster IDs",
    )
    parser.add_argument("--output", required=True, help="Output CSV file path")

    parser.add_argument(
        "--rename-columns",
        nargs="+",
        default=[],
        help="Column names in output CSV, must match number of cluster rounds plus one for member_id and one for initial clusters",
    )
    args = parser.parse_args()

    # Convert to Path objects
    cluster_files = [Path(f) for f in args.cluster_files]
    initial_members_file = Path(args.initial_members)
    output_file = Path(args.output)

    # Validate input files
    for f in cluster_files:
        if not f.exists():
            print(f"Error: Clustering file not found: {f}", file=sys.stderr)
            sys.exit(1)

    if not initial_members_file.exists():
        print(
            f"Error: Initial members file not found: {initial_members_file}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Load initial members
    print(f"Loading initial members from {initial_members_file}", file=sys.stderr)
    initial_members = load_initial_members(
        initial_members_file,
        member_column=args.initial_member_id_column,
        cluster_column=args.initial_member_cluster_column,
    )
    print(f"  Loaded {len(initial_members)} initial members", file=sys.stderr)
    # Process clustering rounds
    associations = process_clustering_rounds(cluster_files, initial_members)

    # Rename columns to match requested names
    if args.rename_columns:
        expected_num_columns = len(cluster_files) + 2  # member_id + initial + rounds
        if len(args.rename_columns) != expected_num_columns:
            print(
                f"Error: Number of rename columns ({len(args.rename_columns)}) does not match expected ({expected_num_columns})",
                file=sys.stderr,
            )
            sys.exit(1)

        # Create mapping from old to new column names
        old_column_names = ["member_id"] + [
            f"cluster_round_{i}" for i in range(len(cluster_files) + 1)
        ]
        rename_mapping = dict(zip(old_column_names, args.rename_columns))

        # Apply renaming
        for assoc in associations:
            for old_col, new_col in rename_mapping.items():
                if old_col in assoc:
                    assoc[new_col] = assoc.pop(old_col)
    # Write output
    write_associations_csv(associations, output_file)

    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()

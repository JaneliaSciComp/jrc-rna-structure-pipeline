"""
Extracts x,y,z coordinates for each residue position in the sequence from metadata CSV.
For positions not present in the structure, fills with NaN.
"""

import argparse
from typing import List
import fsspec
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import parallel_process
from biotite.structure.atoms import AtomArray
from biotite.structure.io import pdbx


def extract_coordinates_for_chain(
    atomarray: AtomArray,
    chain_list: List[str],
    sequences: List[str],
    sequences_expected: List[str],
    atom_name: str = "C1'",
    pdbid: str = None,  # Optional PDB ID for logging
    skip_mismatch_validation: bool = False,
) -> pd.DataFrame:
    """
    Extract coordinates for each position in the sequence.

    Parameters
    ----------
    atomarray : AtomArray
        The atom array containing the structure.
    chain_id : str
        The chain ID to extract coordinates for.
    sequence : str
        The sequence string (e.g., "CAUGUGAC").
    atom_name : str
        The atom name to extract coordinates for (default: C1').

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: resname, resid, x, y, z for each position.
    """
    # Create expected sequence DataFrame
    expected_df = pd.DataFrame(
        {
            "resid": sum(
                [list(range(1, len(sequence) + 1)) for sequence in sequences], []
            ),
            "chain_id": sum(
                [
                    [chain_id] * len(sequence)
                    for chain_id, sequence in zip(chain_list, sequences)
                ],
                [],
            ),
            "expected_resname": sum([list(seq) for seq in sequences_expected], []),
            "resname": sum(
                [list(seq) for seq in sequences], []
            ),  # Will be updated with structure data
        }
    )
    # Filter to the specific chain
    chain_mask = np.isin(atomarray.chain_id, chain_list)
    chain_atoms = atomarray[chain_mask]

    if len(chain_atoms) == 0:
        # Return empty results for all positions if chain not found
        print("Warning PDB:", pdbid, "No atoms found for chain", chain_id)
        return pd.DataFrame()

    # Filter to target atom and create DataFrame from structure
    atom_mask = chain_atoms.atom_name == atom_name
    target_atoms = chain_atoms[atom_mask]

    if len(target_atoms) == 0:
        # No target atoms found, return empty results
        print(
            "Warning PDB:",
            pdbid,
            "No target atoms found for chain",
            chain_list,
            "and atom",
            atom_name,
        )
        return pd.DataFrame()

    # Create DataFrame from structure coordinates
    struct_df = pd.DataFrame(
        {
            "resid": target_atoms.res_id,
            "chain_id": target_atoms.chain_id,
            "struct_resname": target_atoms.res_name,
            "x": target_atoms.coord[:, 0],
            "y": target_atoms.coord[:, 1],
            "z": target_atoms.coord[:, 2],
        }
    )
    # Merge expected sequence with structure coordinates (left join to keep all expected positions)
    merged_df = expected_df.merge(struct_df, on=["chain_id", "resid"], how="left")

    # Validate that residue names match where we have coordinates
    if not skip_mismatch_validation:
        has_coords = merged_df["struct_resname"].notna()
        if has_coords.any():
            mismatches = merged_df.loc[
                has_coords
                & (merged_df["expected_resname"] != merged_df["struct_resname"])
            ]
            if len(mismatches) > 0:
                print(
                    f"Warning PDB: {pdbid} Residue name mismatches in chains {chain_list}:"
                )
                for _, row in mismatches.iterrows():
                    print(
                        f"  resid {row['resid']}: expected '{row['expected_resname']}', got '{row['struct_resname']}'"
                    )
    # Make resid sequential
    merged_df["resid"] = range(1, len(merged_df) + 1)
    return merged_df[["resname", "resid", "x", "y", "z"]]


def process_metadata_csv(
    metadata_csv: Path,
    cif_dir: Path,
    output_csv: Path,
    atom_name: str = "C1'",
    num_processes: int = None,
    file_pattern: str = "{pdb_id}.cif",
    multiple_chains: bool = False,
    skip_mismatch_validation: bool = False,
    limit: int = None,
):
    """
    Process metadata CSV and extract coordinates for each entry.

    Parameters
    ----------
    metadata_csv : Path
        Path to the metadata CSV file.
    cif_dir : Path
        Directory containing CIF files.
    output_csv : Path
        Path to the output CSV file.
    atom_name : str
        The atom name to extract coordinates for (default: C1').
    num_processes : int
        Number of parallel processes to use (default: all CPUs).
    """
    # Read metadata CSV
    metadata_df = pd.read_csv(metadata_csv, keep_default_na=False)

    # Group by pdb_id
    grouped = metadata_df.groupby("pdb_id")

    # Prepare arguments for parallel processing
    process_args = [
        (
            pdb_id,
            group_df,
            cif_dir,
            atom_name,
            file_pattern,
            multiple_chains,
            skip_mismatch_validation,
        )
        for pdb_id, group_df in grouped
    ]

    # Process in parallel
    if num_processes is None:
        import multiprocessing as mp

        num_processes = mp.cpu_count()
    if limit is not None:
        process_args = process_args[:limit]

    results = parallel_process(
        _process_pdb_group,
        process_args,
        num_processes=num_processes,
        desc="Extracting coordinates",
    )
    results = [res for res in results if res is not None and not res.empty]
    # Combine all results into a single DataFrame
    all_results = pd.concat(results, ignore_index=True)

    # Write output CSV
    all_results.to_csv(output_csv, index=False)


def _process_pdb_group(args) -> pd.DataFrame:
    """
    Process a single PDB group (all chains from one PDB file).

    Parameters
    ----------
    args : tuple
        (pdb_id, group_df, cif_dir, atom_name, file_pattern, multiple_chains, skip_mismatch_validation)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ID, pdb_id, chain_id, target_id, resname, resid, x_1, y_1, z_1
    """
    (
        pdb_id,
        group_df,
        cif_dir,
        atom_name,
        file_pattern,
        multiple_chains,
        skip_mismatch_validation,
    ) = args

    # Load the structure once for this PDB
    cif_path = cif_dir / file_pattern.format(pdb_id=pdb_id)
    atomarray = None

    if not cif_path.exists():
        print(f"Warning: CIF file not found for {pdb_id}, filling with NaN")
    else:
        try:
            with fsspec.open(str(cif_path), mode="rt", compression="infer") as fh:
                pdb_file = pdbx.CIFFile.read(fh)
                # Important to use use_author_fields=False to get consistent residue ids
                atomarray = pdbx.get_structure(
                    pdb_file, use_author_fields=False, model=1
                )
        except Exception as e:
            print(f"Warning: Failed to load {pdb_id}: {e}")

    # Process each chain in the group
    chain_results = []
    chain_key = "chain_ids" if multiple_chains else "chain_id"
    if multiple_chains:
        group_df["chain_list"] = group_df[chain_key].str.split(";")
        group_df["sequences"] = group_df["sequences"].str.split(";")
        group_df["sequences_expected"] = group_df["sequences_expected"].str.split(";")
    else:
        group_df["chain_list"] = group_df[chain_key].apply(lambda x: [x])
        group_df["sequences"] = group_df["sequence"].apply(lambda x: [x])
        group_df["sequences_expected"] = group_df["sequence_expected"].apply(
            lambda x: [x]
        )

    for _, row in group_df.iterrows():
        chain_list = row["chain_list"]
        sequences = row["sequences"]
        sequences_expected = [seq.split(":") for seq in row["sequences_expected"]]
        # They should match, if they don't don't use it
        # if len(sequence) != len(sequence_expected):
        #     print(
        #         f"Warning: Sequence length mismatch for {pdb_id} chain {chain_id}: "
        #         f"got {len(sequence)}, expected {len(sequence_expected)}"
        #     )
        #     continue

        if atomarray is not None:
            # Extract coordinates (returns DataFrame with resname, resid, x, y, z)
            chain_df = extract_coordinates_for_chain(
                atomarray,
                chain_list,
                sequences,
                sequences_expected,
                atom_name,
                pdbid=pdb_id,
                skip_mismatch_validation=skip_mismatch_validation,
            )
            target_id = row["target_id"]
            if chain_df is not None and not chain_df.empty:
                # Rename columns and add ID fields
                chain_df = chain_df.rename(columns={"x": "x_1", "y": "y_1", "z": "z_1"})
                chain_df["pdb_id"] = pdb_id
                chain_df["chain_id"] = chain_list[0] if not multiple_chains else None
                chain_df["target_id"] = target_id
                chain_df["ID"] = target_id + "_" + chain_df["resid"].astype(str)

                chain_results.append(chain_df)

    # Combine all chains and reorder columns
    if chain_results:
        result_df = pd.concat(chain_results, ignore_index=True)
        return result_df[
            [
                "ID",
                "pdb_id",
                "chain_id",
                "target_id",
                "resname",
                "resid",
                "x_1",
                "y_1",
                "z_1",
            ]
        ]
    else:
        # If no chains processed, return empty DataFrame with correct columns
        return pd.DataFrame(
            columns=[
                "ID",
                "pdb_id",
                "chain_id",
                "target_id",
                "resname",
                "resid",
                "x_1",
                "y_1",
                "z_1",
            ]
        )


def main():
    parser = argparse.ArgumentParser(
        description="Extract x,y,z coordinates for each residue position from metadata CSV"
    )
    parser.add_argument(
        "metadata_csv", type=Path, help="Path to the input metadata CSV file"
    )
    parser.add_argument(
        "cif_dir", type=Path, help="Directory containing CIF structure files"
    )
    parser.add_argument(
        "--file-pattern",
        type=str,
        default="{pdb_id}.cif",
        help="Glob pattern to match CIF files (default: '{pdb_id}.cif')",
    )
    parser.add_argument("output_csv", type=Path, help="Path to the output CSV file")
    parser.add_argument(
        "--atom-name",
        type=str,
        default="C1'",
        help="Atom name to extract coordinates for (default: C1')",
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=None,
        help="Number of parallel processes to use (default: all CPUs)",
    )
    parser.add_argument(
        "--multiple-chains",
        action="store_true",
        help="Process mulitichain entries (bioassemblies)",
    )
    parser.add_argument(
        "--skip-mismatch-validation",
        action="store_true",
        help="Skip validation of residue name mismatches between expected and structure sequences",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit the number of entries to process",
    )

    args = parser.parse_args()

    process_metadata_csv(
        args.metadata_csv,
        args.cif_dir,
        args.output_csv,
        args.atom_name,
        multiple_chains=args.multiple_chains,
        num_processes=args.num_processes,
        file_pattern=args.file_pattern,
        skip_mismatch_validation=args.skip_mismatch_validation,
        limit=args.limit,
    )
    print(f"Coordinates extracted to {args.output_csv}")


if __name__ == "__main__":
    main()

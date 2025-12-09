"""
Extracts x,y,z coordinates for each residue position in the sequence from metadata CSV.
For positions not present in the structure, fills with NaN.
"""

import argparse
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
    chain_id: str,
    sequence: str,
    sequence_expected: str,
    atom_name: str = "C1'",
    pdbid: str = None,  # Optional PDB ID for logging
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
            "resid": range(1, len(sequence) + 1),
            "expected_resname": list(sequence_expected),
            "resname": list(sequence),  # Will be updated with structure data
        }
    )
    # Filter to the specific chain
    chain_mask = atomarray.chain_id == chain_id
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
            chain_id,
            "and atom",
            atom_name,
        )
        return pd.DataFrame()

    # Create DataFrame from structure coordinates
    struct_df = pd.DataFrame(
        {
            "resid": target_atoms.res_id,
            "struct_resname": target_atoms.res_name,
            "x": target_atoms.coord[:, 0],
            "y": target_atoms.coord[:, 1],
            "z": target_atoms.coord[:, 2],
        }
    )

    # Merge expected sequence with structure coordinates (left join to keep all expected positions)
    merged_df = expected_df.merge(struct_df, on="resid", how="left")

    # Validate that residue names match where we have coordinates
    has_coords = merged_df["struct_resname"].notna()
    if has_coords.any():
        mismatches = merged_df.loc[
            has_coords & (merged_df["expected_resname"] != merged_df["struct_resname"])
        ]
        if len(mismatches) > 0:
            print(f"Warning PDB: {pdbid} Residue name mismatches in chain {chain_id}:")
            for _, row in mismatches.iterrows():
                print(
                    f"  resid {row['resid']}: expected '{row['expected_resname']}', got '{row['struct_resname']}'"
                )
    return merged_df[["resname", "resid", "x", "y", "z"]]


def process_metadata_csv(
    metadata_csv: Path,
    cif_dir: Path,
    output_csv: Path,
    atom_name: str = "C1'",
    num_processes: int = None,
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
    metadata_df = pd.read_csv(metadata_csv)

    # Group by pdb_id
    grouped = metadata_df.groupby("pdb_id")

    # Prepare arguments for parallel processing
    process_args = [
        (pdb_id, group_df, cif_dir, atom_name) for pdb_id, group_df in grouped
    ]

    # Process in parallel
    if num_processes is None:
        import multiprocessing as mp

        num_processes = mp.cpu_count()

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
        (pdb_id, group_df, cif_dir, atom_name)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ID, pdb_id, chain_id, target_id, resname, resid, x_1, y_1, z_1
    """
    pdb_id, group_df, cif_dir, atom_name = args

    # Load the structure once for this PDB
    cif_path = cif_dir / f"{pdb_id}.cif"
    atomarray = None

    if not cif_path.exists():
        print(f"Warning: CIF file not found for {pdb_id}, filling with NaN")
    else:
        try:
            pdb_file = pdbx.CIFFile.read(str(cif_path))
            # Important to use use_author_fields=False to get consistent residue ids
            atomarray = pdbx.get_structure(pdb_file, use_author_fields=False, model=1)
        except Exception as e:
            print(f"Warning: Failed to load {pdb_id}: {e}")

    # Process each chain in the group
    chain_results = []

    for _, row in group_df.iterrows():
        chain_id = row["chain_id"]
        sequence = row["sequence"]
        sequence_expected = row["sequence_expected"].split(":")
        # They should match, if they don't don't use it
        if len(sequence) != len(sequence_expected):
            print(
                f"Warning: Sequence length mismatch for {pdb_id} chain {chain_id}: "
                f"got {len(sequence)}, expected {len(sequence_expected)}"
            )
            continue

        target_id = f"{pdb_id}_{chain_id}"

        if atomarray is not None:
            # Extract coordinates (returns DataFrame with resname, resid, x, y, z)
            chain_df = extract_coordinates_for_chain(
                atomarray,
                chain_id,
                sequence,
                sequence_expected,
                atom_name,
                pdbid=pdb_id,
            )
            if chain_df is not None and not chain_df.empty:
                # Rename columns and add ID fields
                chain_df = chain_df.rename(columns={"x": "x_1", "y": "y_1", "z": "z_1"})
                chain_df["pdb_id"] = pdb_id
                chain_df["chain_id"] = chain_id
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

    args = parser.parse_args()

    process_metadata_csv(
        args.metadata_csv,
        args.cif_dir,
        args.output_csv,
        args.atom_name,
        num_processes=args.num_processes,
    )
    print(f"Coordinates extracted to {args.output_csv}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Extract ligand SMILES from PDB structures.

For each CIF file, extracts unique non-polymer entities (ligands) and looks up
their SMILES strings from a reference file.
"""

from functools import partial
from pathlib import Path
import biotite.structure.io.pdbx as pdbx
import pandas as pd
import argparse


def _get_cif_column_array(category, column_name):
    """
    Safely get a column from a CIF category as an array.

    Parameters
    ----------
    category : CIFCategory
        The CIF category to access
    column_name : str
        The column name (without the leading underscore)

    Returns
    -------
    list or None
        The column values as a list, or None if not found
    """
    if category is None or column_name not in category:
        return None
    try:
        column = category[column_name]
        return column.as_array() if hasattr(column, "as_array") else list(column)
    except Exception:
        return None


def load_smiles_lookup(smiles_file):
    """
    Load SMILES lookup table from tab-delimited file.

    Parameters
    ----------
    smiles_file : Path
        Path to the SMILES file (tab-delimited, SMILES in column 1, comp_id in column 2)

    Returns
    -------
    dict
        Dictionary mapping comp_id to SMILES string
    """
    df = pd.read_csv(
        smiles_file,
        sep="\t",
        header=None,
        keep_default_na=False,
        names=["smiles", "comp_id", "name"],
        dtype=str,
    )
    return dict(zip(df["comp_id"], df["smiles"]))


def extract_ligand_smiles(mmcif_file, smiles_lookup):
    """
    Extract ligand information and SMILES from a CIF file.

    Parameters
    ----------
    mmcif_file : Path
        Path to the mmCIF file
    smiles_lookup : dict
        Dictionary mapping comp_id to SMILES string

    Returns
    -------
    list of dict
        List of records with pdb_id, ligand_comp_id, ligand_SMILES
    """
    pdb_id = mmcif_file.stem

    try:
        cif_file = pdbx.CIFFile.read(str(mmcif_file))
        block = cif_file.block
    except Exception as e:
        print(f"Error reading {mmcif_file}: {e}")
        return []

    # Get non-polymer entity information
    nonpoly_category = block.get("pdbx_entity_nonpoly", None)

    if nonpoly_category is None:
        return []

    # Extract entity names and comp_ids
    names = _get_cif_column_array(nonpoly_category, "name")
    comp_ids = _get_cif_column_array(nonpoly_category, "comp_id")

    if comp_ids is None or len(comp_ids) == 0:
        return []

    # Filter out water (HOH) and get unique comp_ids
    unique_ligands = {}
    for name, comp_id in zip(names, comp_ids):
        if comp_id and comp_id != "HOH" and comp_id not in unique_ligands:
            unique_ligands[comp_id] = name

    # Create records with SMILES lookup
    records = []
    for comp_id, name in unique_ligands.items():
        smiles = smiles_lookup.get(comp_id, None)
        records.append(
            {
                "pdb_id": pdb_id,
                "ligand_comp_id": comp_id,
                "ligand_name": name,
                "ligand_SMILES": smiles,
            }
        )

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Extract ligand SMILES from PDB structures"
    )
    parser.add_argument(
        "input_dir", type=str, help="Path to the directory with mmCIF files"
    )
    parser.add_argument("output_file", type=str, help="Path to output CSV file")
    parser.add_argument(
        "--smiles-file",
        type=str,
        default="data/Components-smiles-stereo-oe_20251221.smi",
        help="Path to tab-delimited SMILES lookup file (default: data/Components-smiles-stereo-oe_20251221.smi)",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Run in debug mode (no parallelization)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit the number of files processed (for debugging)",
    )
    parser.add_argument(
        "--multiprocessing",
        action="store_true",
        help="Use multiprocessing via python multiprocessing instead of Dask",
    )

    args = parser.parse_args()

    # Load SMILES lookup table
    smiles_file = Path(args.smiles_file)
    if not smiles_file.exists():
        print(f"Error: SMILES file not found: {smiles_file}")
        return

    print(f"Loading SMILES lookup table from {smiles_file}...")
    smiles_lookup = load_smiles_lookup(smiles_file)
    print(f"Loaded {len(smiles_lookup)} SMILES entries")
    print(list(smiles_lookup.items())[:10])
    # Get list of CIF files
    input_dir = Path(args.input_dir)
    file_list = list(input_dir.glob("*.cif"))
    if args.limit is not None:
        file_list = file_list[: args.limit]

    print(f"Processing {len(file_list)} mmCIF files...")

    if args.debug:
        # Process sequentially
        all_records = []
        for cif_file in file_list:
            print(f"Processing {cif_file}...")
            records = extract_ligand_smiles(cif_file, smiles_lookup)
            all_records.extend(records)

        df = pd.DataFrame(all_records)

    elif args.multiprocessing:
        # Use Python multiprocessing
        import sys

        sys.path.append(str(Path(__file__).parent.parent))
        from utils import parallel_process

        result = parallel_process(
            partial(extract_ligand_smiles, smiles_lookup=smiles_lookup),
            file_list,
            desc="Extracting ligand SMILES",
        )

        all_records = []
        for records in result:
            all_records.extend(records)

        df = pd.DataFrame(all_records)

    else:
        # Use Dask for parallel processing
        import dask
        import dask.bag as db
        from dask.diagnostics import ProgressBar

        dask.config.set(scheduler="processes", num_workers=22)

        b = db.from_sequence(file_list, npartitions=22)

        with ProgressBar():
            df = (
                b.map(partial(extract_ligand_smiles, smiles_lookup=smiles_lookup))
                .flatten()
                .to_dataframe()
            )

    # Sort and save
    if len(df) > 0:
        df = df.sort_values(by=["pdb_id", "ligand_comp_id"])

        # Combine multiple entries per PDB into semicolon-delimited strings
        df_grouped = (
            df.groupby("pdb_id")
            .agg(
                {
                    "ligand_comp_id": lambda x: ";".join(x),
                    "ligand_SMILES": lambda x: ";".join(
                        [str(s) if pd.notna(s) else "" for s in x]
                    ),
                }
            )
            .reset_index()
        )

        # Select only required columns for output
        output_cols = ["pdb_id", "ligand_comp_id", "ligand_SMILES"]
        df_output = df_grouped[output_cols]

        output_path = Path(args.output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_output.to_csv(output_path, index=False)

        print(f"\nResults saved to {output_path}")
        print(f"Total PDB entries: {len(df_output)}")
        print(f"Total individual ligands: {len(df)}")
        print(f"Unique ligands: {df['ligand_comp_id'].nunique()}")

        # Report ligands without SMILES
        missing_smiles = df[df["ligand_SMILES"].isna()]
        if len(missing_smiles) > 0:
            print(f"\nWarning: {len(missing_smiles)} ligands missing SMILES")
            print(
                f"Unique ligands without SMILES: {missing_smiles['ligand_comp_id'].nunique()}"
            )
    else:
        print("No ligand data extracted")


if __name__ == "__main__":
    main()

"""
Extracts pairwise distances between selected atoms
"""

from collections import defaultdict
import numpy as np
import itertools
from biotite.structure.atoms import AtomArray
from biotite.structure.geometry import distance, index_distance
from biotite.structure.residues import (
    get_atom_name_indices,
)
from biotite.structure.chains import get_chains
from biotite.structure.io import pdbx

atom_sets = {
    "RNA": ("C1'", "C1'"),
    "RNA_protein": ("C1'", "CA"),
    "RNA_DNA": ("C1'", "C1'"),
}


def compute_structuredness(atomarray: AtomArray, atom_set_key):
    """
    Compute pairwise distances between selected atom sets.

    Parameters
    ----------
    atomarray : AtomArray
        The atom array containing the structure.
    atom_set_key : str
        Key to select the atom sets from `atom_sets`.

    Returns
    -------
    structuredness: dict[str, dict[str, float]]
        A dictionary with structuredness metrics per chain.
    """
    # Get the indexes of the atoms per residue
    atom_indices = get_atom_name_indices(atomarray, ("C1'",)).flatten()
    # Get only ones with matching atoms
    atom_indices = atom_indices[atom_indices != -1]
    # Compute product of indexes for pairwise distance calculation
    atom_indices_pairs = np.array(np.meshgrid(atom_indices, atom_indices)).T.reshape(
        -1, 2
    )
    # Remove reverse pairs (i,j) and (j,i)
    atom_indices_pairs = atom_indices_pairs[
        atom_indices_pairs[:, 0] < atom_indices_pairs[:, 1]
    ]
    # Remove self pairs
    atom_indices_pairs = atom_indices_pairs[
        atom_indices_pairs[:, 0] != atom_indices_pairs[:, 1]
    ]
    # Remove duplicate pairs
    atom_indices_pairs = np.unique(atom_indices_pairs, axis=0)

    # Get residues and chains for filtering (we have only one atom per residue here)
    residue_pairs = atomarray.res_id[atom_indices_pairs]
    chain_pairs = atomarray.chain_id[atom_indices_pairs]
    same_chain_mask = chain_pairs[:, 0] == chain_pairs[:, 1]

    # Remove pairs that are from the same chain but closer than 4 residues apart
    residue_sequence_distance = np.abs(residue_pairs[:, 0] - residue_pairs[:, 1])
    valid_pairs_mask = ~(same_chain_mask & (residue_sequence_distance < 4))

    atom_indices_pairs = atom_indices_pairs[valid_pairs_mask]
    residue_pairs = residue_pairs[valid_pairs_mask]
    chain_pairs = chain_pairs[valid_pairs_mask]
    # Calculate distances
    distances = index_distance(atomarray, atom_indices_pairs)

    # Get closest distance per residue for inter and intra chain
    inter_chain_mask = chain_pairs[:, 0] != chain_pairs[:, 1]
    intra_chain_mask = chain_pairs[:, 0] == chain_pairs[:, 1]

    min_inter_distances_per_idx = []
    min_intra_distances_per_idx = []

    # Precompute a mapping from each idx to the rows in atom_indices_pairs where it appears
    idx_to_pair_rows = defaultdict(list)
    for row, (i, j) in enumerate(atom_indices_pairs):
        idx_to_pair_rows[i].append(row)
        idx_to_pair_rows[j].append(row)

    for idx in atom_indices:
        rows = idx_to_pair_rows[idx]
        if rows:
            rows = np.array(rows)
            inter_mask = inter_chain_mask[rows]
            intra_mask = intra_chain_mask[rows]
            inter_distances = distances[rows][inter_mask]
            intra_distances = distances[rows][intra_mask]
        else:
            inter_distances = np.array([])
            intra_distances = np.array([])

        min_inter_distance = (
            np.min(inter_distances) if len(inter_distances) > 0 else np.inf
        )
        min_intra_distance = (
            np.min(intra_distances) if len(intra_distances) > 0 else np.inf
        )
        min_inter_distances_per_idx.append(min_inter_distance)
        min_intra_distances_per_idx.append(min_intra_distance)

    min_inter_distances_per_idx = np.array(min_inter_distances_per_idx)
    min_intra_distances_per_idx = np.array(min_intra_distances_per_idx)
    any_distance_per_idx = np.minimum(
        min_inter_distances_per_idx, min_intra_distances_per_idx
    )
    chain_per_idx = atomarray.chain_id[atom_indices]
    # Provide structuredness metrics per chain
    results = defaultdict(dict)

    for chain in get_chains(atomarray):
        this_chain_mask = chain_per_idx == chain
        if this_chain_mask.sum() == 0:
            results[chain]["inter_chain_structuredness"] = np.nan
            results[chain]["intra_chain_structuredness"] = np.nan
            results[chain]["total_chain_structuredness"] = np.nan
        else:
            results[chain]["inter_chain_structuredness"] = np.mean(
                min_inter_distances_per_idx[this_chain_mask] < 12
            )
            results[chain]["intra_chain_structuredness"] = np.mean(
                min_intra_distances_per_idx[this_chain_mask] < 12
            )
            results[chain]["total_chain_structuredness"] = np.mean(
                any_distance_per_idx[this_chain_mask] < 12
            )

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="Path to the input PDB file")
    parser.add_argument("output_file", type=str, help="Output CSV file for distances")
    parser.add_argument(
        "--pdb_id",
        type=str,
        help="PDB ID of the structure to assign. Uses file name stem if not provided.",
        default=None,
    )
    parser.add_argument(
        "--atom-set",
        type=str,
        choices=atom_sets.keys(),
        default="RNA",
        help="Atom set for distance calculation",
    )
    args = parser.parse_args()

    # Load the structure
    pdb_file = pdbx.CIFFile.read(args.input_file)
    atom_array = pdbx.get_structure(pdb_file, use_author_fields=False, model=1)
    if args.pdb_id is None:
        pdb_id = args.input_file.split("/")[-1].split(".")[0]
    else:
        pdb_id = args.pdb_id

    # Compute pairwise distances
    results = compute_structuredness(atom_array, args.atom_set)
    # Save to CSV
    with open(args.output_file, "w") as f:
        f.write(
            "pdb_id,chain_id,inter_chain_structuredness,intra_chain_structuredness,total_structuredness\n"
        )
        for chain, metrics in results.items():
            f.write(
                f"{pdb_id},{chain},{metrics['inter_chain_structuredness']},{metrics['intra_chain_structuredness']},{metrics['total_chain_structuredness']}\n"
            )
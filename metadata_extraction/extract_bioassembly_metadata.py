#!/usr/bin/env python3
"""
Extract PDB metadata including keywords, composition, and assembly information.

Uses biotite instead of BioPython for CIF file parsing.
Instead of filtering, this script extracts analysis information as metadata
that can be stored and used for downstream filtering decisions.
"""

from functools import partial
import json
from pathlib import Path
from collections import defaultdict
import re
import biotite.structure.io.pdbx as pdbx
import pandas as pd
import numpy as np
import dask.multiprocessing
from tqdm import tqdm
from rna_reference import modified_to_unmodified_nakb
import fsspec


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


def get_fasta_description(
    pdb_id,
    representative_entity_id,
    all_chains,
    entity_id_to_description,
    auth_to_asym_chain,
):
    if len(all_chains) == 1:
        auth_chain = all_chains[0]
        asym_chain = auth_to_asym_chain.get(auth_chain, auth_chain)
        chainid_description = f"Chain {asym_chain}[auth {auth_chain}]"
    else:
        asym_chains = [
            auth_to_asym_chain.get(auth_chain, auth_chain) for auth_chain in all_chains
        ]
        chain_strings = [
            f"{asym}[auth {auth}]" for asym, auth in zip(asym_chains, all_chains)
        ]
        chainid_description = f"Chains {', '.join(chain_strings)}"

    entity_description = entity_id_to_description.get(representative_entity_id, "null")
    description = f"{pdb_id}_{representative_entity_id}|{chainid_description}|{entity_description}|"
    auth_to_asym_chain
    return description


def extract_bioassembly_metadata(metadata_group, input_dir, file_pattern):
    """Extract assembly information from bioassembly mmCIF file.

    Parameters
    ----------
    block : CIFBlock
        The mmCIF block to extract from
    entity_types : dict, optional
        Mapping of entity type names to (column_suffix, display_name) tuples.
        Defaults to DEFAULT_ENTITY_TYPES if None. These are used to identify which ones to count explicitly in bioassembly
    """
    pdb_id = metadata_group[0]
    metadata = metadata_group[1]  # .groupby("entity_id").first().reset_index()
    metadata_entity_to_chain = (
        metadata.groupby("entity_id")["auth_chain_id"].apply(list).to_dict()
    )

    metadata_chain_to_sequence = metadata.set_index("auth_chain_id")[
        "sequence"
    ].to_dict()

    mmcif_file = Path(input_dir) / file_pattern.format(pdb_id=pdb_id)
    # Use fsspec to transparently handle compressed mmCIF files
    print(f"Processing {mmcif_file}")
    with fsspec.open(str(mmcif_file), mode="rt", compression="infer") as fh:
        cif_file = pdbx.CIFFile.read(fh)
    block = cif_file.block
    # Check RNA entities
    entity_poly_category = block.get("entity_poly", None)

    # Get entities and how many copies of each in the assembly
    entity_ids = list(_get_cif_column_array(entity_poly_category, "entity_id"))
    types = _get_cif_column_array(entity_poly_category, "type")
    entity_sequences = _get_cif_column_array(
        entity_poly_category, "pdbx_seq_one_letter_code_can"
    )
    strands = np.char.split(
        _get_cif_column_array(entity_poly_category, "pdbx_strand_id"), sep=","
    )
    n_strands = [len(s) for s in strands]
    # Map entity id to the list of strands
    entity_to_strand = {eid: strand for eid, strand in zip(entity_ids, strands)}
    strand_to_entity = {}
    for eid, strand_list in entity_to_strand.items():
        for strand in strand_list:
            strand_to_entity[strand] = eid
    # Get enttity_id to description mapping using entity block in cif
    entity_category = block.get("entity", None)
    entity_description_ids = _get_cif_column_array(entity_category, "id")
    entity_descriptions = _get_cif_column_array(entity_category, "pdbx_description")
    entity_id_to_description = dict(zip(entity_description_ids, entity_descriptions))

    # Auth
    # Get mapping using the pdbx_poly_seq_scheme from pdb chain_id (asym_id) to auth_chain_id (pdb_strand_id)
    poly_scheme_category = block.get("pdbx_poly_seq_scheme", None)
    asym_ids = _get_cif_column_array(poly_scheme_category, "asym_id")
    pdb_strand_ids = _get_cif_column_array(poly_scheme_category, "pdb_strand_id")
    auth_to_asym_chain = dict(zip(pdb_strand_ids, asym_ids))

    # Map strands(chains) to auth_chain_id in metadata, remove redundancy by removing '-1','-2' suffixes
    strand_to_entity = {}
    unique_strands = {}
    for eid, strand_list in zip(entity_ids, strands):
        strand_set = [re.sub(r"-\d+$", "", strand) for strand in strand_list]
        unique_strands[eid] = strand_set

    # Map entity_id to sequences from metadata using auth_chain_id matching
    entity_to_sequence = {}
    skipped_entities = []
    any_match = False
    for eid, strand_set in unique_strands.items():
        has_match = False
        for strand in strand_set:
            seq = metadata_chain_to_sequence.get(strand, None)
            entity_to_sequence[eid] = seq
            if seq is not None:
                has_match = True
                any_match = True
                break
        if not has_match:
            skipped_entities.append(eid)

    if not any_match:
        print(
            f"Warning: No matching sequences found for any entities in bioassembly for {pdb_id}"
        )
        return None

    df = pd.DataFrame(
        list(entity_to_sequence.items()), columns=["entity_id", "sequence"]
    )
    df["n_strands"] = n_strands
    df["strands"] = strands
    df = df.sort_values(
        by="sequence"
    )  # Sort by sequence to have consistent order between entries

    # Group by sequence to get total count of each unique sequence in assembly
    sequence_group = df.groupby("sequence")
    sequence_counts = sequence_group["n_strands"].sum().to_dict()
    sequence_to_entity_ids = sequence_group["entity_id"].apply(list).to_dict()
    sequence_to_strands = (
        sequence_group["strands"].apply(lambda x: sum(x, [])).to_dict()
    )
    # Construct the stochiometry of the assembly based on unique sequences
    # Assign each unique sequence to one chain id (arbitrarily the first one)
    stoichiometry_list = []
    merged_sequence = ""
    processed_entity_ids = set()

    all_sequences = {}
    sequence_list = []
    auth_chain_ids = []
    asym_chain_ids = []

    for sequence, count in sequence_counts.items():
        entity_ids_for_seq = sequence_to_entity_ids[sequence]
        processed_entity_ids.update(entity_ids_for_seq)
        representative_entity_id = sorted(entity_ids_for_seq)[0]
        all_chains = sorted(sequence_to_strands[sequence])
        representative_chain_id = all_chains[0]
        stoichiometry_list.append(f"{representative_chain_id}:{count}")
        description = get_fasta_description(
            pdb_id,
            representative_entity_id,
            all_chains,
            entity_id_to_description,
            auth_to_asym_chain,
        )
        all_sequences[description] = sequence
        merged_sequence += sequence * count
        sequence_list.extend([sequence] * count)

        auth_chain_ids.extend(all_chains)
        asym_chain_ids.extend(
            [auth_to_asym_chain.get(chain_id) for chain_id in all_chains]
        )

    stoichiometry = ";".join(stoichiometry_list)

    # Create a all_sequences fasta representation, where we record all sequences in the assembly including proteins
    for entity_id in skipped_entities:
        entity_sequence = entity_sequences[entity_ids.index(entity_id)]
        entity_chains = entity_to_strand[entity_id]
        representative_chain_id = sorted(entity_chains)[0]
        sequence = entity_sequence.replace("\n", "").replace(" ", "").strip()
        description = get_fasta_description(
            pdb_id,
            entity_id,
            entity_chains,
            entity_id_to_description,
            auth_to_asym_chain,
        )
        all_sequences[description] = sequence

    assembly_info = {
        "pdb_id": pdb_id,
        "target_id": pdb_id,
        "stoichiometry": stoichiometry,
        "sequence": merged_sequence,
        "sequences": ";".join(sequence_list),
        "sequences_expected": ";".join([":".join(seq) for seq in sequence_list]),
        "all_sequences": "\n".join(
            f">{desc}\n{seq}" for desc, seq in all_sequences.items()
        ),
        "auth_chain_ids": ";".join(auth_chain_ids),
        "chain_ids": ";".join(asym_chain_ids),
        "temporal_cutoff": metadata["temporal_cutoff"].values[0],
        "title": metadata["title"].values[0],
    }
    return assembly_info


if __name__ == "__main__":
    import sys

    sys.path.append(str(Path(__file__).parent.parent))
    from utils import parallel_process
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract metadata from PDB structures for filtering analysis using biotite"
    )
    parser.add_argument(
        "input_dir", type=str, help="Path to the directory with mmCIF files"
    )
    parser.add_argument(
        "metadata_file",
        type=str,
        help="Path to input metadata file",
    )

    parser.add_argument(
        "output_file",
        type=str,
        help="Path to output file with extracted bioassemblymetadata",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run in debug mode (no parallelization)",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.cif",
        help="Pattern to match mmCIF files in the input directory",
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

    input_dir = Path(args.input_dir)
    # front_cols = [
    #     "target_id",
    #     "pdb_id",
    #     "chain_id",
    #     "auth_chain_id",
    #     "sequence",
    # ]
    # Get metadata for each chain
    metadata = pd.read_csv(args.metadata_file, keep_default_na=False, dtype=str)
    metadata_groups = metadata.groupby("pdb_id")
    if args.limit is not None:
        metadata_groups = list(metadata_groups)[: args.limit]
    print(
        f"Processing {len(metadata_groups)} mmCIF files for bioassembly metadata extraction"
    )
    if not args.debug and args.multiprocessing:
        result = parallel_process(
            partial(
                extract_bioassembly_metadata,
                input_dir=input_dir,
                file_pattern=args.pattern,
            ),
            metadata_groups,
            desc="Extracting metadata (biotite)",
        )
        r = pd.concat([pd.DataFrame(r) for r in result], axis=0)
        r.to_csv(args.output_file, date_format="%Y-%m-%d", index=False)
    else:
        result = pd.DataFrame()
        result = []
        for metadata_group in tqdm(metadata_groups):
            metadata = extract_bioassembly_metadata(
                metadata_group, input_dir=input_dir, file_pattern=args.pattern
            )
            if metadata is not None:
                result.append(metadata)
        result = pd.DataFrame(result).sort_values(by=["temporal_cutoff", "target_id"])
        result.to_csv(args.output_file, date_format="%Y-%m-%d", index=False)

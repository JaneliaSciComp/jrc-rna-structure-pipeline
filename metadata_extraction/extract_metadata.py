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
from rna_reference import modified_to_unmodified_nakb


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


def extract_keywords(mmcif_file, metadata, keywords):
    """Checks if keyword is in description/title/keywords fields."""
    metadata = metadata.lower()
    keywords = {kw[0].lower(): kw[1] if kw[1] else kw[0] for kw in keywords}
    result = {}
    for keyword, store_name in keywords.items():
        has_keyword = keyword.lower() in metadata
        result[store_name] = has_keyword
        if not has_keyword:
            # Check mmCIF file as well
            with open(mmcif_file) as f:
                mmcif_data = f.read().lower()
            has_keyword = keyword in mmcif_data
            result[store_name] = has_keyword
    return result


def extract_composition_metadata(block):
    """Extract composition information from mmCIF file."""
    composition_info = {
        "rna_fraction": None,
        "rna_mass": None,
        "total_mass": None,
    }

    # Get entity information
    entity_category = block.get("entity", None)
    entity_poly_category = block.get("entity_poly", None)

    if entity_category is None:
        raise ValueError("No entity category found")

    # Extract entity information
    entity_ids = _get_cif_column_array(entity_category, "id")
    entity_formula_weights = _get_cif_column_array(entity_category, "formula_weight")
    entity_pdbx_number_of_molecules = _get_cif_column_array(
        entity_category, "pdbx_number_of_molecules"
    )

    if len(entity_ids) == 0:
        raise ValueError("No entity IDs found")

    # Calculate mass by entity ID
    mass_by_id = {}
    for k, weight, copies in zip(
        entity_ids, entity_formula_weights, entity_pdbx_number_of_molecules
    ):
        weight_float = float(weight) if weight and weight != "?" else None
        copies_int = int(copies) if copies and copies != "?" else None
        if weight_float and copies_int:
            mass_by_id[k] = weight_float * copies_int

    # Get entity polymer types and calculate mass by type
    if entity_poly_category is not None:
        entity_poly_entity_ids = _get_cif_column_array(
            entity_poly_category, "entity_id"
        )
        entity_poly_types = _get_cif_column_array(entity_poly_category, "type")
        # These are chain ids

        mass_by_type = defaultdict(float)
        for entity_id, poly_type in zip(entity_poly_entity_ids, entity_poly_types):
            if entity_id in mass_by_id:
                mass_by_type[poly_type] += mass_by_id[entity_id]

        # Calculate RNA fraction
        rna_mass = mass_by_type.get("polyribonucleotide", 0.0)
        total_mass = sum(mass_by_type.values())

        composition_info["rna_mass"] = float(rna_mass)
        composition_info["total_mass"] = float(total_mass)

        if total_mass > 0:
            composition_info["rna_fraction"] = float(rna_mass) / float(total_mass)

    return composition_info


def extract_assembly_metadata(block):
    """Extract assembly information from mmCIF file."""
    assembly_info = {
        "num_assemblies": 0,
        "assembly_defined": False,
        "has_symmetry_operators": False,
        "rna_entities": 0,
        "rna_copies": 0,
    }
    # Get assembly generation information
    assembly_gen_category = block.get("pdbx_struct_assembly_gen", None)

    if assembly_gen_category is None:
        return assembly_info

    # Extract assembly information
    assemblies = _get_cif_column_array(assembly_gen_category, "assembly_id")
    oper_expressions = _get_cif_column_array(assembly_gen_category, "oper_expression")
    assemblies_chain_ids = _get_cif_column_array(assembly_gen_category, "asym_id_list")

    if len(assemblies) == 0:
        return assembly_info

    assembly_info["assembly_defined"] = True
    assembly_info["num_assemblies"] = len(set(assemblies))

    # Check for symmetry operators
    has_symmetry = any(op != "1" for op in oper_expressions)
    assembly_info["has_symmetry_operators"] = has_symmetry

    # # Store assembly details
    # for asm_id, oper_expr, chain_list in zip(
    #     assemblies, oper_expressions or [], assemblies_chain_ids or []
    # ):
    #     assembly_info["assemblies"].append(
    #         {
    #             "assembly_id": str(asm_id),
    #             "operator_expression": str(oper_expr),
    #             "chains": [c.strip() for c in str(chain_list).split(",")]
    #             if chain_list
    #             else [],
    #         }
    #     )

    # Check RNA entities
    entity_poly_category = block.get("entity_poly", None)
    struct_asym_id_category = block.get("struct_asym", None)
    asym_id_to_entity_id = dict(
        zip(
            _get_cif_column_array(struct_asym_id_category, "id"),
            _get_cif_column_array(struct_asym_id_category, "entity_id"),
        )
    )

    # Get first assembly right now and check if we have RNA entities
    asym_ids = set(str(assemblies_chain_ids[0]).split(","))
    assembly_entity_ids = [asym_id_to_entity_id[asym_id] for asym_id in asym_ids]

    entity_ids = _get_cif_column_array(entity_poly_category, "entity_id")
    types = _get_cif_column_array(entity_poly_category, "type")
    entity_id_to_type = {eid: etype for eid, etype in zip(entity_ids, types)}
    rna_entities = [
        entity_id
        for entity_id in assembly_entity_ids
        if entity_id_to_type.get(entity_id) == "polyribonucleotide"
    ]
    assembly_info["rna_entities"] = len(set(rna_entities))
    assembly_info["rna_copies"] = len(rna_entities)
    return assembly_info


def extract_general_metadata(metadata):
    """
    Extract general information from metadata json file. Such as title and release date.
    """
    try:
        resolution = metadata["rcsb_entry_info"]["resolution_combined"][0]
    except KeyError:
        resolution = None
    general_info = {
        "temporal_cutoff": metadata["rcsb_accession_info"]["initial_release_date"],
        "resolution": resolution,
        "method": metadata["rcsb_entry_info"]["experimental_method"],
        "title": metadata["struct"]["title"],
    }
    return general_info


def extract_all_metadata(mmcif_file, keywords=[]):
    """Extract all metadata for a mmCIF file."""
    pdb_id = mmcif_file.stem
    cif_file = pdbx.CIFFile.read(str(mmcif_file))
    block = cif_file.block

    json_metadata_file = mmcif_file.with_name(mmcif_file.stem + ".metadata.json")
    with open(json_metadata_file) as f:
        metadata_str = f.read()
        metadata = json.loads(metadata_str)
    all_metadata = {
        "pdb_id": pdb_id,
        "mmcif_file": str(mmcif_file),
        "keywords": extract_keywords(mmcif_file, metadata_str, keywords),
        "composition": extract_composition_metadata(block),
        "assembly": extract_assembly_metadata(block),
        "general": extract_general_metadata(metadata),
    }

    # Flatten the metadata dictionary
    keywords = dict(**{f"keyword_{k}": v for k, v in all_metadata["keywords"].items()})
    composition = dict(
        **{f"composition_{k}": v for k, v in all_metadata["composition"].items()}
    )
    assembly = dict(**{f"assembly_{k}": v for k, v in all_metadata["assembly"].items()})
    flat = dict(
        pdb_id=all_metadata["pdb_id"],
        mmcif_file=all_metadata["mmcif_file"],
        **keywords,
        **composition,
        **assembly,
        **all_metadata["general"],
    )

    # Get records by chain
    chain_metadata = extract_chain_metadata(block)
    metadata_records = []
    for chain_id, metadata in chain_metadata.items():
        metadata.update(flat)
        metadata_records.append(metadata)
    return metadata_records


def extract_chain_metadata(block):
    """Extract per-chain metadata from mmCIF file."""
    chain_metadata = {}
    struct_asym_category = block.get("struct_asym", None)
    entity_poly_category = block.get("entity_poly", None)

    if struct_asym_category is None:
        return chain_metadata

    asym_ids = _get_cif_column_array(struct_asym_category, "id")
    entity_ids = _get_cif_column_array(struct_asym_category, "entity_id")
    chain_to_entity = dict(zip(asym_ids, entity_ids))

    entity_poly_entity_ids = _get_cif_column_array(entity_poly_category, "entity_id")
    entity_poly_types = _get_cif_column_array(entity_poly_category, "type")
    seqs = _get_cif_column_array(entity_poly_category, "pdbx_seq_one_letter_code_can")
    seqs_full = _get_cif_column_array(entity_poly_category, "pdbx_seq_one_letter_code")
    nonstandard_residues = _get_cif_column_array(entity_poly_category, "nstd_monomer")
    entity_to_type = dict(zip(entity_poly_entity_ids, entity_poly_types))
    entity_to_seq = dict(zip(entity_poly_entity_ids, seqs))
    entity_to_seq_full = dict(zip(entity_poly_entity_ids, seqs_full))
    entity_to_nonstandard = dict(zip(entity_poly_entity_ids, nonstandard_residues))
    entity_to_nonstandard = {k: (v == "yes") for k, v in entity_to_nonstandard.items()}

    for chain_id, entity_id in chain_to_entity.items():
        type_str = entity_to_type.get(entity_id, "unknown")
        if type_str == "polyribonucleotide":
            full_seq = entity_to_seq_full.get(entity_id, "")
            full_seq = full_seq.replace("\n", "").replace(" ", "").strip()
            # Split and use separator that is easier to process later, remove parentheses
            full_seq_list = re.findall(r"(?<=\()[^)]*(?=\))|[A-Z]", full_seq)
            full_seq = ":".join(full_seq_list)
            canonical_seq_tmp = entity_to_seq.get(entity_id, "")
            canonical_seq_tmp = (
                canonical_seq_tmp.replace("\n", "").replace(" ", "").strip()
            )
            # Split and use separator that is easier to process later, remove parentheses (but should not be any here)
            canonical_seq_list = re.findall(
                r"(?<=\()[^)]*(?=\))|[A-Z]", canonical_seq_tmp
            )
            # Check if all are single letter
            canonical_seq = ":".join(canonical_seq_list)
            # Makes sure the list are equal length
            assert len(full_seq_list) == len(canonical_seq_list)
            unmapped_canonical = False
            mapped_seq_list = canonical_seq_list.copy()

            if "X" in canonical_seq:
                unmapped_canonical = True
                # If residues are not mapped, use NAKB for X in canonical sequence
                for i, res in enumerate(canonical_seq_list):
                    if res == "X":
                        monomer_id = full_seq_list[i]
                        new_monomer_id = modified_to_unmodified_nakb.get(
                            monomer_id, "X"
                        )
                        mapped_seq_list[i] = new_monomer_id
            unmapped_nakb = "X" in mapped_seq_list
            mapped_seq = "".join(mapped_seq_list)
            is_mapped_single_letter = all(len(res) == 1 for res in mapped_seq_list)
            assert is_mapped_single_letter, (
                f"Mapped sequence has non-single letter residues: {mapped_seq_list}"
            )
            undefined_residues = "N" in mapped_seq_list
            unexpected_residues = (
                len(set(mapped_seq_list) - {"A", "U", "G", "C", "X", "N"}) > 0
            )

            chain_metadata[chain_id] = {
                "chain_id": chain_id,
                "entity_id": entity_id,
                "unmapped_canonical": unmapped_canonical,
                "unmapped_nakb": unmapped_nakb,
                "undefined_residues": undefined_residues,
                "unexpected_residues": unexpected_residues,
                "nonstandard_residues": entity_to_nonstandard.get(entity_id, False),
                "sequence": mapped_seq,
                "canonical_sequence": canonical_seq,
                "full_sequence": full_seq,
                # "type": entity_to_type.get(entity_id, "unknown"),
            }

    # Get mapping between label_asym_id and auth_asym_id and extract sequence
    # This are specified in _pdbx_poly_seq_scheme as asym_id and  pdb_strand_id as auth_asym_id
    poly_seq_scheme_category = block.get("pdbx_poly_seq_scheme", None)

    label_to_auth = {}
    if poly_seq_scheme_category is not None:
        label_asym_ids = _get_cif_column_array(poly_seq_scheme_category, "asym_id")
        mon_ids = _get_cif_column_array(poly_seq_scheme_category, "mon_id")
        pdb_mon_ids = _get_cif_column_array(poly_seq_scheme_category, "pdb_mon_id")
        auth_asym_ids = _get_cif_column_array(poly_seq_scheme_category, "pdb_strand_id")
        # Get unique pairs
        label_auth_asym_ids = np.array([label_asym_ids, auth_asym_ids])
        unique_pairs = np.unique(label_auth_asym_ids, axis=1)
        label_to_auth = dict(zip(unique_pairs[0], unique_pairs[1]))

        # Get sequence from mon_id
        for chain_id in chain_metadata.keys():
            chain_mask = label_asym_ids == chain_id
            sequence_expected = ":".join(mon_ids[chain_mask])
            sequence_observed = ":".join(pdb_mon_ids[chain_mask])

            # These are used for experiment
            chain_metadata[chain_id]["sequence_expected"] = sequence_expected
            # These are observed / modeled
            chain_metadata[chain_id]["sequence_observed"] = sequence_observed

            # Check if they match, if not flag
            sequence_match = sequence_expected == sequence_observed
            chain_metadata[chain_id]["missing_residues"] = not sequence_match

        for chain_id, metadata in chain_metadata.items():
            metadata["auth_chain_id"] = label_to_auth.get(chain_id, None)

    return chain_metadata


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
        "output_file",
        type=str,
        help="Path to output JSON file with extracted metadata",
    )
    parser.add_argument(
        "--format",
        choices=["json", "jsonl"],
        default="json",
        help="Output format: json (single file with array) or jsonl (one object per line)",
    )

    parser.add_argument(
        "--keyword",
        required=True,
        action="append",
        type=str,
        help="Keyword to filter by if two values provided per keyword then search for first but store as second value",
        nargs="+",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run in debug mode (no parallelization)",
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
    file_list = list(input_dir.glob("*.cif"))
    if args.limit is not None:
        file_list = file_list[: args.limit]

    if not args.debug:
        if args.multiprocessing:
            result = parallel_process(
                partial(extract_all_metadata, keywords=args.keyword),
                file_list,
                desc="Extracting metadata (biotite)",
            )
            r = pd.concat([pd.DataFrame(r) for r in result], axis=0)
            # Proper format
            r["temporal_cutoff"] = pd.to_datetime(
                r["temporal_cutoff"], errors="coerce"
            ).dt.strftime("%Y-%m-%d")

            r.set_index("pdb_id").to_csv(args.output_file, date_format="%Y-%m-%d")
        else:
            # Use Dask bag for parallel processing
            import dask.bag as db
            from dask.diagnostics import ProgressBar
            # from dask.distributed import Client, LocalCluster

            # cluster = LocalCluster(n_workers=20, threads_per_worker=1)
            # c = Client(cluster)
            dask.config.set(scheduler="processes", num_workers=22)
            print(f"Processing {len(file_list)} mmCIF files...")
            b = db.from_sequence(file_list, npartitions=22)

            b.map(
                partial(extract_all_metadata, keywords=args.keyword)
            ).flatten().to_dataframe().set_index("pdb_id").to_csv(
                args.output_file, date_format="%Y-%m-%d"
            )  # , engine="pyarrow", compression="gzip")

    else:
        result = pd.DataFrame()
        for cif_file in file_list:
            print(f"Processing {cif_file}...")
            metadata = extract_all_metadata(cif_file, keywords=args.keyword)
            result = pd.concat([result, pd.DataFrame(metadata)], axis=0)
        result = result.set_index("pdb_id")
        # Proper format
        result["temporal_cutoff"] = pd.to_datetime(
            result["temporal_cutoff"], errors="coerce"
        ).dt.strftime("%Y-%m-%d")
        result.to_csv(args.output_file, date_format="%Y-%m-%d")

    # print(f"Processing {len(file_list)} mmCIF files...")
    # print(args.keyword)
    # metadata_list = parallel_process(
    #     partial(extract_all_metadata, keywords=args.keyword),
    #     file_list,
    #     desc="Extracting metadata (biotite)",
    # )

    # output_path = Path(args.output_file)
    # output_path.parent.mkdir(parents=True, exist_ok=True)

    # if args.format == "json":
    #     with open(output_path, "w") as f:
    #         json.dump(metadata_list, f, indent=2)
    #     print(f"Metadata saved to {output_path}")
    # else:  # jsonl
    #     with open(output_path, "w") as f:
    #         for metadata in metadata_list:
    #             f.write(json.dumps(metadata) + "\n")
    #     print(f"Metadata saved to {output_path} (JSONL format)")

    # print(f"Total files processed: {len(metadata_list)}")

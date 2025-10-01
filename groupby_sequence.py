from collections import defaultdict


def groupby_sequence(data):
    grouped_data = {}
    seq_to_id = {}
    for i in range(len(data["sequence"])):
        seq = data["sequence"][i]
        xyz = data["xyz"][i]
        cif_file = data["data_cif_files"][i]
        pdb_id = cif_file.split("/")[-1].split(".")[0]
        if seq not in seq_to_id:
            seq_to_id[seq] = pdb_id
            ref_id = pdb_id
            grouped_data[ref_id] = defaultdict(list)
        else:
            ref_id = seq_to_id[seq]
        grouped_data[ref_id]["sequence"].append(seq)
        grouped_data[ref_id]["data_cif_files"].append(cif_file)
        grouped_data[ref_id]["pdb_ids"].append(pdb_id)
        grouped_data[ref_id]["xyz"].append(xyz)

    return grouped_data


if __name__ == "__main__":
    import argparse
    import pickle
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        help="Pickle file with data to be grouped by sequence",
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Output directory for the grouped data to be saved",
    )

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_file, "rb") as f:
        data = pickle.load(f)

    grouped_data = groupby_sequence(data)

    for group_id, entries in grouped_data.items():
        print(f"Group {group_id} has {len(entries['pdb_ids'])} entries.")
        output_file = output_dir / f"{group_id}_grouped.pkl"
        with open(output_file, "wb") as f:
            pickle.dump(entries, f)

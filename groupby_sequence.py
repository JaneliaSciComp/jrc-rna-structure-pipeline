from collections import defaultdict
import numpy as np

class PartialSequenceMapping(dict):
    """A dictionary that allows partial key matching.

    Retrieves the value for a key if the key is a substring of any existing key or contains existing keys
    with a minimum match ratio.
    """

    def __init__(self, *args, min_match_ratio=0.8, **kwargs):
        super().__init__(*args, **kwargs)
        self.min_match_ratio = min_match_ratio
        if self.min_match_ratio >= 1:
            self.get_match_item = self._get_exact_item
            self.__contains__ = super().__contains__

    def _get_exact_item(self, key):
        value = super().__getitem__(key)
        return key, 0, value

    def get_match_item(self, key):
        # First try exact match
        try:
            return self._get_exact_item(key)
        except KeyError:
            # Then try partial match with minimum ratio
            # Search for keys that are substrings or superstrings of the given key and return the position
            for k in self.keys():
                # Check length ratios first to avoid unnecessary computations
                ratio = len(key) / len(k)
                if (
                    ratio < 1 and ratio >= self.min_match_ratio
                ):  # key is potential substring
                    pos = k.find(key)
                    if pos != -1:
                        return k, pos, super().__getitem__(k)
                elif (
                    ratio > 1 and (1 / ratio) >= self.min_match_ratio
                ):  # key is potential superstring
                    pos = key.find(k)
                    if pos != -1:
                        return k, -pos, super().__getitem__(k)
        raise KeyError(
            f"No matching key found for '{key}' with min_match_ratio {self.min_match_ratio}"
        )

def get_message(
    msg,
    ref_len,
    L,
    ref_id,
    ref_pdb_id,
    entry_idx,
    pdb_id,
    position,
    seq,
    ref_seq,
    **kwargs,
):
    msg_info = (
        f"{msg}"
        f" {ref_len}, current sequence length: {L})"
        f" (ref_id: {ref_id} ({ref_pdb_id}), entry_idx: {entry_idx} ({pdb_id}))"
        f" (position: {position})"
        f" (sequence: {seq})"
        f" (reference sequence: {ref_seq})"
    )
    return msg_info

def collect_grouped_idx(sequences, min_match_ratio=0.9):
    """
    Collects grouped sequence IDs based on partial matching with a minimum match ratio.

    Args:
        sequences (list): A list of sequences to group.
        min_match_ratio (float): Minimum ratio for partial matching.

    Returns:
        dict: A dictionary mapping reference IDs to grouped sequence IDs and their positions.
    """
    grouped_ids = {}
    seq_to_id = PartialSequenceMapping(min_match_ratio=min_match_ratio)
    # Collect sequences and coordinates that group by partial matching
    for i in range(len(sequences)):
        seq = sequences[i]
        try:
            matched_seq, position, ref_id = seq_to_id.get_match_item(seq)
        except KeyError:
            # No match found, add new entry
            seq_to_id[seq] = i
            ref_id = i
            grouped_ids[ref_id] = {i: 0}
        else:
            # Match found, use existing group, if position is positive new sequence is substring
            # so keep current ref_id, else if negative new sequence is superstring, so update ref_id
            # and reference sequence

            # We have the following cases:
            if position >= 0 and position + len(seq) <= len(matched_seq):
                #    New sequence is substring of matched_seq, keep ref_id
                #    ....
                #     ..
                grouped_ids[ref_id][i] = position
            elif position < 0 and position + len(seq) >= len(matched_seq):
                #    New sequence is superstring of matched_seq, update ref_id to i
                #     ..
                #    ....
                current_group = grouped_ids.pop(ref_id)
                # Shift all existing positions by -position
                for existing_id in current_group:
                    current_group[existing_id] -= position
                # Add new reference
                current_group[i] = 0
                grouped_ids[i] = current_group

                ref_id = i
                seq_to_id[seq] = i
                del seq_to_id[matched_seq]
                assert len(seq) > len(matched_seq), (
                    "New reference sequence should be longer than previous matched sequence"
                )
                print(
                    f"Updated reference sequence L:{len(matched_seq)} for group {ref_id} to new sequence L:{len(seq)} (superstring case)"
                    f"\nPrevious sequence: {matched_seq}\n     New sequence: {seq}"
                )
            elif position >= 0 and position + len(seq) > len(matched_seq):
                #    New sequence is partial substring but exceeds matched_seq length,
                #    Keep current reference and merge sequences
                #    ....
                #      ....
                grouped_ids[ref_id][i] = position
                new_ref_seq = ["-"] * (position + len(seq))

                new_ref_seq[0 : len(matched_seq)] = list(matched_seq)
                # This is exact match, does not matter which sequence we take
                new_ref_seq[position : position + len(seq)] = list(seq)
                new_ref_seq_str = "".join(new_ref_seq)
                assert "-" not in new_ref_seq_str, "Merged reference sequence has gaps"
                # Update the reference sequence in seq_to_id
                seq_to_id[new_ref_seq_str] = ref_id
                del seq_to_id[matched_seq]
                assert len(new_ref_seq_str) > len(matched_seq), (
                    "New reference sequence should be longer than previous matched sequence"
                )
                print(
                    f"Updated reference sequence L:{len(matched_seq)} for group {ref_id} to new sequence L:{len(new_ref_seq_str)} (substring extending end)"
                    f"\nPrevious sequence: {matched_seq}\n     New sequence: {new_ref_seq_str}"
                )

            elif position < 0 and position + len(seq) < len(matched_seq):
                #    New sequence is partial substring but shorter than matched_seq length, merge sequences
                #      ....
                #    ....
                current_group = grouped_ids.pop(ref_id)
                # Shift all existing positions by -position
                for existing_id in current_group:
                    current_group[existing_id] -= position
                # Add new reference
                current_group[i] = 0
                grouped_ids[i] = current_group
                new_ref_seq = ["-"] * (-position + len(matched_seq))

                new_ref_seq[-position:] = list(matched_seq)
                # This is exact match, does not matter which sequence we take
                new_ref_seq[0 : len(seq)] = list(seq)
                new_ref_seq_str = "".join(new_ref_seq)
                assert "-" not in new_ref_seq_str, "Merged reference sequence has gaps"
                # Update the reference sequence in seq_to_id
                seq_to_id[new_ref_seq_str] = ref_id
                del seq_to_id[matched_seq]
                assert len(new_ref_seq_str) > len(matched_seq), (
                    "New reference sequence should be longer than previous matched sequence"
                )
                print(
                    f"Updated reference sequence L:{len(matched_seq)} for group {ref_id} to new sequence L:{len(new_ref_seq_str)} (substring extending start)"
                    f"\nPrevious sequence: {matched_seq}\n     New sequence: {new_ref_seq_str}"
                )
            else:
                assert False, "Unhandled case in partial sequence matching"
    return grouped_ids, seq_to_id

def groupby_sequence(data, min_match_ratio=0.9):
    grouped_ids, seq_to_id = collect_grouped_idx(
        data["sequences"], min_match_ratio=min_match_ratio
    )
    # Now we have grouped_ids with groupings, need pad the cordinates to match the longest sequence in each group
    grouped_data = {}
    for ref_id, entries in grouped_ids.items():
        grouped_data[ref_id] = {
            "sequence": [],
            "data_cif_files": [],
            "pdb_ids": [],
            "xyz": [],
        }
        ref_seq = next((k for k, v in seq_to_id.items() if v == ref_id), None)
        ref_len = len(ref_seq)
        print("Reconstructing group for ref_id:", ref_id, "ref_seq:", ref_seq)
        ref_pdb_id = data["data_cif_files"][ref_id].split("/")[-1].split(".")[0]

        for entry_idx, position in entries.items():
            seq = data["sequence"][entry_idx]
            xyz = data["xyz"][entry_idx]
            cif_file = data["data_cif_files"][entry_idx]
            pdb_id = cif_file.split("/")[-1].split(".")[0]
            L = len(seq)
            msg = "The reference sequence should always be the longest one"
            assert L <= ref_len, get_message(**locals())
            #     "The reference sequence should always be the longest one (reference length: "
            #     msg,
            #     f"{ref_len}, current sequence length: {L})"
            #     f" (ref_id: {ref_id} ({ref_pdb_id}), entry_idx: {entry_idx} ({pdb_id}))"
            #     f" (position: {position})"
            #     f" (sequence: {seq})"
            #     f" (reference sequence: {data['sequence'][ref_id]})"
            # )
            assert position >= 0, "The position should be non-negative"
            msg = "The sequence plus position should not exceed reference length"
            assert position + L <= ref_len, get_message(**locals())
            if len(seq) < ref_len:
                # Get reference sequence
                # Pad the current sequence with gaps '-' to match the reference
                new_seq = ["-"] * ref_len
                new_seq[position : position + L] = list(seq)
                new_seq_str = "".join(new_seq)
                print
                print("Previous sequence:", seq, "\npadded sequence:  ", new_seq_str)

                assert len(new_seq_str) == ref_len, get_message(**locals())
                msg = "Unexpected length of sequence after padding"
                assert len(new_seq_str.replace("-", "")) == L, get_message(**locals())
                seq = new_seq_str

                # Now pad the coordinates
                ref_xyz = data["xyz"][ref_id]
                # Dict should be the same for all residues in the reference
                atom_dict = ref_xyz[0]
                new_xyz = [
                    {
                        group: np.full_like(atom_dict[group], np.nan)
                        for group in atom_dict
                    }
                    for _ in range(ref_len)
                ]
                new_xyz[position : position + L] = xyz
                msg = "Unexpected length of coordinates after padding"
                assert len(new_xyz) == ref_len, get_message(**locals())
                xyz = new_xyz

            grouped_data[ref_id]["sequence"].append(seq)
            grouped_data[ref_id]["data_cif_files"].append(cif_file)
            grouped_data[ref_id]["pdb_ids"].append(pdb_id)
            grouped_data[ref_id]["xyz"].append(xyz)
        grouped_data[ref_id]["ref_sequence"] = ref_seq

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
    parser.add_argument(
        "--min_match_ratio",
        type=float,
        default=1,
        help="Minimum match ratio for partial sequence matching",
    )
    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    min_match_ratio = args.min_match_ratio
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_file, "rb") as f:
        data = pickle.load(f)

    grouped_data = groupby_sequence(data, min_match_ratio=min_match_ratio)

    for group_id, entries in grouped_data.items():
        print(f"Group {group_id} has {len(entries['pdb_ids'])} entries.")
        # data[group_id]['pdb_id']
        cif_file = data["data_cif_files"][group_id]
        group_pdb_id = cif_file.split("/")[-1].split(".")[0]
        output_file = output_dir / f"{group_pdb_id}_grouped.pkl"
        with open(output_file, "wb") as f:
            pickle.dump(entries, f)
        # Save the sequences as a FASTA file
        fasta_output_path = output_dir / f"{group_pdb_id}.fasta"
        with open(fasta_output_path, "w") as fasta_file:
            for seq, pdb_id in zip(entries["sequence"], entries["pdb_ids"]):
                fasta_file.write(f">{pdb_id}\n")
                fasta_file.write(f"{seq}\n")
            fasta_file.write(f">reference\n")
            fasta_file.write(f"{entries['ref_sequence']}\n")
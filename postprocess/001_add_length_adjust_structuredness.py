# Corrects the structuredness calculated before 3.2.0rc5 to take into account the unresolved residues
# The structuredness is calculated as a fraction of resolved/observed residues and the adjusted version is calculated as a fraction of all expected (observed/unobserved) residues
import argparse
import pandas as pd


def main(input_file: str, output_file: str):
    # Read CSV
    df = pd.read_csv(input_file, keep_default_na=False, na_values=[""])
    modified = False
    added_fields = []

    # Track missing source columns to warn the user
    missing_sources = []

    # Add diagnostics / length columns if missing
    if "length" not in df.columns:
        if "sequence" in df.columns:
            df["length"] = df["sequence"].str.len()
            modified = True
            added_fields.append("length")
        else:
            missing_sources.append("sequence")
    if "length_observed" not in df.columns:
        if "sequence_observed" in df.columns:
            seq_observed = df["sequence_observed"]
            # The observed sequence has '?' as unobserved residues, removing them to count only observed residues
            length_observed = seq_observed.str.split(":").apply(
                lambda x: len([res for res in x if res != "?"])
            )
            df["length_observed"] = length_observed
            modified = True
            added_fields.append("length_observed")
        else:
            missing_sources.append("sequence_observed")
            length_observed = None
    else:
        length_observed = df["length_observed"]

    if "length_expected" not in df.columns:
        if "sequence_expected" in df.columns:
            seq_expected = df["sequence_expected"]
            length_expected = seq_expected.str.split(":").str.len()
            df["length_expected"] = length_expected
            modified = True
            added_fields.append("length_expected")
        else:
            missing_sources.append("sequence_expected")
            length_expected = None
    else:
        length_expected = df["length_expected"]

    if "fraction_observed" not in df.columns:
        if (length_observed is not None) and (length_expected is not None):
            fraction_observed = length_observed / length_expected
            df["fraction_observed"] = fraction_observed
            modified = True
            added_fields.append("fraction_observed")
        else:
            missing_sources.append("length_observed/length_expected")
            fraction_observed = None
    else:
        fraction_observed = df["fraction_observed"]

    # Add adjusted structuredness columns only if they don't already exist
    if "inter_chain_structuredness_adjusted" not in df.columns:
        if "inter_chain_structuredness" in df.columns:
            if fraction_observed is not None:
                df["inter_chain_structuredness_adjusted"] = (
                    df["inter_chain_structuredness"] * fraction_observed
                )
                modified = True
                added_fields.append("inter_chain_structuredness_adjusted")
            else:
                print(
                    "Warning: cannot compute inter_chain_structuredness_adjusted because fraction_observed is missing"
                )
        else:
            print(
                "Warning: source column 'inter_chain_structuredness' not found; skipping adjusted field"
            )
    if "intra_chain_structuredness_adjusted" not in df.columns:
        if "intra_chain_structuredness" in df.columns:
            if fraction_observed is not None:
                df["intra_chain_structuredness_adjusted"] = (
                    df["intra_chain_structuredness"] * fraction_observed
                )
                modified = True
                added_fields.append("intra_chain_structuredness_adjusted")
            else:
                print(
                    "Warning: cannot compute intra_chain_structuredness_adjusted because fraction_observed is missing"
                )
        else:
            print(
                "Warning: source column 'intra_chain_structuredness' not found; skipping adjusted field"
            )

    if "total_structuredness_adjusted" not in df.columns:
        if "total_structuredness" in df.columns:
            if fraction_observed is not None:
                df["total_structuredness_adjusted"] = (
                    df["total_structuredness"] * fraction_observed
                )
                modified = True
                added_fields.append("total_structuredness_adjusted")
            else:
                print(
                    "Warning: cannot compute total_structuredness_adjusted because fraction_observed is missing"
                )
        else:
            print(
                "Warning: source column 'total_structuredness' not found; skipping adjusted field"
            )

    # Report what changed and warnings
    if len(missing_sources) > 0:
        print(
            "Warning: missing source columns or values:",
            ", ".join(sorted(set(missing_sources))),
        )

    if modified:
        print("Added/modified columns:", ", ".join(added_fields))
        df.to_csv(output_file, float_format="%.3f", index=False)
        print(f"Wrote updated file to: {output_file}")
    else:
        print("No changes made; nothing written.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Adjust structuredness values to account for unresolved residues"
    )
    parser.add_argument(
        "input_file",
        help="Input metadata CSV with structuredness",
    )
    parser.add_argument(
        "output_file",
        help="Output CSV path for adjusted structuredness",
    )
    args = parser.parse_args()
    main(args.input_file, args.output_file)

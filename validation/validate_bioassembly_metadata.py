#!/usr/bin/env python3
"""
Validate bioassembly metadata CSV file.

This script validates that:
1. The sequence field is properly reconstructed using stoichiometry and all_sequences
2. Each chain's copies in stoichiometry matches the number of chains in FASTA description
3. The sequence is a correct concatenation: chain_id_sequence1*copies1 + chain_id_sequence2*copies2 + ...
"""

import csv
from functools import partial
from pathlib import Path
import re
import sys
from typing import Dict, List, Tuple
from dataclasses import dataclass
import traceback

try:
    from biotite import sequence as bioseq

    has_biotite = True
except ImportError:
    has_biotite = False
# Import parse_fasta from the tools/fasta directory
sys.path.insert(0, str(Path(__file__).parent / "tools/fasta"))
from chain_parser import parse_fasta

csv.field_size_limit(sys.maxsize)


@dataclass
class ValidationResult:
    """Result of validation for a single row."""

    pdb_id: str
    is_valid: bool
    errors: List[str]
    warnings: List[str] = None


def parse_stoichiometry(stoichiometry: str) -> List[Tuple[str, int]]:
    """
    Parse stoichiometry string into list of (chain_id, copies) tuples.

    Args:
        stoichiometry: String like "A:2" or "A:1;B:1" or "B:60;C:60"

    Returns:
        List of tuples: [("A", 2)] or [("A", 1), ("B", 1)]
    """
    if not stoichiometry:
        return []

    result = []
    for part in stoichiometry.split(";"):
        chain_id, copies = part.split(":")
        result.append((chain_id.strip(), int(copies.strip())))

    return result


def validate_row(
    pdb_id: str, stoichiometry: str, sequence: str, all_sequences: str
) -> ValidationResult:
    """
    Validate a single row of the CSV.

    Args:
        pdb_id: PDB ID
        stoichiometry: Stoichiometry string
        sequence: Expected concatenated sequence
        all_sequences: FASTA formatted sequences

    Returns:
        ValidationResult object
    """
    errors = []
    warnings = []
    # Parse stoichiometry
    try:
        stoich_list = parse_stoichiometry(stoichiometry)
    except Exception as e:
        errors.append(f"Failed to parse stoichiometry: {e}\n{traceback.format_exc()}")
        return ValidationResult(pdb_id, False, errors)

    # Parse FASTA
    try:
        fasta_dict = parse_fasta(all_sequences)
    except Exception as e:
        errors.append(f"Failed to parse FASTA: {e}\n{traceback.format_exc()}")
        return ValidationResult(pdb_id, False, errors)

    # Validate each chain in stoichiometry
    reconstructed_sequence = []

    for chain_id, copies in stoich_list:
        # Check if chain exists in FASTA
        if chain_id not in fasta_dict:
            errors.append(f"Chain '{chain_id}' from stoichiometry not found in FASTA")
            continue

        chain_sequence, chain_list = fasta_dict[chain_id]

        # Validate that copies matches the number of chains in FASTA description
        expected_copies = len(chain_list)
        if copies != expected_copies:
            warnings.append(
                f"Chain '{chain_id}': stoichiometry says {copies} copies, "
                f"but FASTA description lists {expected_copies} chains: {', '.join(chain_list)}"
            )

        # Add to reconstructed sequence
        reconstructed_sequence.append(chain_sequence * copies)

    # Reconstruct the full sequence
    reconstructed = "".join(reconstructed_sequence)

    # Compare with provided sequence
    if reconstructed != sequence:
        # Allow mismatch with X in reconstructed sequence
        matches = [i == j for i, j in zip(sequence, reconstructed) if j != "X"]
        if all(matches) and len(sequence) == len(reconstructed):
            pass  # Acceptable mismatch due to 'X'
        else:
            # Calculate alignment using biotite
            # RNASequence = partial(bioseq.GeneralSequence, bioseq.Alphabet("ACGU"))
            # seq1 = RNASequence(sequence)
            # seq2 = RNASequence(reconstructed)
            if has_biotite:
                try:
                    sequence_t = sequence.replace("U", "T")
                    reconstructed_t = reconstructed.replace("U", "T").replace("X", "N")
                    seq1 = bioseq.NucleotideSequence(sequence_t)
                    seq2 = bioseq.NucleotideSequence(reconstructed_t)
                    alignments = bioseq.align.align_optimal(
                        seq1,
                        seq2,
                        matrix=bioseq.align.SubstitutionMatrix.std_nucleotide_matrix(),
                        gap_penalty=-5,
                        terminal_penalty=False,
                    )
                    best_alignment = str(alignments[0]).replace("T", "U")
                except Exception as e:
                    best_alignment = f"Error computing alignment: {e}"
            else:
                best_alignment = "Biotite library not installed; alignment unavailable."
            errors.append(
                f"Sequence mismatch:\n"
                f"  Expected length: {len(sequence)}\n"
                f"  Reconstructed length: {len(reconstructed)}\n"
                f"  Expected:      {sequence}\n"
                f"  Reconstructed: {reconstructed}\n"
                f"  Is partial match: {sequence in reconstructed or reconstructed in sequence}\n"
                f"  Alignment:\n{best_alignment}\n"
            )

    is_valid = len(errors) == 0
    return ValidationResult(pdb_id, is_valid, errors, warnings)


def validate_csv_file(
    filepath: str, verbose: bool = False
) -> Tuple[int, int, List[ValidationResult]]:
    """
    Validate entire CSV file.

    Args:
        filepath: Path to CSV file
        verbose: If True, print details for all rows; if False, only print errors

    Returns:
        Tuple of (total_rows, valid_rows, list_of_failed_validations)
    """
    total_rows = 0
    valid_rows = 0
    failed_validations = []

    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            total_rows += 1
            if "target_id" not in row:
                target_id = row["pdb_id"]
            else:
                target_id = row["target_id"]

            stoichiometry = row["stoichiometry"]
            sequence = row["sequence"]
            all_sequences = row["all_sequences"]

            result = validate_row(target_id, stoichiometry, sequence, all_sequences)

            if result.is_valid:
                valid_rows += 1
                if verbose:
                    print(f"✓ {target_id}: VALID")
                    for warning in result.warnings:
                        print(f"  - WARNING: {warning}")
                    print()
            else:
                failed_validations.append(result)
                print(f"✗ {target_id}: INVALID")
                for error in result.errors:
                    print(f"  - {error}")
                for warning in result.warnings:
                    print(f"  - WARNING: {warning}")
                print()

    return total_rows, valid_rows, failed_validations


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate bioassembly metadata CSV file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s bioassembly_metadata.csv
  %(prog)s bioassembly_metadata.csv --verbose
  %(prog)s bioassembly_metadata.csv --summary-only
        """,
    )
    parser.add_argument("csv_file", help="Path to CSV file to validate")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print details for all rows (not just errors)",
    )
    parser.add_argument(
        "-s",
        "--summary-only",
        action="store_true",
        help="Only print summary statistics",
    )

    args = parser.parse_args()

    print(f"Validating {args.csv_file}...")
    print()

    # Validate
    if args.summary_only:
        # Suppress individual error output
        import io
        from contextlib import redirect_stdout

        with redirect_stdout(io.StringIO()):
            total, valid, failed = validate_csv_file(args.csv_file, args.verbose)
    else:
        total, valid, failed = validate_csv_file(args.csv_file, args.verbose)

    # Print summary
    print("=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"Total rows validated: {total}")
    print(f"Valid rows: {valid}")
    print(f"Invalid rows: {len(failed)}")

    if total > 0:
        success_rate = (valid / total) * 100
        print(f"Success rate: {success_rate:.2f}%")

    if failed:
        print()
        print(f"Failed PDB IDs: {', '.join([r.pdb_id for r in failed])}")

    # Exit with appropriate code
    sys.exit(0 if len(failed) == 0 else 1)


if __name__ == "__main__":
    main()

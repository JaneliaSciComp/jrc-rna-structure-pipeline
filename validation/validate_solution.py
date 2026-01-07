#!/usr/bin/env python3
"""
Validate solution CSV file against metadata CSV file.

This script validates that:
1. Each row in the solutions file has an ID field matching pattern {target_id}_{sequence_idx}
2. The resname in the solution matches the sequence_idx letter from the sequence field in metadata
3. All target_ids in solutions exist in the metadata file
4. The sequence_idx (resid) is within the bounds of the sequence length
"""

import csv
import sys
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
import traceback

csv.field_size_limit(sys.maxsize)


@dataclass
class ValidationResult:
    """Result of validation for solutions."""

    target_id: str
    is_valid: bool
    errors: List[str]
    warnings: List[str] = None
    total_residues: int = 0
    validated_residues: int = 0


def load_metadata(metadata_file: str) -> Dict[str, str]:
    """
    Load metadata file and extract target_id -> sequence mapping.

    Args:
        metadata_file: Path to metadata CSV file

    Returns:
        Dictionary mapping target_id to sequence string
    """
    metadata = {}

    with open(metadata_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            target_id = row["target_id"]
            sequence = row["sequence"]
            metadata[target_id] = sequence

    return metadata


def parse_solution_id(solution_id: str) -> Tuple[str, int]:
    """
    Parse solution ID into target_id and sequence_idx.

    Args:
        solution_id: ID string like "17RA_A_1" or "1A1T_A_15"

    Returns:
        Tuple of (target_id, sequence_idx)
        Example: ("17RA_A", 1) or ("1A1T_A", 15)
    """
    # Split by underscore and take the last part as sequence_idx
    parts = solution_id.rsplit("_", 1)
    if len(parts) != 2:
        raise ValueError(f"Invalid ID format: {solution_id}")

    target_id = parts[0]
    try:
        sequence_idx = int(parts[1])
    except ValueError:
        raise ValueError(f"Invalid sequence index in ID: {solution_id}")

    return target_id, sequence_idx


def validate_solutions(
    solutions_file: str, metadata: Dict[str, str], verbose: bool = False
) -> Tuple[Dict[str, ValidationResult], int, int]:
    """
    Validate solutions file against metadata.

    Args:
        solutions_file: Path to solutions CSV file
        metadata: Dictionary mapping target_id to sequence
        verbose: If True, print details for all rows

    Returns:
        Tuple of (validation_results_by_target, total_rows, valid_rows)
    """
    # Track validation results per target_id
    results_by_target = defaultdict(
        lambda: {
            "errors": [],
            "warnings": [],
            "total_residues": 0,
            "validated_residues": 0,
            "seen_positions": set(),  # Track which positions we've seen
        }
    )

    # Track which target_ids we've seen
    seen_targets = set()

    total_rows = 0
    valid_rows = 0

    with open(solutions_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            total_rows += 1
            solution_id = row["ID"]
            resname = row["resname"]
            resid = int(row["resid"])

            try:
                # Parse the solution ID
                target_id, sequence_idx = parse_solution_id(solution_id)
                seen_targets.add(target_id)

                # Get the target's data
                target_data = results_by_target[target_id]
                target_data["total_residues"] += 1

                # Check if target exists in metadata
                if target_id not in metadata:
                    if (
                        f"Target '{target_id}' not found in metadata"
                        not in target_data["errors"]
                    ):
                        target_data["errors"].append(
                            f"Target '{target_id}' not found in metadata"
                        )
                    continue

                sequence = metadata[target_id]

                # Validate sequence_idx is within bounds (1-indexed)
                if sequence_idx < 1 or sequence_idx > len(sequence):
                    target_data["errors"].append(
                        f"Row {solution_id}: sequence_idx {sequence_idx} out of bounds "
                        f"(sequence length: {len(sequence)})"
                    )
                    continue

                # Validate that resid matches sequence_idx
                if resid != sequence_idx:
                    target_data["warnings"].append(
                        f"Row {solution_id}: resid ({resid}) does not match "
                        f"sequence_idx ({sequence_idx})"
                    )

                # Validate that resname matches the sequence character
                expected_resname = sequence[sequence_idx - 1]  # Convert to 0-indexed
                if resname != expected_resname:
                    target_data["errors"].append(
                        f"Row {solution_id}: resname mismatch - expected '{expected_resname}' "
                        f"at position {sequence_idx}, got '{resname}'"
                    )
                    continue

                # Track this position as seen
                target_data["seen_positions"].add(sequence_idx)

                # This row is valid
                target_data["validated_residues"] += 1
                valid_rows += 1

            except ValueError as e:
                # Error parsing the ID
                if "Invalid ID format" in str(e) or "Invalid sequence index" in str(e):
                    # We don't know the target_id, so create a generic error entry
                    error_key = f"_parse_error_{solution_id}"
                    if error_key not in results_by_target:
                        results_by_target[error_key]["errors"].append(
                            f"Failed to parse solution ID '{solution_id}': {e}"
                        )
                        results_by_target[error_key]["total_residues"] = 1
                else:
                    raise
            except Exception as e:
                print(f"Unexpected error processing row {solution_id}: {e}")
                print(traceback.format_exc())

    # Convert to ValidationResult objects and check for missing positions
    validation_results = {}
    for target_id, data in results_by_target.items():
        # Check if any positions are missing
        if target_id in metadata:
            sequence = metadata[target_id]
            expected_positions = set(range(1, len(sequence) + 1))
            seen_positions = data["seen_positions"]
            missing_positions = expected_positions - seen_positions

            if missing_positions:
                # Sort and format missing positions
                missing_sorted = sorted(missing_positions)
                if len(missing_sorted) <= 10:
                    missing_str = ", ".join(map(str, missing_sorted))
                else:
                    missing_str = f"{', '.join(map(str, missing_sorted[:10]))} ... and {len(missing_sorted) - 10} more"

                data["errors"].append(
                    f"Missing {len(missing_positions)} position(s) out of {len(sequence)}: {missing_str}"
                )

        is_valid = len(data["errors"]) == 0
        validation_results[target_id] = ValidationResult(
            target_id=target_id,
            is_valid=is_valid,
            errors=data["errors"],
            warnings=data["warnings"],
            total_residues=data["total_residues"],
            validated_residues=data["validated_residues"],
        )

    return validation_results, total_rows, valid_rows


def check_metadata_coverage(
    metadata: Dict[str, str], seen_targets: Set[str]
) -> List[str]:
    """
    Check if all targets in metadata have corresponding solutions.

    Args:
        metadata: Dictionary of target_id -> sequence
        seen_targets: Set of target_ids found in solutions file

    Returns:
        List of target_ids in metadata but not in solutions
    """
    metadata_targets = set(metadata.keys())
    missing_in_solutions = metadata_targets - seen_targets
    return sorted(list(missing_in_solutions))


def check_solution_coverage(
    metadata: Dict[str, str], seen_targets: Set[str]
) -> List[str]:
    """
    Check if all targets in solution have corresponding metadata.

    Args:
        metadata: Dictionary of target_id -> sequence
        seen_targets: Set of target_ids found in solutions file

    Returns:
        List of target_ids in metadata but not in solutions
    """
    metadata_targets = set(metadata.keys())
    missing_in_metadata = seen_targets - metadata_targets
    return sorted(list(missing_in_metadata))


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate solution CSV file against metadata CSV file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s solutions.csv metadata.csv
  %(prog)s solutions.csv metadata.csv --verbose
  %(prog)s solutions.csv metadata.csv --summary-only
        """,
    )
    parser.add_argument("solutions_file", help="Path to solutions CSV file")
    parser.add_argument("metadata_file", help="Path to metadata CSV file")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print details for all targets (not just errors)",
    )
    parser.add_argument(
        "-s",
        "--summary-only",
        action="store_true",
        help="Only print summary statistics",
    )
    parser.add_argument(
        "--check-coverage",
        action="store_true",
        help="Check if all metadata targets have solutions",
    )

    args = parser.parse_args()

    print(f"Loading metadata from {args.metadata_file}...")
    metadata = load_metadata(args.metadata_file)
    print(f"Loaded {len(metadata)} targets from metadata")
    print()

    print(f"Validating {args.solutions_file}...")
    print()

    # Validate solutions
    if args.summary_only:
        # Suppress individual error output
        import io
        from contextlib import redirect_stdout

        with redirect_stdout(io.StringIO()):
            validation_results, total_rows, valid_rows = validate_solutions(
                args.solutions_file, metadata, args.verbose
            )
    else:
        validation_results, total_rows, valid_rows = validate_solutions(
            args.solutions_file, metadata, args.verbose
        )

    # Print results per target
    if not args.summary_only:
        for target_id in sorted(validation_results.keys()):
            result = validation_results[target_id]

            if result.is_valid:
                if args.verbose:
                    print(f"✓ {target_id}: VALID")
                    print(
                        f"  Validated {result.validated_residues}/{result.total_residues} residues"
                    )
                    if result.warnings:
                        for warning in result.warnings:
                            print(f"  - WARNING: {warning}")
                    print()
            else:
                print(f"✗ {target_id}: INVALID")
                print(
                    f"  Validated {result.validated_residues}/{result.total_residues} residues"
                )
                for error in result.errors:
                    print(f"  - {error}")
                if result.warnings:
                    for warning in result.warnings:
                        print(f"  - WARNING: {warning}")
                print()

    # Count valid/invalid targets
    valid_targets = sum(1 for r in validation_results.values() if r.is_valid)
    invalid_targets = len(validation_results) - valid_targets

    # Check coverage if requested
    seen_targets = set(validation_results.keys())
    missing_targets = []
    if args.check_coverage:
        missing_targets = check_metadata_coverage(metadata, seen_targets)
        if missing_targets:
            print("=" * 70)
            print("METADATA COVERAGE CHECK")
            print("=" * 70)
            print(
                f"Found {len(missing_targets)} targets in metadata without solutions:"
            )
            for target in missing_targets[:20]:  # Show first 20
                print(f"  - {target}")
            if len(missing_targets) > 20:
                print(f"  ... and {len(missing_targets) - 20} more")
            print()
        missing_targets_solution = check_solution_coverage(metadata, seen_targets)
        if missing_targets_solution:
            print("=" * 70)
            print("SOLUTION COVERAGE CHECK")
            print("=" * 70)
            print(
                f"Found {len(missing_targets_solution)} targets in solution without metadata:"
            )
            for target in missing_targets_solution[:20]:  # Show first 20
                print(f"  - {target}")
            if len(missing_targets_solution) > 20:
                print(f"  ... and {len(missing_targets_solution) - 20} more")
            print()

    # Print summary
    print("=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    print(f"Total solution rows validated: {total_rows}")
    print(f"Valid rows: {valid_rows}")
    print(f"Invalid rows: {total_rows - valid_rows}")

    if total_rows > 0:
        success_rate = (valid_rows / total_rows) * 100
        print(f"Row success rate: {success_rate:.2f}%")

    print()
    print(f"Total targets in solutions: {len(validation_results)}")
    print(f"Valid targets: {valid_targets}")
    print(f"Invalid targets: {invalid_targets}")

    if len(validation_results) > 0:
        target_success_rate = (valid_targets / len(validation_results)) * 100
        print(f"Target success rate: {target_success_rate:.2f}%")

    if args.check_coverage:
        print()
        print(f"Targets in metadata: {len(metadata)}")
        print(f"Targets with solutions: {len(seen_targets)}")
        print(f"Targets without solutions: {len(missing_targets)}")
        print(f"Targets in solution without metadata: {len(missing_targets_solution)}")

    if invalid_targets > 0:
        print()
        invalid_target_ids = [
            r.target_id for r in validation_results.values() if not r.is_valid
        ]
        print(f"Invalid target IDs: {', '.join(sorted(invalid_target_ids)[:10])}")
        if len(invalid_target_ids) > 10:
            print(f"  ... and {len(invalid_target_ids) - 10} more")

    # Exit with appropriate code
    sys.exit(0 if invalid_targets == 0 else 1)


if __name__ == "__main__":
    main()

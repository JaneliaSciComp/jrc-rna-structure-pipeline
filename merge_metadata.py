import duckdb
import sys
from pathlib import Path
import argparse
import glob


def get_read_function(file_pattern: str) -> str:
    """Determine the appropriate DuckDB read function based on file extension."""
    # Get the first file to determine extension
    files = glob.glob(file_pattern)
    if not files:
        # If no glob match, use the pattern itself
        sample_file = file_pattern
    else:
        sample_file = files[0]

    ext = Path(sample_file).suffix.lower()
    if ext == ".json":
        return "read_json_auto", ""
    elif ext == ".csv":
        return "read_csv_auto", ", nullstr=nan"
    else:
        raise ValueError(f"Unsupported file extension: {ext}. Use .json or .csv")


def merge_sequences_and_metadata(
    sequences_file: str,
    metadata_file: str,
    output_file: str = None,
    exclude_columns: list = ["mmcif_file"],
    join_columns: list = ["pdb_id"],
    extract_pdb_id: bool = True,
):
    """
    Merge sequences with metadata using DuckDB.

    Args:
        sequences_file: Path or glob pattern to CSV/JSON file(s) with sequences
        metadata_file: Path or glob pattern to JSON/CSV file(s) with metadata
        output_file: Optional output file path (default: stdout)
        exclude_columns: Columns to exclude from the output
        join_columns: Columns to use for join condition
        extract_pdb_id: Whether to extract pdb_id from target_id if not present
    """
    try:
        # Initialize DuckDB connection
        conn = duckdb.connect(":memory:")

        # Determine read functions based on file extensions
        seq_read_func, seq_options = get_read_function(sequences_file)
        meta_read_func, meta_options = get_read_function(metadata_file)

        # Load sequences table
        conn.execute(f"""
            CREATE TABLE sequences AS
            SELECT * FROM {seq_read_func}('{sequences_file}' {seq_options})
        """)

        # Check if pdb_id column exists and add it if needed
        columns = conn.execute("SELECT * FROM sequences LIMIT 0").description
        column_names = [col[0] for col in columns]

        if (
            "pdb_id" not in column_names
            and extract_pdb_id
            and "target_id" in column_names
        ):
            conn.execute("""
                ALTER TABLE sequences ADD COLUMN pdb_id VARCHAR
            """)
            conn.execute("""
                UPDATE sequences SET pdb_id = substring(target_id, 1, 4)
            """)

        # Load metadata table
        conn.execute(f"""
            CREATE TABLE metadata AS
            SELECT * FROM {meta_read_func}('{metadata_file}' {meta_options})
        """)

        # Build exclude clause
        exclude_clause = ""
        if exclude_columns:
            exclude_clause = f"EXCLUDE ({', '.join(exclude_columns)})"

        # Build join condition
        join_condition = ", ".join(join_columns)

        # Merge and save
        conn.execute(f"""
            COPY (
                SELECT 
                    * {exclude_clause}
                FROM sequences
                LEFT JOIN metadata USING ({join_condition})
                ORDER BY temporal_cutoff, target_id 
            ) TO '{output_file}' (FORMAT CSV, HEADER)
        """)
        print(f"Results saved to {output_file}")

        conn.close()

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge sequences with metadata using DuckDB. Supports CSV and JSON files with glob patterns."
    )
    parser.add_argument(
        "sequences_file", help="Path or glob pattern to CSV/JSON file(s) with sequences"
    )
    parser.add_argument(
        "metadata_file", help="Path or glob pattern to JSON/CSV file(s) with metadata"
    )
    parser.add_argument("output_file", help="Output file path for merged data")
    parser.add_argument(
        "--exclude-columns",
        nargs="+",
        default=["mmcif_file"],
        help="Columns to exclude from the output",
    )
    parser.add_argument(
        "--join-columns",
        nargs="+",
        default=["pdb_id"],
        help="Columns to use for join condition (default: pdb_id)",
    )
    parser.add_argument(
        "--no-extract-pdb-id",
        action="store_true",
        help="Disable automatic extraction of pdb_id from target_id",
    )

    args = parser.parse_args()

    sequences_file = args.sequences_file
    metadata_file = args.metadata_file
    output_file = args.output_file

    # Check if input files exist (only for non-glob patterns)
    if "*" not in sequences_file and "?" not in sequences_file:
        if not Path(sequences_file).exists():
            print(
                f"Error: Sequences file '{sequences_file}' does not exist.",
                file=sys.stderr,
            )
            sys.exit(1)
    if "*" not in metadata_file and "?" not in metadata_file:
        if not Path(metadata_file).exists():
            print(
                f"Error: Metadata file '{metadata_file}' does not exist.",
                file=sys.stderr,
            )
            sys.exit(1)

    merge_sequences_and_metadata(
        sequences_file,
        metadata_file,
        output_file,
        exclude_columns=args.exclude_columns,
        join_columns=args.join_columns,
        extract_pdb_id=not args.no_extract_pdb_id,
    )

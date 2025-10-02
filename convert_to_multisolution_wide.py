import pandas as pd


def to_wide(result_solution, select=None):
    """
    Convert long format solution dataframe to wide format and save to CSV, keeping only specified atoms.

    Keeps all atoms if select is None, otherwise filters to atoms starting with the given prefixes.
    """
    # Columns to keep that must be present
    keep = ["group_id", "resname", "resid", "solution_id"]
    columns = set(result_solution.columns)
    missing = set(keep) - columns
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Create selected only version
    coordinate_mask = result_solution.columns.str.endswith(("_x", "_y", "_z"))
    if select is None:
        selected_columns = result_solution.columns[coordinate_mask]
    else:
        selected_columns = result_solution.columns[
            result_solution.columns.str.startswith(tuple(select)) & coordinate_mask
        ]
    if len(selected_columns) == 0:
        raise ValueError(f"No columns found for selected atoms: {select}")

    # # Non-coordinate columns
    # other_columns = list(result_solution.columns[~coordinate_mask])
    # result_solution_selected = result_solution[other_columns + selected_columns]

    # Use atom coordinates and solution_id as multiindex for columns
    result_solution_selected = (
        result_solution[keep + list(selected_columns)]
        .set_index(keep)
        .unstack("solution_id")
    )

    # Map multilevel column names to {x,y,z}_index
    result_solution_selected.columns = [
        f"{i}_{j}" if j != "" else f"{i}" for i, j in result_solution_selected.columns
    ]

    # Sort by resid within each group
    result_solution_selected = result_solution_selected.sort_values(
        by=["group_id", "resid"]
    )

    # Reconstruct ID
    result_solution_selected = result_solution_selected.reset_index()
    result_solution_selected["ID"] = (
        result_solution_selected["group_id"].astype(str)
        + "_"
        + (result_solution_selected["resid"]).astype(str)
    )

    # Remove C1'_ prefix from column names
    result_solution_selected.columns = [
        col.replace("C1'_", "") for col in result_solution_selected.columns
    ]

    # Move ID to front
    cols = list(result_solution_selected.columns)
    cols.remove("ID")
    cols = ["ID"] + cols
    result_solution_selected = result_solution_selected[cols]
    return result_solution_selected


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file", type=str, help="Path to the input CSV file in long format"
    )
    parser.add_argument("output_file", type=str, help="Output CSV file in wide format")

    parser.add_argument(
        "--select",
        nargs="+",
        type=str,
        default=None,
        help="Select atoms starting with this prefix",
    )

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    result_solution = pd.read_csv(input_file)
    result_solution_selected = to_wide(result_solution, select=args.select)

    result_solution_selected.to_csv(
        output_file,
        float_format="%.3f",
        index=False,
    )

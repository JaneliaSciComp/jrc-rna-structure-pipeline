def to_wide(result_solution):
    """
    Convert long format solution dataframe to wide format and save to CSV.
    """
    import pandas as pd

    # Create C1' only version
    c1_columns = [col for col in result_solution.columns if "C1'" in col]

    # Non-coordinate columns
    other_columns = list(
        result_solution.columns[
            ~result_solution.columns.str.endswith(("_x", "_y", "_z"))
        ]
    )
    result_solution_c1 = result_solution[other_columns + c1_columns]

    # Use atom coordinates and solution_id as multiindex for columns
    keep = ["group_id", "resname", "resid", "solution_id"]
    result_solution_c1 = (
        result_solution_c1[keep + c1_columns].set_index(keep).unstack("solution_id")
    )

    # Map multilevel column names to {x,y,z}_index
    result_solution_c1.columns = [
        f"{i}_{j}" if j != "" else f"{i}" for i, j in result_solution_c1.columns
    ]

    # Sort by resid within each group
    result_solution_c1 = result_solution_c1.sort_values(by=["group_id", "resid"])

    # Reconstruct ID
    result_solution_c1 = result_solution_c1.reset_index()
    result_solution_c1["ID"] = (
        result_solution_c1["group_id"].astype(str)
        + "_"
        + (result_solution_c1["resid"]).astype(str)
    )

    # Remove C1'_ prefix from column names
    result_solution_c1.columns = [
        col.replace("C1'_", "") for col in result_solution_c1.columns
    ]

    # Move ID to front
    cols = list(result_solution_c1.columns)
    cols.remove("ID")
    cols = ["ID"] + cols
    result_solution_c1 = result_solution_c1[cols]
    return result_solution_c1


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file", type=str, help="Path to the input CSV file in long format"
    )
    parser.add_argument(
        "output_file_prefix", type=str, help="Prefix for the output CSV files"
    )
    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_file_prefix = args.output_file_prefix
    result_solution = pd.read_csv(input_file)
    result_solution_c1 = to_wide(result_solution)

    result_solution_c1.to_csv(
        f"{output_file_prefix}_solution_C1prime.csv",
        float_format="%.3f",
        index=False,
    )

import pandas as pd

def reduce_xyz_precision(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    
    # Find all columns containing "x", "y", or "z" (case-insensitive)
    xyz_cols = [col for col in df.columns if any(axis in col.lower() for axis in ["x", "y", "z"])]
    
    # Round those columns to 2 decimal places
    df[xyz_cols] = df[xyz_cols].round(2)
    
    # Save to CSV with minimal formatting to reduce size
    df.to_csv(output_csv, index=False, float_format='%.2f')

# Example usage
reduce_xyz_precision("train_solution_allatom.csv", "train_solution_allatom.csv")

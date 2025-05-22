import pandas as pd
import numpy as np
import qnorm

def read_data(input_file):
    """Reads the input file and sets 'GeneID' as the index."""
    print(f"Loading file: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    if 'GeneID' not in df.columns:
        print("ERROR: 'GeneID' column not found!")
        return None

    # Set 'GeneID' as index for easier handling
    df.set_index('GeneID', inplace=True)
    return df

def normalize_by_std_dev(df):
    """Normalizes the TPM values by dividing each gene by its standard deviation across tissues.
    If the standard deviation is 0, set the values to 0
    """
    # Calculate the standard deviation for each gene across tissues
    std_dev = df.std(axis=1, skipna=True)

    # Replace standard deviation of 0 with 1 for calculation (so values are not divided by NaN)
    std_dev[std_dev == 0] = 1  # Set std_dev of 0 to 1 to avoid division by NaN

    # Normalize TPM values by dividing by the standard deviation of each gene
    normalized_df = df.div(std_dev, axis=0)

    return normalized_df

def quantile_normalization(df):

    normalized_quantile_df = qnorm.quantile_normalize((df))

    return normalized_quantile_df

def log2_transformation(df):
	
	return df.applymap(lambda x: np.log2(x + 1) if pd.notna(x) and x >= 0 else np.nan)


def calculate_tissue_averages(df):
    """Calculates the average value for each tissue by grouping replicates."""
    tissue_columns = df.columns
    tissue_groups = {}

    # Group columns by tissue type (e.g., 'heart', 'test') by removing the replicate part
    for col in tissue_columns:
        tissue_name = '_'.join(col.split('_')[:-1])  # Extract tissue name from 'heart_rep1' -> 'heart'

        if tissue_name not in tissue_groups:
            tissue_groups[tissue_name] = []
        tissue_groups[tissue_name].append(col)

    # Create a new dataframe to store tissue averages
    averaged_df = pd.DataFrame(index=df.index)

    # For each tissue group, calculate the average of the replicates
    for tissue_name, columns in tissue_groups.items():
        averaged_df[tissue_name] = df[columns].mean(axis=1, skipna=True)

    return averaged_df

def normalize_tpm(input_file, output_file):
    """Main function that reads the data, applies normalizations, and saves the result."""
    df = read_data(input_file)
    if df is None:
        return

    # Normalize by standard deviation
    normalized_df = normalize_by_std_dev(df)
    print(normalized_df)

    # Apply quantile normalization
    quantile_normalized_df = quantile_normalization(normalized_df)
    print(quantile_normalized_df)

    # Apply log2 transformation
    log_transformed_df = log2_transformation(quantile_normalized_df)
    print(log_transformed_df)

    # Calculate tissue averages
    averaged_df = calculate_tissue_averages(log_transformed_df)

    # Reset 'GeneID' as a column and save the result to a new file
    averaged_df.reset_index(inplace=True)
    averaged_df.to_csv(output_file, sep='\t', index=False)
    print(f"Averaged tissue values saved: {output_file}")

if __name__ == "__main__":
    input_file_rat = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/kallisto_outputs/combined_mapped.tsv"
    output_file_rat = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/kallisto_outputs/combined_mapped_normalized.tsv"
    
    normalize_tpm(input_file_rat, output_file_rat)


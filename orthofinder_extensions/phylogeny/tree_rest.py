import pandas as pd
import os
import argparse
import pickle


def load_orthogroups(file_path):
    """
    Loads the orthogroups table from a file.
    Expects a tab-separated values (TSV) file.
    """
    try:
        return pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None


def find_start_index(df, start_orthogroup):
    """
    Finds the starting index of the orthogroup in the dataframe.
    Returns the index of the row immediately following the specified orthogroup.
    """
    start_index = df[df.iloc[:, 0] == start_orthogroup].index
    return start_index[0] + 1 if not start_index.empty else None


def extract_gene_to_species(row, columns):
    """
    Extracts a mapping of genes to species from a single dataframe row.
    Assumes genes are listed in each column as comma-separated values.
    """
    gene_to_species = {}
    for column in columns:
        genes = row[column]
        if pd.notna(genes):  # Check if the cell is not empty
            for gene in genes.split(", "):  # Split genes by comma
                gene_to_species[gene] = column  # Map each gene to its species (column name)
    return gene_to_species


def classify_all_genes(df, start_index, score_dict):
    """
    Iterates over all rows of the dataframe starting from `start_index`
    and classifies the genes into CO, SIN, SOU, and IN relationships.
    """
    co_results, sin_results, sou_results, in_results = [], [], [], []

    for _, row in df.iloc[start_index:].iterrows():
        orthogroup_id = row.iloc[0]  # The ID of the current orthogroup
        gene_to_species = extract_gene_to_species(row, df.columns[1:])
        co_data, sin_data, sou_data, in_data = classify_genes(gene_to_species, score_dict, orthogroup_id)
        co_results.extend(co_data)
        sin_results.extend(sin_data)
        sou_results.extend(sou_data)
        in_results.extend(in_data)

    return co_results, sin_results, sou_results, in_results


def classify_genes(gene_to_species, score_dict, orthogroup_id):
    """
    Classifies genes into:
    - Conserved Orthologs (CO)
    - Special In-Paralogs (SIN)
    - Special Out-Paralogs (SOU)
    - In-Paralogs (IN)

    Classification is based on the number of genes and species in the orthogroup.
    """
    # Map species to their corresponding genes
    species_to_genes = {}
    for gene, species in gene_to_species.items():
        species_to_genes.setdefault(species, []).append(gene)

    co_data, sin_data, sou_data, in_data = [], [], [], []

    # len(gene_to_species) means just the number of genes in that row -> {'FBpp0148700': 'dgri', 'FBpp0161065': 'dmoj', 'FBpp0227963': 'dvir'}
    # len(species_to_genes) means just the number of species in that row -> {'dgri': ['FBpp0148700'], 'dmoj': ['FBpp0161065'], 'dvir': ['FBpp0227963']}

    # Case 1: Exactly 2 genes and 2 species -> CO
    if len(gene_to_species) == 2 and len(species_to_genes) == 2:
        genes = list(gene_to_species.keys())
        co_data.append([orthogroup_id, genes[0], gene_to_species[genes[0]], genes[1], gene_to_species[genes[1]]])

    # Case 2: Exactly 2 genes and 1 species -> IN
    elif len(gene_to_species) == 2 and len(species_to_genes) == 1:
        genes = list(gene_to_species.keys())
        in_data.append([orthogroup_id, genes[0], gene_to_species[genes[0]], genes[1], gene_to_species[genes[1]]])

    # Case 3: Exactly 3 genes and 1 species -> IN
    elif len(gene_to_species) == 3 and len(species_to_genes) == 1:
        genes = list(gene_to_species.keys())
        for i in range(len(genes)):
            for j in range(i+1,len(genes)):
                in_data.append([orthogroup_id, genes[i], gene_to_species[genes[i]], genes[j], gene_to_species[genes[j]]])


    # Case 4: Exactly 3 genes and 2 species -> 1 CO , 1 SIN, 1 SOU
    elif len(gene_to_species) == 3 and len(species_to_genes) == 2:
        # Identify the species with 2 genes and the one with 1 gene
        species_with_two_genes = [s for s, genes in species_to_genes.items() if len(genes) == 2][0]
        species_with_one_gene = [s for s, genes in species_to_genes.items() if len(genes) == 1][0]

        genes_in_two = species_to_genes[species_with_two_genes]
        gene_in_one = species_to_genes[species_with_one_gene][0]

        # Determine the best scoring gene pair
        best_pair, best_score = None, float("-inf")
        for gene in genes_in_two:
            pair = (gene, gene_in_one) if (gene, gene_in_one) in score_dict else (gene_in_one, gene)
            score = score_dict.get(pair, float("-inf"))
            if score > best_score:
                best_score = score
                best_pair = pair

        # Classify the genes based on relationships
        co_data.append(
            [orthogroup_id, best_pair[0], gene_to_species[best_pair[0]], best_pair[1], gene_to_species[best_pair[1]]])
        other_gene = [g for g in genes_in_two if g != best_pair[0]][0]
        sin_data.append(
            [orthogroup_id, best_pair[0], gene_to_species[best_pair[0]], other_gene, gene_to_species[other_gene]])
        sou_data.append(
            [orthogroup_id, other_gene, gene_to_species[other_gene], best_pair[1], gene_to_species[best_pair[1]]])

    # Case 4: Exactly 3 genes and 3 species -> CO
    elif len(gene_to_species) == 3 and len(species_to_genes) == 3:
        genes = list(gene_to_species.keys())
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                co_data.append(
                    [orthogroup_id, genes[i], gene_to_species[genes[i]], genes[j], gene_to_species[genes[j]]])

    return co_data, sin_data, sou_data, in_data


def append_to_file(df, file_path):
    """
    Appends the DataFrame to the file, ensuring headers are written only if the file is new.
    """
    # If the file already exists, append data without writing the header
    if os.path.exists(file_path):
        df.to_csv(file_path, sep='\t', header=False, index=False, mode='a')
    else:
        # If the file doesn't exist, write data with header
        df.to_csv(file_path, sep='\t', header=True, index=False, mode='w')

    # Danach kompletten DataFrame einlesen
    full_df = pd.read_csv(file_path, sep='\t')

    # Exakte Zeilenduplikate entfernen
    full_df = full_df.drop_duplicates()

    # Optional: Datei überschreiben mit der bereinigten Version
    full_df.to_csv(file_path, sep='\t', index=False)

    return full_df

def main():

    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Classify orthogroups and append results to files.")

    # Arguments for the dictionary (harmonic_dict) and file paths
    parser.add_argument("--arg1", type=str, required=True, help="File path for CO results.")
    parser.add_argument("--arg2", type=str, required=True, help="File path for IN results.")
    parser.add_argument("--arg3", type=str, required=True, help="File path for SIN results.")
    parser.add_argument("--arg4", type=str, required=True, help="File path for SOU results.")
    parser.add_argument("--start_index", type=str, required=True,
                        help="The orthogroup index from which to start processing.")
    parser.add_argument("--df", type=str, required=True,
                        help="JSON string representation of the DataFrame.")

    # Parse arguments
    args = parser.parse_args()

    # Load the harmonic_dict from the passed argument (JSON)
    harmonic_dict_path = "results/harmonic_means_by_id.pkl"

    # File paths for appending results
    co_path = args.arg1
    in_path = args.arg2
    sin_path = args.arg3
    sou_path = args.arg4
    orthogroups_df = pd.read_json(args.df, orient="split")
    start_orthogroup = args.start_index

    if os.path.exists(harmonic_dict_path):
        with open(harmonic_dict_path, 'rb') as f:
            harmonic_dict = pickle.load(f)
        print("Loaded harmonic_dict from pickle file.")


    # Find the starting index for processing
    start_index = find_start_index(orthogroups_df, start_orthogroup)
    if start_index is None:
        print(f"Orthogroup '{start_orthogroup}' not found.")
        return

    # Perform classifications
    co_results, sin_results, sou_results, in_results = classify_all_genes(orthogroups_df, start_index, harmonic_dict)

    # Create DataFrames for the results
    co_df = pd.DataFrame(co_results, columns=["Orthogroup", "Gene1", "Species1", "Gene2", "Species2"])
    sin_df = pd.DataFrame(sin_results, columns=["Orthogroup", "Gene1", "Species1", "Gene2", "Species2"])
    sou_df = pd.DataFrame(sou_results, columns=["Orthogroup", "Gene1", "Species1", "Gene2", "Species2"])
    in_df = pd.DataFrame(in_results, columns=["Orthogroup", "Gene1", "Species1", "Gene2", "Species2"])

    # Append results to the specified file paths
    append_to_file(co_df, co_path)
    append_to_file(sin_df, sin_path)
    append_to_file(sou_df, sou_path)
    append_to_file(in_df, in_path)


if __name__ == "__main__":
    main()
    print("finished with Tree_rest.py")

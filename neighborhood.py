import pandas as pd
import re
import os
import pickle
import networkx as nx
import csv


# 1. Load GTF files and extract relevant information
def load_gtf_genes(gtf_file):
    """
    Loads a GTF file and extracts relevant gene information.

    Parameters:
    gtf_file (str): Path to the GTF file.

    Returns:
    pd.DataFrame: A DataFrame containing the filtered gene data, sorted by chromosome and start position.
    """
    try:
        col_names = ["chromosome", "source", "feature_type", "start", "end", "score", "strand", "frame", "attribute"]
        genes_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=col_names)
        genes_df = genes_df[genes_df['feature_type'] == 'gene']  # Filter for gene entries only
        genes_df['gene_id'] = genes_df['attribute'].apply(lambda x: re.search(r'GeneID:(\d+)', x).group(1))
        genes_df['gene_symbol'] = genes_df['attribute'].apply(lambda x: re.search(r'gene_id "([^"]+)"', x).group(1))
        genes_df = genes_df[['chromosome', 'start', 'end', 'strand', 'gene_id', 'gene_symbol']]
        genes_df = genes_df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
        #print(genes_df)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {gtf_file} does not exist.")
    return genes_df


def read_blast(blast_file):
    """
    Reads a BLAST results file and returns a DataFrame.

    Parameters:
    blast_file (str): Path to the BLAST results file.

    Returns:
    pd.DataFrame: A DataFrame containing the BLAST results.
    """
    """columns = ['query_id', 'subject_id', 'percent_identity', 'alignment_length',
               'mismatches', 'gap_opens', 'query_start', 'query_end',
               'subject_start', 'subject_end', 'e_value', 'bit_score']"""

    columns = ['query_id', 'subject_id', 'query_start','query_end', 'qlen' , 'subject_start' ,  'subject_end' , 'slen' , 'percent_identity' , 'e_value' , 'bit_score' ]

    blast_df = pd.read_csv(
        blast_file,
        sep='\t',
        header=None,
        names=columns,
        on_bad_lines='error',
        dtype={0: str, 1: str}  # Spalte 0 und 1 als string einlesen
    )

    return blast_df


def create_relationship_dictionary(blast_df):
    """
    Creates a dictionary of relationships from the BLAST DataFrame.

    Parameters:
    blast_df (pd.DataFrame): A DataFrame containing BLAST results.

    Returns:
    dict: A dictionary where keys are query IDs and values are lists of tuples with subject IDs and percent identity.
    """
    relationship_dict = {}
    for _, row in blast_df.iterrows():
        query_id = str(row['query_id'])
        subject_id = str(row['subject_id'])
        percent_identity = float(row['percent_identity'])  # explizit in Python float umwandeln
        #print(query_id, subject_id, percent_identity)
        if query_id != subject_id:  # Avoid self-relations
            relationship_dict.setdefault(query_id, []).append((subject_id, percent_identity))
            relationship_dict.setdefault(subject_id, []).append((query_id, percent_identity))
    return relationship_dict


def count_groups_by_chromosome(genes_df, species_name):
    """
    Counts the number of genes per chromosome and filters chromosomes with more than 10 genes.

    Parameters:
    genes_df (pd.DataFrame): A DataFrame containing gene information.

    Returns:
    dict: A dictionary where keys are chromosomes and values are lists of gene IDs.
    """
    total_genes = len(genes_df)
    # Group by the 'chromosome' column and count the number of genes in each group
    chromosome_count = genes_df.groupby('chromosome').size()

    # Convert the result to a DataFrame for better visualization
    chromosome_count_df = chromosome_count.reset_index(name='gene_count')
    chromosomes_more_than_10_genes = chromosome_count_df[chromosome_count_df['gene_count'] > 10]
    num_chromosomes_more_than_10_genes = chromosomes_more_than_10_genes.shape[0]
    avg_genes_more_than_10 = chromosomes_more_than_10_genes['gene_count'].mean()
    filtered_chromosomes = chromosomes_more_than_10_genes['chromosome'].tolist()
    filtered_genes = genes_df[genes_df['chromosome'].isin(filtered_chromosomes)]
    filtered_genes = filtered_genes.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    # print(filtered_genes)
    gene_dict = create_gene_dictionary(filtered_genes)
    # print(gene_dict)

    total_groups = chromosome_count_df.shape[0]
    print(f"\nSpecies: {species_name}")
    print(f"Total of genes: {total_genes} ")
    print(f"Total of chromosomes: {total_groups}")
    print(f"Number of genes per chromosome:\n{chromosome_count_df}")
    print(f"Chromosomes with more than 10 genes: {num_chromosomes_more_than_10_genes}")
    print(f"Number of genes per chromosome:\n{chromosomes_more_than_10_genes}")
    print(f"Average genes per chromosome (for chromosomes with >10 genes): {avg_genes_more_than_10:.2f}")
    #print(gene_dict)
    return gene_dict


def create_gene_dictionary(genes_df):
    """
    Creates a dictionary where the key is the chromosome and the value is a list of gene IDs.

    Parameters:
    genes_df (pd.DataFrame): A DataFrame containing gene information.

    Returns:
    dict: A dictionary where keys are chromosomes and values are lists of gene IDs.
    """
    gene_dict = {}
    for _, row in genes_df.iterrows():
        chromosome = row['chromosome']
        gene_id = str(row['gene_id'])

        if chromosome not in gene_dict:
            gene_dict[chromosome] = []
        gene_dict[chromosome].append(gene_id)

    return gene_dict


def save_relationship_dictionary(relationship_dict, file_path='relationship_dict.pkl'):
    """
    Saves the relationship dictionary to a file.

    Parameters:
    relationship_dict (dict): A dictionary containing the gene relationships.
    file_path (str): Path to the file where the dictionary will be saved.
    """
    with open(file_path, 'wb') as file:
        pickle.dump(relationship_dict, file)


def load_relationship_dictionary(file_path):
    """
    Loads the relationship dictionary from a file.

    Parameters:
    file_path (str): Path to the file where the dictionary is saved.

    Returns:
    dict: A dictionary containing the gene relationships.
    """
    with open(file_path, 'rb') as file:
        relationship_dict = pickle.load(file)
    return relationship_dict


def create_neighborhood_dictionary(gene_dicts, tandem_dict, num_neighbors=10):
    """
    Creates a neighborhood dictionary where each gene is associated with its neighboring genes,
    skipping tandem genes unless they are the representative gene.

    Parameters:
    gene_dicts (dict): A dictionary where keys are species and values are dictionaries of genes by chromosome.
    tandem_dict (dict): A dictionary where each gene maps to a set of its tandem group members and the representative gene.
    num_neighbors (int): The number of neighboring genes to consider on each side (default is 10).

    Returns:
    dict: A dictionary where keys are gene IDs and values are lists of neighboring gene IDs.
    """
    half_neighbors = int(num_neighbors / 2)
    neighborhood_dict = {}
    count = 0

    for species, genes_by_chromosome in gene_dicts.items():
        for chromosome, genes in genes_by_chromosome.items():

            # Get neighbors for each gene
            for index, gene in enumerate(genes):
                neighborhood = []

                # Check if the gene is a tandem gene and if it's representative
                if gene in tandem_dict:
                    representative = tandem_dict[gene]["representative_gene"]
                    if gene != representative:
                        continue

                # Get left-side neighbors
                left_count = 0
                for i in range(1, len(genes)):
                    if left_count >= half_neighbors:
                        break
                    neighbor_index = index - i
                    if neighbor_index >= 0:
                        neighbor_gene = genes[neighbor_index]
                        if neighbor_gene in tandem_dict:
                            rep = tandem_dict[neighbor_gene]["representative_gene"]
                            if neighbor_gene != rep:
                                # print("ignorando1")
                                # count+=1
                                continue
                            else:
                                count += 1
                        neighborhood.append(neighbor_gene)
                        left_count += 1

                # Get right-side neighbors
                right_count = 0
                for i in range(1, len(genes)):
                    if right_count >= half_neighbors:
                        break
                    neighbor_index = index + i
                    if neighbor_index < len(genes):
                        neighbor_gene = genes[neighbor_index]
                        if neighbor_gene in tandem_dict:
                            rep = tandem_dict[neighbor_gene]["representative_gene"]
                            if neighbor_gene != rep:
                                # print("ignorando2")

                                continue
                            else:
                                count += 1
                        neighborhood.append(neighbor_gene)
                        right_count += 1

                # Add the neighborhood to the dictionary
                neighborhood_dict[gene] = neighborhood
    # print(count)
    print(len(neighborhood_dict))
    return neighborhood_dict


def load_tandem_dict(tandem_file):
    """
    Loads tandem gene information from a file into a dictionary, including the representative gene.

    Parameters:
    tandem_file (str): Path to the file containing tandem gene information.

    Returns:
    dict: A dictionary where each gene maps to a dictionary with tandem group members and the representative gene.
    """
    tandem_dict = {}

    if not os.path.exists(tandem_file):
        raise FileNotFoundError(f"The file {tandem_file} was not found.")

    with open(tandem_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        try:
            next(reader)  # Skip the header
        except StopIteration:
            return tandem_dict  # Return an empty dictionary if the file is empty

        for row in reader:
            genes = row[0].split(", ")  # Split genes in the tandem group
            species = row[1]  # Capture species (if needed later)
            representative_gene = row[2]  # Capture the representative gene

            if len(genes) < 2:
                raise ValueError(f"Invalid tandem group format: {row[0]}")

            for gene in genes:
                tandem_dict[gene] = {
                    "tandem_group": set(genes),
                    "representative_gene": representative_gene
                }

    return tandem_dict


def count_shared_genes(gene, gene_related, neighborhood_1, neighborhood_2, relationship_dict, csv_all_edges,
                       csv_unique_edges):
    """
    Counts the maximum number of unique gene pairs between two neighborhoods
    using the Hopcroft-Karp algorithm for bipartite matching.

    Parameters:
    neighborhood_1 (list): A list of genes in the first neighborhood.
    neighborhood_2 (list): A list of genes in the second neighborhood.
    relationship_dict (dict): A dictionary containing gene relationships.

    Returns:
    int: The maximum number of unique gene pairs between the two neighborhoods.
    """
    count = 0
    conjunto_vecindad_2 = set(neighborhood_2)

    # Create a bipartite graph
    B = nx.Graph()
    B.add_nodes_from(neighborhood_1, bipartite=0)  # Conjunto izquierdo
    B.add_nodes_from(neighborhood_2, bipartite=1)  # Conjunto derecho

    # Agregar aristas basadas en el diccionario de relaciones
    for gene1 in neighborhood_1:
        related_genes = relationship_dict.get(gene1, [])
        for gene2 in related_genes:
            related_gene = gene2[0]
            if related_gene in conjunto_vecindad_2 and gene1 not in conjunto_vecindad_2:  # Asegurar que gene2 está en el vecindario 2
                B.add_edge(gene1, related_gene)

    if len(B.edges) == 0:
        return 0, 0  # No hay aristas, no puede haber emparejamiento

    # Encontrar el máximo emparejamiento
    try:
        matching = nx.bipartite.maximum_matching(B, top_nodes=neighborhood_2)
    except:
        print(B.edges)
        print(neighborhood_2)
    # Contar solo pares únicos
    max_unique_pairs = len(matching) // 2

    return max_unique_pairs, len(B.edges)


def analyze_genes(gene_dicts, gene_to_species, neighborhoods, relationship_dict, output_file, csv_all_edges,
                  csv_unique_edges, batch_size=10000):
    """
    Analyzes genes, calculates shared neighbors, and writes the results to a file.

    Parameters:
    gene_dicts (dict): A dictionary of genes by species and chromosome.
    gene_to_species (dict): A dictionary mapping gene IDs to species names.
    neighborhoods (dict): A dictionary of neighborhoods by gene ID.
    relationship_dict (dict): A dictionary of gene relationships.
    output_file (str): The path to the output file where results will be saved.
    batch_size (int): The number of records to write at a time (default is 10000).
    """
    records = []  # List to store records temporarily
    saved_pairs = set()  # To avoid duplicates
    count = 0
    print(f'relationsghip_dict : {relationship_dict}')
    with open(output_file, 'a') as f:
        for species, chromosomes in gene_dicts.items():
            print(f"Analyzing species {species}")
            for chromosome, genes in chromosomes.items():
                for gene in genes:
                    max_shared_genes_per_species = {}
                    neighborhood = neighborhoods.get(gene, [])

                    # Get directly related genes from the relationship dictionary
                    gene_relations = relationship_dict.get(gene, [])
                    #print(f'gene relations {gene_relations}')

                    # Process each related gene
                    for relation in gene_relations:
                        related_gene = relation[0]
                        percent_identity = relation[1]
                        related_species = gene_to_species.get(related_gene, 'Unknown species')

                        # Skip processing the same gene
                        if related_gene == gene:
                            continue

                        # Get the neighborhood of the related gene
                        related_neighborhood = neighborhoods.get(related_gene, [])
                        if len(related_neighborhood) == 0:
                            continue
                        shared_genes, all_shared_genes = count_shared_genes(gene, related_gene, neighborhood,
                                                                            related_neighborhood, relationship_dict,
                                                                            csv_all_edges, csv_unique_edges)
                        # print(shared_genes)
                        # Update the maximum shared genes by species
                        if related_species not in max_shared_genes_per_species:
                            max_shared_genes_per_species[related_species] = (
                            related_gene, shared_genes, all_shared_genes)
                        else:
                            current_gene, current_count, current_all_count = max_shared_genes_per_species[
                                related_species]
                            if shared_genes > current_count:
                                max_shared_genes_per_species[related_species] = (
                                related_gene, shared_genes, all_shared_genes)
                            if all_shared_genes > current_all_count and shared_genes < current_count:
                                count += 1
                            # if related_gene in ["FBgn0146602","FBgn0140504"] and gene =="FBgn0031558": #or all_shared_genes>20:
                            # print(gene, related_gene,related_species,all_shared_genes,shared_genes)
                            # print(gene, related_gene,related_species,current_gene,current_all_count,current_count,"\n")

                    # Add non-duplicate records and write in batches
                    for species_aux, (gene_aux, count, count_aux) in max_shared_genes_per_species.items():
                        if count > 0 and species != species_aux:
                            current_pair = (gene, species, gene_aux, species_aux)
                            reverse_pair = (gene_aux, species_aux, gene, species)
                            if current_pair not in saved_pairs and reverse_pair not in saved_pairs:
                                records.append(f"{gene}\t{species}\t{gene_aux}\t{species_aux}\t{count}\n")
                                # print(records)
                                saved_pairs.add(current_pair)

                    # Write to file if the batch size is reached
                    if len(records) >= batch_size:
                        print("saving...")
                        f.writelines(records)
                        records = []  # Clear list after writing

        # Write any remaining records that were not saved
        if records:
            f.writelines(records)
    # print(count)


def check_or_create_file(file, header):
    """
    Checks if a file exists; if not, creates it and adds the header.

    Parameters:
    file (str): Path to the file.
    header (str): Header to write to the file if it is created.
    """
    # Check if the file exists
    if not os.path.isfile(file):
        # If it doesn't exist, create it and add the header
        with open(file, 'w') as f:
            f.write(header + '\n')


def main():
    """
    Main function that loads data, processes it, and analyzes genes.
    """
    # Load the GTF files for the species
    genes_df_Canis = load_gtf_genes(r'/storage/EasyVectorOmics/FastQ_GSE125483_JK/gtf/gtf_Canis.gtf').reset_index(drop=True)
    genes_df_Macaca = load_gtf_genes(r'/storage/EasyVectorOmics/FastQ_GSE125483_JK/gtf/gtf_Macaca.gtf').reset_index(drop=True)
    genes_df_Mus = load_gtf_genes(r'/storage/EasyVectorOmics/FastQ_GSE125483_JK/gtf/gtf_Mus.gtf').reset_index(drop=True)
    genes_df_Rat = load_gtf_genes(r'/storage/EasyVectorOmics/FastQ_GSE125483_JK/gtf/gtf_Rat.gtf').reset_index(drop=True)
    # genes_df_dvir = load_gtf_genes('material/gtf_files/dvir-all-r1.07.gtf').reset_index(drop=True)
    # genes_df_dsec = load_gtf_genes('material/gtf_files/dsec-all-r1.3.gtf').reset_index(drop=True)
    # genes_df_dpse = load_gtf_genes('material/gtf_files/dpse-all-r3.04.gtf').reset_index(drop=True)
    # genes_df_dper = load_gtf_genes('material/gtf_files/dper-all-r1.3.gtf').reset_index(drop=True)
    # genes_df_dmoj = load_gtf_genes('material/gtf_files/dmoj-all-r1.04.gtf').reset_index(drop=True)
    # genes_df_dmel = load_gtf_genes('material/gtf_files/dmel-all-r6.24.gtf').reset_index(drop=True)
    # genes_df_dgri = load_gtf_genes('material/gtf_files/dgri-all-r1.05.gtf').reset_index(drop=True)
    # genes_df_dere = load_gtf_genes('material/gtf_files/dere-all-r1.05.gtf').reset_index(drop=True)

    species = {
        'Canis_Lupus': genes_df_Canis,
        'Macaca': genes_df_Macaca,
        'Mus': genes_df_Mus,
        'Rat': genes_df_Rat
        # 'dgri': genes_df_dgri,
        # 'dmoj': genes_df_dmoj,
        # 'dper': genes_df_dper,
        # 'dpse': genes_df_dpse,
        # 'dsec': genes_df_dsec,
        # 'dvir': genes_df_dvir,
        # 'dana': genes_df_dana,
        # 'dwil': genes_df_dwil
    }

    # Create gene dictionaries per species
    gene_dicts = {species_name: count_groups_by_chromosome(genes_df, species_name) for species_name, genes_df in
                  species.items()}

    # Create gene-to-species mapping
    gene_to_species = {gene: species_name for species_name, genes_df in species.items() for gene in genes_df['gene_id']}

    tandem_file = "/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tandems/Tandems_output.tsv"
    tandem_dict = load_tandem_dict(tandem_file)
    #print(tandem_dict)

    # Create neighborhoods
    neighborhoods = create_neighborhood_dictionary(gene_dicts, tandem_dict)
    print(len(neighborhoods))

    path_relationship_dict = 'material/clean_relationships.pkl'

    csv_all_edges = r"D:\EasyVectorOmics\TEST\results/all_edges.csv"
    csv_unique_edges = r"D:\EasyVectorOmics\TEST\results/unique_edges.csv"

    # Comprobar si ya existe el archivo de diccionario
    if os.path.exists(path_relationship_dict):
        relationship_dict = load_relationship_dictionary(path_relationship_dict)
    else:
        blast_df = read_blast(r'/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/blast/blast_results.tsv')  # Cargar datos de BLAST
        relationship_dict = create_relationship_dictionary(blast_df)
        #save_relationship_dictionary(relationship_dict, path_relationship_dict)

    # Check or create output file
    header_orth = "Gene1\tSpecies1\tGene2\tSpecies2\tCount"
    orthologues_nbh_path = r"/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/neighborhood/neighborhood_results.tsv"

    check_or_create_file(orthologues_nbh_path, header_orth)

    # Analyze genes and save the results
    analyze_genes(gene_dicts, gene_to_species, neighborhoods, relationship_dict, orthologues_nbh_path, csv_all_edges,
                  csv_unique_edges)


if __name__ == "__main__":
    main()

import re
from Bio import Phylo
import pandas as pd
import os
import pickle
import csv
from collections import defaultdict
import subprocess
import tempfile
import sys
# List of species that can appear in the tree
species_list = []
allowed_orthologs = None


def check_or_create_file(file, header):
    # Check if the file exists
    if not os.path.isfile(file):
        # If it doesn't exist, create it and add the header
        with open(file, 'w') as f:
            f.write(header + '\n')



def validate_args(args):
    """
    Validates command-line arguments.

    Args:
        args (list): Command-line arguments.

    Returns:
        tuple: GEO ID and output directory if valid, otherwise exits the program.
    """
    if len(args) != 2:
        print("\nInvalid Arguments")
        print("\nUsage: python script.py <method>")
        print("\nmethod options: 'standard' , 'majority' , 'whitelist'")
        sys.exit(1)

    method = args[1]


    return method

def write_file(classifications, co_path, in_path, sin_path, sout_path, out_path):
    """
    Schreibt die Klassifikationen in CSV-Dateien. Jede Klassifikation wird in eine separate Datei geschrieben,
    ohne die bestehenden Daten zu überschreiben.

    :param classifications: Dictionary mit den Klassifikationen.
    :param orthogroup: Die ID der Orthogruppe (wird ggf. in die Ausgabe geschrieben).
    :param co_path: Pfad zur CSV-Datei für 'conserved_orthologs'.
    :param in_path: Pfad zur CSV-Datei für 'inparalogs'.
    :param sin_path: Pfad zur CSV-Datei für 'special_inparalogs'.
    :param out_path: Pfad zur CSV-Datei für 'outparalogs'.
    :param sout_path: Pfad zur CSV-Datei für 'special_outparalogs'.
    """

    # Helper function to write data to a CSV file (in append mode)
    def append_to_csv(file_path, data):
        # "Open the file in append mode ('a') if it already exists
        with open(file_path, 'a', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            # If the file is empty (i.e., the header is missing), then write the header.
            if f.tell() == 0:
                writer.writerow(['Orthogroup', 'Gene1', 'Species1', 'Gene2', 'Species2'])  # Header line
            for entry in data:
                writer.writerow(entry)


    # write the different classifications in their respective files
    append_to_csv(co_path, classifications['conserved_orthologs'])
    append_to_csv(in_path, classifications['inparalogs'])
    append_to_csv(sin_path, classifications['special_inparalogs'])
    append_to_csv(out_path, classifications['outparalogs'])
    append_to_csv(sout_path, classifications['special_outparalogs'])


def extract_protein_combinations(directory):
    protein_dict = {}

    # Iterate through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".tsv"):
            file_path = os.path.join(directory, filename)

            # Open and read the TSV file
            with open(file_path, "r") as tsv_file:
                tsv_reader = csv.reader(tsv_file, delimiter="\t")

                # Extract header and format the name of the 3rd column
                header = next(tsv_reader)
                species_name_col3 = header[2].replace(".", "_") + "_"

                for row in tsv_reader:
                    # Ensure the row has at least four columns
                    if len(row) >= 4:
                        # Format the 2nd column's value
                        species_name_col2 = row[1].replace(".", "_") + "_"

                        # Extract and clean proteins from columns 3 and 4
                        proteins_col3 = [protein.strip() for protein in row[2].split(",")]
                        proteins_col4 = [protein.strip() for protein in row[3].split(",")]

                        # Generate all unique combinations of proteins between column 3 and column 4
                        for protein1 in proteins_col3:
                            for protein2 in proteins_col4:
                                # Add formatted prefixes
                                prefixed_protein1 = species_name_col3 + protein1
                                prefixed_protein2 = species_name_col2 + protein2

                                # Store the combination in the dictionary
                                combination = (prefixed_protein1, prefixed_protein2)
                                protein_dict[combination] = None  # Value is irrelevant

    return protein_dict


def call_tree_rest_script(str_arg1, str_arg2, str_arg3, str_arg4, df_arg, index_str):
    """
    Calls Tree_rest.py with the given arguments.

    :param dictionary_arg: A Python dictionary to be passed as an argument.
    :param str_arg1: First string argument (e.g. file path for CO results).
    :param str_arg2: Second string argument (e.g. file path for IN results).
    :param str_arg3: Third string argument (e.g. file path for SIN results).
    :param str_arg4: Fourth string argument (e.g. file path for SOU results).
    :param df_arg: DataFrame to be passed as an argument.
    :param index_str: The index string for the orthogroup (as a string).
    """

    # save the dataframe as temp file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as temp_df_file:
        df_arg.to_json(temp_df_file, orient="split")
        df_file_path = temp_df_file.name

    # Build the command
    command = [
        "python3", "tree_rest.py",
        "--arg1", str_arg1,
        "--arg2", str_arg2,
        "--arg3", str_arg3,
        "--arg4", str_arg4,
        "--df", df_file_path,  # Übergabe des Dateipfads
        "--start_index", index_str
    ]

    # Run the script
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running Tree_rest.py: {e}")


#--------------------Section_1----------------------------------------------------------#


""" This section is for the Tree Algorithm to find the conserved orthologs only"""

# Function to parse Newick tree structure from a string
def parse_newick_tree(newick_file):
    """
    Parses a Newick formatted tree string.

    Args:
        newick_str (str): A string representing a Newick tree.

    Returns:
        Phylo.BaseTree.Tree: A parsed tree object.

    Function:
        Converts a Newick formatted string into a tree object using the Bio.Phylo library.
    """
    #newick_file = io.StringIO(newick_str)
    #return Phylo.read(newick_file, "newick")

    tree = Phylo.read(newick_file, "newick")
    return tree


# Function to extract species name from an entry
def get_species_name(entry_name):
    """
    Extracts the species name from a given entry.

    Args:
        entry_name (str): The name of the entry, typically a leaf node name in the tree.

    Returns:
        str or None: The species name if found in the entry, otherwise None.

    Function:
        Searches for known species identifiers in the entry name using regex.
    """
    for species in species_list:
        if re.search(species, entry_name):
            return species
    return None


# Function to collect all leaf nodes (terminals) in a subtree
def get_all_leafs(clade):
    """
    Collects all leaf nodes (terminals) of a given clade.

    Args:
        clade (Phylo.BaseTree.Clade): A clade object representing part of the tree.

    Returns:
        list: A list of leaf nodes (Phylo.BaseTree.Clade) in the given clade.

    Function:
        Uses the get_terminals() method to retrieve all leaf nodes of the clade.
    """

    return clade.get_terminals()


# Function to create connections between two lists of leaf nodes
def add_connections_as_individual_entries(left_leafs, right_leafs, scoring_matrix):
    """
    Creates connections between leaf nodes from two lists based on species and scoring.

    Args:
        left_leafs (list): A list of leaf nodes (Phylo.BaseTree.Clade) from the left subtree.
        right_leafs (list): A list of leaf nodes (Phylo.BaseTree.Clade) from the right subtree.
        scoring_matrix (dict): A dictionary containing scores for leaf node pairs.

    Returns:
        list: A list of tuples, each representing a pair of connected leaf node names (str).

    Function:
        Matches species between the left and right leaf lists, either directly or using scores to find the best pair.
    """
    connections = []
    left_species_map = {}
    right_species_map = {}

    # Collect species from the left subtree
    for left_leaf in left_leafs:
        left_species = get_species_name(left_leaf.name)
        if left_species not in left_species_map:
            left_species_map[left_species] = []
        left_species_map[left_species].append(left_leaf)

    # Collect species from the right subtree
    for right_leaf in right_leafs:
        right_species = get_species_name(right_leaf.name)
        if right_species not in right_species_map:
            right_species_map[right_species] = []
        right_species_map[right_species].append(right_leaf)

    # Create connections based on species mapping and scoring
    for left_species, left_clade_list in left_species_map.items():
        for right_species, right_clade_list in right_species_map.items():
            if right_species != left_species:
                if len(left_clade_list) == 1 and len(right_clade_list) == 1:
                    connections.append((left_clade_list[0].name, right_clade_list[0].name))
                else:
                    best_pair = None
                    best_score = float('-inf')
                    for left_clade in left_clade_list:
                        for right_clade in right_clade_list:
                            pair_score = calculate(left_clade, right_clade, scoring_matrix)
                            if pair_score > best_score:
                                best_score = pair_score
                                best_pair = (left_clade.name, right_clade.name)
                    if best_pair:
                        connections.append(best_pair)

    return connections


def add_connections_as_individual_entries_majority_rule(left_leafs, right_leafs, scoring_matrix):
    """
    Creates connections between leaf nodes from two lists based on species and scoring.

    Args:
        left_leafs (list): A list of leaf nodes (Phylo.BaseTree.Clade) from the left subtree.
        right_leafs (list): A list of leaf nodes (Phylo.BaseTree.Clade) from the right subtree.
        scoring_matrix (dict): A dictionary containing scores for leaf node pairs.

    Returns:
        list: A list of tuples, each representing a pair of connected leaf node names (str).

    Function:
        Matches species between the left and right leaf lists, either directly or using scores to find the best pair.
    """
    connections = []
    left_species_map = {}
    right_species_map = {}

    pair_counts_left = defaultdict(lambda: defaultdict(int))
    pair_counts_right = defaultdict(lambda: defaultdict(int))

    # Collect species from the left subtree
    for left_leaf in left_leafs:
        left_species = get_species_name(left_leaf.name)
        if left_species not in left_species_map:
            left_species_map[left_species] = []
        left_species_map[left_species].append(left_leaf)

    # Collect species from the right subtree
    for right_leaf in right_leafs:
        right_species = get_species_name(right_leaf.name)
        if right_species not in right_species_map:
            right_species_map[right_species] = []
        right_species_map[right_species].append(right_leaf)

    # Create connections based on species mapping and scoring
    for left_species, left_clade_list in left_species_map.items():
        for right_species, right_clade_list in right_species_map.items():
            if left_species != right_species :
                if len(left_clade_list) == 1 and len(right_clade_list) == 1 :
                    connections.append((left_clade_list[0].name, right_clade_list[0].name))
                else:
                    best_pair = None
                    best_score = float('-inf')
                    for left_clade in left_clade_list:
                        for right_clade in right_clade_list:
                            pair_score = calculate(left_clade, right_clade, scoring_matrix)
                            if pair_score > best_score :
                                best_score = pair_score
                                best_pair = (left_clade.name, right_clade.name)
                    if best_pair:
                        left_species_id = best_pair[0]
                        right_species_id = best_pair[1]
                        pair_counts_left[left_species][left_species_id] += 1
                        pair_counts_right[right_species][right_species_id] += 1


    max_values_dict_left = {
        species: max(gen_dict, key=gen_dict.get)
        for species, gen_dict in pair_counts_left.items()
    }

    max_values_dict_right = {
        species: max(gen_dict, key=gen_dict.get)
        for species, gen_dict in pair_counts_right.items()
    }

    for left_species, left_gen in max_values_dict_left.items():
        for right_species, right_gen in max_values_dict_right.items():
            connections.append((left_gen, right_gen))


    return connections



def add_connections_as_individual_entries_whitelist(left_leafs, right_leafs, scoring_matrix):
    """
    Creates connections between leaf nodes from two lists based on species and scoring.

    Args:
        left_leafs (list): A list of leaf nodes (Phylo.BaseTree.Clade) from the left subtree.
        right_leafs (list): A list of leaf nodes (Phylo.BaseTree.Clade) from the right subtree.
        scoring_matrix (dict): A dictionary containing scores for leaf node pairs.

    Returns:
        list: A list of tuples, each representing a pair of connected leaf node names (str).

    Function:
        Matches species between the left and right leaf lists, either directly or using scores to find the best pair.
    """
    connections = []
    left_species_map = {}
    right_species_map = {}
    print(allowed_orthologs)

    # Collect species from the left subtree
    for left_leaf in left_leafs:
        left_species = get_species_name(left_leaf.name)
        if left_species not in left_species_map:
            left_species_map[left_species] = []
        left_species_map[left_species].append(left_leaf)

    # Collect species from the right subtree
    for right_leaf in right_leafs:
        right_species = get_species_name(right_leaf.name)
        if right_species not in right_species_map:
            right_species_map[right_species] = []
        right_species_map[right_species].append(right_leaf)

    # Create connections based on species mapping and scoring
    for left_species, left_clade_list in left_species_map.items():
        for right_species, right_clade_list in right_species_map.items():
            if left_species != right_species:
                if len(left_clade_list) == 1 and len(right_clade_list) == 1:
                    pair = (left_clade_list[0].name, right_clade_list[0].name)
                    # Check if the pair is in scoring_matrix (normal or reversed order)
                    if pair in allowed_orthologs or (pair[1], pair[0]) in allowed_orthologs:
                        connections.append(pair)
                else:
                    best_pair = None
                    best_score = float('-inf')
                    for left_clade in left_clade_list:
                        for right_clade in right_clade_list:
                            pair = (left_clade.name, right_clade.name)
                            reverse_pair = (right_clade.name, left_clade.name)
                            if pair in allowed_orthologs or reverse_pair in allowed_orthologs:
                                pair_score = calculate(left_clade, right_clade, scoring_matrix)
                                if pair_score > best_score:
                                    best_score = pair_score
                                    best_pair = pair
                    if best_pair:
                        connections.append(best_pair)

    return connections



# Function to calculate the score for a pair of leaf nodes
def calculate(left_leaf, right_leaf, scoring_matrix):
    """
    Calculates a score for a pair of leaf nodes using a scoring matrix.

    Args:
        left_leaf (Phylo.BaseTree.Clade): The left leaf node.
        right_leaf (Phylo.BaseTree.Clade): The right leaf node.
        scoring_matrix (dict): A dictionary with scores for leaf pairs.

    Returns:
        float: The score for the given leaf pair, or 0 if not found in the scoring matrix.

    Function:
        Uses the scoring matrix to retrieve the harmonic mean for the leaf pair.
    """
    left_leaf = str(left_leaf.name).strip()
    right_leaf = str(right_leaf.name).strip()
    harmonic_mean = scoring_matrix.get((left_leaf, right_leaf), 0)

    if harmonic_mean == 0:
        harmonic_mean = scoring_matrix.get((right_leaf, left_leaf), 0)

    # If still not found, write to a file and raise an exception

    return harmonic_mean

# Main function to process the tree and collect connections
def process_tree(tree, orthogroup, scoring_matrix,duplication_dict, method):
    """Processes the given tree to extract conserved ortholog connections.

    Args:
        tree (Phylo.BaseTree.Tree): A parsed tree object.
        orthogroup (str): The identifier for the orthogroup being processed.
        scoring_matrix (dict): A dictionary containing scores for leaf node pairs.

    Returns:
        list: A list of tuples, each representing a conserved ortholog connection.

    Function:
        Traverses the tree, skipping duplication nodes, and collects connections between subtrees.
    """
    result_connections = []

    # Mapping methods to functions
    method_function_map = {
        "standard": add_connections_as_individual_entries,
        "majority": add_connections_as_individual_entries_majority_rule,
        "whitelist": add_connections_as_individual_entries_whitelist,
    }

    add_connections_function = method_function_map.get(method)

    def traverse(clade):
        if clade.name in duplication_dict.get(orthogroup, []):
            for child in clade.clades:
                if not child.is_terminal():
                    traverse(child)
            return

        left_leafs = get_all_leafs(clade.clades[0]) if len(clade.clades) > 0 else []
        right_leafs = get_all_leafs(clade.clades[1]) if len(clade.clades) > 1 else []

        if left_leafs and right_leafs:
            connections = add_connections_function(left_leafs, right_leafs, scoring_matrix)
            result_connections.extend(connections)

        for child in clade.clades:
            if not child.is_terminal():
                traverse(child)

    traverse(tree.root)
    return result_connections

#-------------------------------Section_2----------------------------------------------#

""" This Section contains the algorithm to find the other relationships:
    CO, IN , SIN, OUT, SOU using the Tree results """


# Funktion zur Klassifikation von Genen
def classify_genes(csv_df, conserved_orthologs, target_orthogroup):


    def is_conserved_gene(gene):
        """Prüfen, ob ein Gen in einem der Tupel als Substring vorkommt."""

        return any(gene in element for tup in conserved_orthologs for element in tup)


    # Klassifikation starten
    classifications = {
        'inparalogs': [],
        'special_inparalogs': [],
        'outparalogs': [],
        'special_outparalogs': [],
        'conserved_orthologs': [],
    }

    # Hole die Zeile aus dem DataFrame, die die Ziel-Orthogruppe enthält
    target_row = csv_df[csv_df['Orthogroup'] == target_orthogroup]

    if not target_row.empty:
        # Extrahiere die Spaltennamen (Spezies)
        species_names = target_row.columns[1:]

        # Extrahiere die Gene aus der entsprechenden Zeile
        genes_in_orthogroup = defaultdict(list)
        for species, genes in zip(species_names, target_row.iloc[0, 1:]):
            if pd.notna(genes):  # Leere Zellen ignorieren
                gene_list = [gene.strip() for gene in genes.split(',')]
                genes_in_orthogroup[species].extend(gene_list)


    # Konservierte Orthologe mit Speziesinformationen füllen
    for tup in conserved_orthologs:
        protein_id1, species1 = None, None
        protein_id2, species2 = None, None


        # Bearbeite den ersten String im Tupel
        for species, ids in genes_in_orthogroup.items():
            for id in ids:
                if protein_id1 is None and id in tup[0]:  # Nur speichern, wenn noch nichts gefunden
                    protein_id1 = id
                    species1 = species

        # Bearbeite den zweiten String im Tupel
        for species, ids in genes_in_orthogroup.items():
            for id in ids:
                if protein_id2 is None and id in tup[1]:  # Nur speichern, wenn noch nichts gefunden
                    protein_id2 = id
                    species2 = species

        # Wenn beide Gene gefunden wurden, füge sie zur Liste hinzu

        classifications['conserved_orthologs'].append((target_orthogroup, protein_id1, species1, protein_id2, species2))

    # Innerhalb einer Spezies klassifizieren
    for species, genes_in_species in genes_in_orthogroup.items():
        for idx1, gene1 in enumerate(genes_in_species):
            for gene2 in genes_in_species[idx1 + 1:]:
                pair = (target_orthogroup, gene1, species, gene2, species)
                if not is_conserved_gene(gene1) and not is_conserved_gene(gene2):
                    classifications['inparalogs'].append(pair)
                else:
                    classifications['special_inparalogs'].append(pair)

    # Zwischen Spezies klassifizieren
    species_pairs = list(genes_in_orthogroup.keys())
    for idx2, species1 in enumerate(species_pairs):
        for species2 in species_pairs[idx2 + 1:]:
            for gene1x in genes_in_orthogroup[species1]:
                for gene2x in genes_in_orthogroup[species2]:
                    pair = (target_orthogroup, gene1x, species1, gene2x, species2)
                    if not is_conserved_gene(gene1x) and not is_conserved_gene(gene2x):
                        classifications['outparalogs'].append(pair)
                    elif not any((gene1x in t[0] and gene2x in t[1]) or (gene1x in t[1] and gene2x in t[0]) for t in
                                 conserved_orthologs):
                        classifications['special_outparalogs'].append(pair)

    return classifications


def main():
    global allowed_orthologs
    method = validate_args(sys.argv)
    orthologs_directory = "proteom/OrthoFinder/Results_Apr03/Orthologues"
    pickle_file = "results/harmonic_means_by_prefix.pkl"
    Duplications_file = "proteom/OrthoFinder/Results_Apr03/Gene_Duplication_Events/Duplications.tsv"
    tree_files = "proteom/OrthoFinder/Results_Apr03/Resolved_Gene_Trees"
    df_orthogroups = pd.read_csv('proteom/OrthoFinder/Results_Apr03/Orthogroups/Orthogroups.tsv', sep='\t')
    species_names = df_orthogroups.columns[1:]
    species_list.extend(species_names)


    # Load the harmonic_dict from pickle file if it exists
    if os.path.exists(pickle_file):
        with open(pickle_file, 'rb') as f:
            harmonic_dict = pickle.load(f)
        print("Loaded harmonic_dict from pickle file.")

    if method == 'whitelist':
        allowed_orthologs = extract_protein_combinations(orthologs_directory)


    df = pd.read_csv(Duplications_file, sep='\t', header=None)
    Duplications_df = df.groupby(0)[2].apply(list).to_dict()



    header_orth = "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2"
    header_para = "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2"


    th = 2

    folder = "whitelist"  # oder "majority", "whitelist", etc.

    co_path = f"results/tree/{folder}/complete_tree_conserved_orthologs_{th}.tsv"
    in_path = f"results/tree/{folder}/complete_tree_in_paralogs_{th}.tsv"
    sin_path = f"results/tree/{folder}/complete_tree_special_in_paralogs_{th}.tsv"
    sout_path = f"results/tree/{folder}/complete_tree_special_out_paralogs_{th}.tsv"
    out_path = f"results/tree/{folder}/complete_tree_out_paralogs_{th}.tsv"






    check_or_create_file(co_path, header_orth)
    check_or_create_file(in_path, header_para)
    check_or_create_file(sin_path, header_para)
    check_or_create_file(sout_path, header_para)
    check_or_create_file(out_path, header_para)

    sorted_files = sorted(os.listdir(tree_files))

    cont = 0

    for filename in sorted_files:
        if cont % 100 == 0:
            print(f"Processed {cont} files")
        cont += 1
        print(cont)
        file_path = os.path.join(tree_files, filename)
        tree = parse_newick_tree(file_path)
        orthogroup = filename.split("_")[0]
        result_conserved_orthologs = process_tree(tree, orthogroup, harmonic_dict, Duplications_df, method)
        result_conserved_orthologs = sorted(result_conserved_orthologs, key=lambda x: x[0])
        #for entry in result_conserved_orthologs:
        #    print(entry)
        result = classify_genes(df_orthogroups, result_conserved_orthologs, orthogroup)
        write_file(result, co_path, in_path, sin_path, sout_path, out_path)

    call_tree_rest_script(co_path, in_path, sin_path, sout_path, df_orthogroups, orthogroup)
    print("finished")


if __name__ == "__main__":
    main()

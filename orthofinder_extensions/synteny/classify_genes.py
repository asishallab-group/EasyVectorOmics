import pandas as pd
import os
import csv
from collections import defaultdict

def check_or_create_file(file, header):
    # Check if the file exists, if not create it and add header
    if not os.path.isfile(file):
        with open(file, 'w') as f:
            f.write(header + '\n')

def classify_genes(csv_df, conserved_orthologs, target_orthogroup):


    def is_conserved_gene(gene):
        """Prüfen, ob ein Gen in einem der Tupel als Substring vorkommt."""

        return any(gene in element for tup in conserved_orthologs for element in tup)


    # Klassifikation starten
    classifications = {
        'inparalogs': [],
        'source_copy_inparalogs': [],
        'outparalogs': [],
        'source_copy_outparalogs': [],
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
        # print(tup)
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
    # print(classifications['conserved_orthologs'])
    # Innerhalb einer Spezies klassifizieren
    for species, genes_in_species in genes_in_orthogroup.items():
        for idx1, gene1 in enumerate(genes_in_species):
            for gene2 in genes_in_species[idx1 + 1:]:
                pair = (target_orthogroup, gene1, species, gene2, species)
                if not is_conserved_gene(gene1) and not is_conserved_gene(gene2):
                    classifications['inparalogs'].append(pair)
                else:
                    classifications['source_copy_inparalogs'].append(pair)

    # print(classifications['inparalogs'])
    # print(classifications['source_copy_inparalogs'])

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
                        classifications['source_copy_outparalogs'].append(pair)

    return classifications


def write_file(classifications, co_path, in_path, sin_path, sout_path, out_path):
    """
    Write the classifications into separate CSV files.
    """

    def append_to_csv(file_path, data):
        """Helper function to append data to a CSV file."""
        with open(file_path, 'a', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            if f.tell() == 0:  # Check if the file is empty to write header
                writer.writerow(['Orthogroup', 'Gene1', 'Species1', 'Gene2', 'Species2'])
            writer.writerows(data)

    # Escribir las clasificaciones en sus respectivos archivos
    append_to_csv(co_path, classifications['conserved_orthologs'])
    append_to_csv(in_path, classifications['inparalogs'])
    append_to_csv(sin_path, classifications['source_copy_inparalogs'])
    append_to_csv(out_path, classifications['outparalogs'])
    append_to_csv(sout_path, classifications['source_copy_outparalogs'])


# Cargar los datos
df_orthogroups = pd.read_csv('proteom/OrthoFinder/Results_Apr03/Orthogroups/Orthogroups.tsv', sep='\t')
df = pd.read_csv('results/neighborhood/neighborhood_results_formatted.tsv', sep='\t')

# Crear una lista de tuplas para conserved_orthologs
result_conserved_orthologs = defaultdict(list)
for _, row in df.iterrows():
    result_conserved_orthologs[row['Orthogroup']].append((row['Gene1'], row['Gene2']))

# Definir rutas de archivo
co_path = "results/neighborhood/neighborhood_conserved_orthologs.tsv"
in_path = "results/neighborhood/neighborhood_in_paralogs.tsv"
sin_path = "results/neighborhood/neighborhood_source_copy_in_paralogs.tsv"
sout_path = "results/neighborhood/neighborhood_source_copy_out_paralogs.tsv"
out_path = "results/neighborhood/neighborhood_out_paralogs.tsv"

# Verificar o crear archivos
check_or_create_file(co_path, "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2")
check_or_create_file(in_path, "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2")
check_or_create_file(sin_path, "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2")
check_or_create_file(sout_path, "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2")
check_or_create_file(out_path, "Orthogroup\tGene1\tSpecies1\tGene2\tSpecies2")

# Clasificar genes por Orthogroup
for orthogroup in df_orthogroups['Orthogroup'].unique():
    print("Procesando orthogroup ", orthogroup)
    result = classify_genes(df_orthogroups, result_conserved_orthologs[orthogroup], orthogroup)
    # print(result)
    write_file(result, co_path, in_path, sin_path, sout_path, out_path)


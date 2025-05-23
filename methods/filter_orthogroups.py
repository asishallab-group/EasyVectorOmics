import pandas as pd

# Load the files
ortologs_df = pd.read_csv("material/orthologs_filtered_dmel_dsec_dsim.tsv", sep="\t")
paralogs_df = pd.read_csv("material/paralogs_filtered_dmel_dsec_dsim.tsv", sep="\t")
all_orthogroups_df = pd.read_csv("material/Orthogroups.tsv", sep="\t")

# Create a set with all genes present in the paralogs file
paralog_genes = set(paralogs_df['Gene1']).union(set(paralogs_df['Gene2']))
orth_genes = set(ortologs_df['Gene1']).union(set(ortologs_df['Gene2']))

# Create a dictionary that maps each gene to its associated paralogs
gene_to_paralogs = {}
for _, row in paralogs_df.iterrows():
    # Add paralogs in both directions
    gene_to_paralogs.setdefault(row['Gene1'], set()).add(row['Gene2'])
    gene_to_paralogs.setdefault(row['Gene2'], set()).add(row['Gene1'])

gene_to_orth = {}
for _, row in ortologs_df.iterrows():
    # Add orthologs in both directions
    gene_to_orth.setdefault(row['Gene1'], set()).add(row['Gene2'])
    gene_to_orth.setdefault(row['Gene2'], set()).add(row['Gene1'])

# Group orthologs by Orthogroup
filtered_orthogroups = []

# Iterate over orthogroups only once
for og, group_df in ortologs_df.groupby("Orthogroup"):
    # Get the genes in the orthogroup
    genes_in_group = set(group_df['Gene1']).union(set(group_df['Gene2']))

    # Count how many orthologs have at least one related paralog
    orthologs_with_paralogs = 0
    for gene in genes_in_group:
        # Check if the gene has related paralogs
        related_paralogs = gene_to_paralogs.get(gene, set())
        only_par = related_paralogs - orth_genes
        
        # If there's at least one related paralog, count this ortholog
        if len(only_par) > 0:
            orthologs_with_paralogs += 1
    
    # If there are at least 4 orthologs with related paralogs, include the orthogroup
    if orthologs_with_paralogs >= 4:
        filtered_orthogroups.append(og)

# Print how many orthogroups meet the condition
print(f"Number of filtered orthogroups: {len(filtered_orthogroups)}")

# Filter the original ortholog DataFrame using the selected orthogroups
filtered_orth = ortologs_df[ortologs_df['Orthogroup'].isin(filtered_orthogroups)]
filtered_par = paralogs_df[paralogs_df['Orthogroup'].isin(filtered_orthogroups)]
filtered_families = all_orthogroups_df[all_orthogroups_df['Orthogroup'].isin(filtered_orthogroups)]

# Save results
filtered_orth.to_csv("material/filtered_orthologs.tsv", sep="\t", index=False)
filtered_par.to_csv("material/filtered_paralogs.tsv", sep="\t", index=False)
filtered_families.to_csv("material/filtered_families.tsv", sep="\t", index=False)

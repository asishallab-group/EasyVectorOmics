import pandas as pd

# 1. Einlesen der Dateien
map_df = pd.read_csv("D:/EasyVectorOmics/TEST/map_gen_protein_species_orthogroup.tsv", sep="\t", usecols=["GeneID", "ProteinID"], dtype=str)
tandem_df = pd.read_csv("D:/EasyVectorOmics/TEST/Tandems_output.tsv", sep="\t", usecols=["genids"], dtype=str)
dist_df = pd.read_csv(r"D:\EasyVectorOmics\TEST\results_OP\majority\distances_df.tsv", sep="\t", dtype=str)


# 2. Normalize ProteinID/GeneID



# 3. Mapping ProteinID → GeneID
protein_to_gene = dict(zip(map_df["ProteinID"], map_df["GeneID"]))

# 4. Erstelle Menge aller Tandem-GeneIDs
tandem_genes = set()
for entry in tandem_df["genids"]:
    genes = entry.split(",")
    tandem_genes.update(g.strip() for g in genes)

# 5. Nur paraloge betrachten
paralog_mask = dist_df["gene_type"] == "paralog"
paralog_df = dist_df[paralog_mask].copy()

# 6. Mappe ProteinID zu GeneID (und prüfen auf Fehlschläge)
paralog_df["gene_id"] = paralog_df["gene_id"].str.strip()
paralog_df["GeneID"] = paralog_df["gene_id"].map(protein_to_gene)

# Debug-Ausgabe: Wie viele ProteinIDs konnten nicht gemappt werden?
unmapped = paralog_df["GeneID"].isnull().sum()
print(f"{unmapped} ProteinIDs konnten nicht zu GeneIDs gemappt werden.")

# 7. Tandemprüfung
paralog_df["is_tandem"] = paralog_df["GeneID"].isin(tandem_genes)

# 8. Neue Spalte initialisieren mit False
dist_df["is_tandem"] = False

# 9. Setze True für passende gene_id-Einträge
dist_df.loc[paralog_df.index, "is_tandem"] = paralog_df["is_tandem"].values

# 10. Ausgabe, ob es überhaupt ein True gibt
if dist_df["is_tandem"].any():
    count = dist_df["is_tandem"].sum()
    print(f"{count} Tandem-Paraloge erkannt (is_tandem == True).")
else:
    print("Kein Tandem-Paralog erkannt.")

# 11. Speichern
dist_df.to_csv(r"D:\EasyVectorOmics\TEST\results_OP\majority\distances_tandem.tsv", sep="\t", index=False)

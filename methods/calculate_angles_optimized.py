import pandas as pd
import numpy as np

# ======================
# CONFIGURATION SECTION
# ======================
input_file = "results/distances_tandem.tsv"
centroid_file = "results/orthogroup_centroids"
pairwise_file = "results/paralog_pairwise_distances.tsv"

output_df_file = "results/distances_tandem_angle.tsv"
output_expr_file = "results/centroid_angle.tsv"
output_pairwise_file = "results/pairwise_distances_angle.tsv"


# ======================
# UTILITY FUNCTIONS
# ======================

def calculate_angle(vector, unit_diagonal):
    """
    Calculate the angle between a vector and the space diagonal unit vector.

    Args:
        vector (np.array): Input vector
        unit_diagonal (np.array): Unit vector of space diagonal

    Returns:
        float: Angle in degrees or np.nan if vector is zero
    """
    if np.linalg.norm(vector) == 0:
        return np.nan
    cos_theta = np.dot(vector, unit_diagonal) / np.linalg.norm(vector)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))


def project_onto_diagonal(x_vec, d_unit):
    """Project vector x onto the space diagonal defined by unit vector d_unit."""
    proj_scalar = np.dot(x_vec, d_unit)
    return proj_scalar * d_unit


def get_orthogonal_component(x_vec, d_unit):
    """Get the component of vector x orthogonal to the space diagonal."""
    return x_vec - project_onto_diagonal(x_vec, d_unit)


def angle_between_vectors(a, b):
    """Calculate the angle between two vectors in degrees."""
    dot = np.dot(a, b)
    norm_product = np.linalg.norm(a) * np.linalg.norm(b)
    if norm_product == 0:
        return np.nan
    cos_theta = np.clip(dot / norm_product, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))


def normalized_change_vector(a, b):
    """
    Calculate normalized change between two vectors.
    Returns element-wise absolute differences divided by maximum difference.
    """
    #a = np.asarray(a)
    #b = np.asarray(b)
    #delta = b - a
    #abs_delta = np.abs(delta)
    #max_val = np.max(abs_delta)
    #return np.zeros_like(delta) if max_val == 0 else abs_delta / max_val
    #-------------------------------------------------
    a = np.asarray(a)
    b = np.asarray(b)  #delta = b - a
    #len_o = np.linalg.norm(a)
    #len_p = np.linalg.norm(b)
    #change_vector = delta / (len_o + len_p)
    #return abs(change_vector)
    assert a.shape[0] == 13 and b.shape[0] == 13


    shiftvector = a-b
    norm = np.linalg.norm(shiftvector)
    if norm == 0:
        return shiftvector
    return shiftvector / norm




# ======================
# MAIN PROCESSING
# ======================

def main():
    # Load input files
    df = pd.read_csv(input_file, sep="\t")
    expr_df = pd.read_csv(centroid_file, sep="\t")
    pairwise_df = pd.read_csv(pairwise_file, sep="\t")

    # Identify expression columns for genes (between gene_id and Orthogroup)
    gene_expr_start = df.columns.get_loc("gene_id") + 1
    gene_expr_end = df.columns.get_loc("Orthogroup")
    gene_expr_columns = df.columns[gene_expr_start:gene_expr_end]

    # Unit vector for space diagonal
    unit_diagonal = np.ones(len(gene_expr_columns)) / np.sqrt(len(gene_expr_columns))
    print(unit_diagonal)

    # 1. Calculate basic angle to space diagonal for genes if not already present
    if "angle" not in df.columns:
        df["angle"] = df[gene_expr_columns].apply(
            lambda row: calculate_angle(row.values.astype(float), unit_diagonal),
            axis=1
        )

    # 2. Calculate angle for centroids
    expr_columns = expr_df.columns[1:]  # All columns except 'Orthogroup'
    print(expr_columns)
    unit_diagonal_expr = np.ones(len(expr_columns)) / np.sqrt(len(expr_columns))
    print(unit_diagonal_expr)
    expr_df["angle_diagonal"] = expr_df[expr_columns].apply(
        lambda row: calculate_angle(row.values.astype(float), unit_diagonal_expr),
        axis=1
    )

    # 3. Calculate orthogonal components
    expr_df["orth_component"] = expr_df[expr_columns].apply(
        lambda row: get_orthogonal_component(row.values.astype(float), unit_diagonal_expr),
        axis=1
    )

    df["orth_component"] = df[gene_expr_columns].apply(
        lambda row: get_orthogonal_component(row.values.astype(float), unit_diagonal_expr),
        axis=1
    )

    # 4. Calculate angles between genes and their centroids
    centroid_orth_lookup = expr_df.set_index("Orthogroup")["orth_component"]

    df["orth_angle_to_centroid_orth"] = df.apply(
        lambda row: angle_between_vectors(
            row["orth_component"],
            centroid_orth_lookup.get(row["Orthogroup"], np.nan)
        ),
        axis=1
    )

    # 5. Process pairwise data
    orth_component_lookup = df.set_index("gene_id")["orth_component"]

    # Add orthogonal angles between gene pairs
    pairwise_df["orth_angle"] = pairwise_df.apply(
        lambda row: angle_between_vectors(
            orth_component_lookup.get(row["Gene1"], None),
            orth_component_lookup.get(row["Gene2"], None)
        ),
        axis=1
    )

    # Add normalized change vectors centroid always first argument for centroid - gene
    df["orth_change_to_centroid_normalized"] = df.apply(
        lambda row: normalized_change_vector(
            centroid_orth_lookup.get(row["Orthogroup"], None),
            row["orth_component"]
        ),
        axis=1
    )
    # ortholog is always the first argument that results in ortholog - paralog
    pairwise_df["orth_relative_change_normalized"] = pairwise_df.apply(
        lambda row: normalized_change_vector(
            orth_component_lookup.get(row["Gene1"], None),
            orth_component_lookup.get(row["Gene2"], None)
        ) if row["gene1_type"] == "ortholog" else
        normalized_change_vector(
            orth_component_lookup.get(row["Gene2"], None),
            orth_component_lookup.get(row["Gene1"], None)
        ),
        axis=1
    )

    # Add angle information to pairwise data
    angle_lookup = df.set_index("gene_id")["angle"]
    pairwise_df["angle_gene1"] = pairwise_df["Gene1"].map(angle_lookup)
    pairwise_df["angle_gene2"] = pairwise_df["Gene2"].map(angle_lookup)

    # Calculate angle differences based on gene types
    def compute_angle_difference(row):
        a1, a2 = row["angle_gene1"], row["angle_gene2"]
        t1, t2 = row["gene1_type"], row["gene2_type"]

        if t1 == "ortholog" and t2 == "ortholog":
            return abs(a1 - a2)
        elif t1 == "ortholog" and t2 == "paralog":
            return a1 - a2
        elif t1 == "paralog" and t2 == "ortholog":
            return a2 - a1
        return np.nan

    pairwise_df["angle_difference"] = pairwise_df.apply(compute_angle_difference, axis=1)

    # Calculate angle difference between genes and their centroids
    centroid_angle_lookup = expr_df.set_index("Orthogroup")["angle_diagonal"]
    df["angle_difference_to_diagonal_to_centroid"] = df.apply(
        lambda row: centroid_angle_lookup.get(row["Orthogroup"], np.nan) - row["angle"],
        axis=1
    )

    # Save all results
    expr_df.to_csv(output_expr_file, sep="\t", index=False)
    pairwise_df.to_csv(output_pairwise_file, sep="\t", index=False)
    df.to_csv(output_df_file, sep="\t", index=False)

    print(f"✔ Saved centroid angles to: {output_expr_file}")
    print(f"✔ Saved pairwise distances with gene angles to: {output_pairwise_file}")
    print(f"✔ Saved gene vs. centroid angle differences to: {output_df_file}")


if __name__ == "__main__":
    main()
import pandas as pd
import numpy as np
import pickle

def get_harmonic_mean(df):
    """
    Berechnet den harmonischen Mittelwert für BLAST-Daten und speichert zwei Dictionaries.
    :param df: Pandas DataFrame mit den BLAST-Daten.
    :return: Tuple aus zwei Dictionaries:
             1. Dictionary mit (qid, sid) als Schlüssel und harmonischem Mittelwert als Wert.
             2. Dictionary mit (query_protein_prefix, target_protein_prefix) als Schlüssel und harmonischem Mittelwert als Wert.
    """
    results_by_id = {}
    results_by_prefix = {}

    # Extrahiere die benötigten Spalten
    qid = df[0]
    sid = df[1]
    qstart = df[2].to_numpy()
    qend = df[3].to_numpy()
    qlen = df[4].to_numpy()
    sstart = df[5].to_numpy()
    send = df[6].to_numpy()
    slen = df[7].to_numpy()
    pident = df[8].to_numpy()
    query_protein_prefix = df[12]
    target_protein_prefix = df[13]

    # Berechnung des Overlaps
    denominator = qlen + slen
    if np.any(denominator == 0):
        raise ValueError("Division durch Null im Overlap-Berechnungsprozess!")

    overlap = ((qend - qstart) + (send - sstart)) / denominator * 100

    # Berechnung des harmonischen Mittelwerts
    harmonic = 2 * (pident * overlap) / (pident + overlap)

    # Ergebnisse in den Dictionaries speichern
    for q, s, h, qp, tp in zip(qid, sid, harmonic, query_protein_prefix, target_protein_prefix):
        # Dictionary mit (qid, sid) als Schlüssel
        results_by_id[(q, s)] = h

        # Dictionary mit (query_protein_prefix, target_protein_prefix) als Schlüssel
        results_by_prefix[(qp, tp)] = h

    return results_by_id, results_by_prefix

# Hauptprogramm
if __name__ == "__main__":
    # BLAST-Datei einlesen
    blast_file = "results/blast_results_mapped.txt"  # Pfad zur BLAST-Datei
    output_pickle_by_id = "results/harmonic_means_by_id.pkl"  # Pickle-Datei für (qid, sid)
    output_pickle_by_prefix = "results/harmonic_means_by_prefix.pkl"  # Pickle-Datei für (query_protein_prefix, target_protein_prefix)

    try:
        # DataFrame aus BLAST-Datei erstellen
        df = pd.read_csv(blast_file, header=None, sep='\s+')

        # Harmonische Mittelwerte berechnen
        results_by_id, results_by_prefix = get_harmonic_mean(df)

        # Beide Dictionaries separat als Pickle-Dateien speichern
        with open(output_pickle_by_id, "wb") as f:
            pickle.dump(results_by_id, f)

        with open(output_pickle_by_prefix, "wb") as f:
            pickle.dump(results_by_prefix, f)

        print(f"Dictionaries erfolgreich gespeichert:")
        print(f"  - Ergebnisse nach ID in '{output_pickle_by_id}'")
        print(f"  - Ergebnisse nach Prefix in '{output_pickle_by_prefix}'")
    except Exception as e:
        print(f"Fehler: {e}")

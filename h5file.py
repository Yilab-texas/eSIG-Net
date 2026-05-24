import math
import numpy as np
from itertools import product
from Bio import SeqIO
import h5py

def aac_feature(seq):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    seq = seq.upper()
    length = len(seq)
    features = [seq.count(aa) / length if length > 0 else 0 for aa in amino_acids]
    return features

def get_conjoint_triads():
    groups = {
        'A': '0', 'G': '0', 'V': '0',
        'I': '1', 'L': '1', 'F': '1', 'P': '1',
        'Y': '2', 'M': '2', 'T': '2', 'S': '2',
        'H': '3', 'N': '3', 'Q': '3', 'W': '3',
        'R': '4', 'K': '4',
        'D': '5', 'E': '5',
        'C': '6'
    }
    return groups

def ct_feature(seq):
    seq = seq.upper()
    groups = get_conjoint_triads()
    triads = [''.join(p) for p in product('0123456', repeat=3)]
    triad_dict = {k: 0 for k in triads}

    encoded = [groups.get(res, 'X') for res in seq]
    encoded = [e for e in encoded if e != 'X']
    for i in range(len(encoded) - 2):
        triad = ''.join(encoded[i:i+3])
        if triad in triad_dict:
            triad_dict[triad] += 1
    total = sum(triad_dict.values())
    return [triad_dict[k] / total if total > 0 else 0 for k in triads]

PCPNS = ['H1', 'H2', 'NCI', 'P1', 'P2', 'SASA', 'V']
AAPCPVS = {
    'A': {'H1': 0.62, 'H2': -0.5, 'NCI': 0.007187, 'P1': 8.1, 'P2': 0.046, 'SASA': 1.181, 'V': 27.5},
    'C': {'H1': 0.29, 'H2': -1.0, 'NCI': -0.036610, 'P1': 5.5, 'P2': 0.128, 'SASA': 1.461, 'V': 44.6},
    'D': {'H1': -0.90, 'H2': 3.0, 'NCI': -0.023820, 'P1': 13.0, 'P2': 0.105, 'SASA': 1.587, 'V': 40.0},
    'E': {'H1': 0.74, 'H2': 3.0, 'NCI': 0.006802, 'P1': 12.3, 'P2': 0.151, 'SASA': 1.862, 'V': 62.0},
    'F': {'H1': 1.19, 'H2': -2.5, 'NCI': 0.037552, 'P1': 5.2, 'P2': 0.290, 'SASA': 2.228, 'V': 115.5},
    'G': {'H1': 0.48, 'H2': 0.0, 'NCI': 0.179052, 'P1': 9.0, 'P2': 0.000, 'SASA': 0.881, 'V': 0.0},
    'H': {'H1': -0.40, 'H2': -0.5, 'NCI': -0.010690, 'P1': 10.4, 'P2': 0.230, 'SASA': 2.025, 'V': 79.0},
    'I': {'H1': 1.38, 'H2': -1.8, 'NCI': 0.021631, 'P1': 5.2, 'P2': 0.186, 'SASA': 1.810, 'V': 93.5},
    'K': {'H1': -1.50, 'H2': 3.0, 'NCI': 0.017708, 'P1': 11.3, 'P2': 0.219, 'SASA': 2.258, 'V': 100.0},
    'L': {'H1': 1.06, 'H2': -1.8, 'NCI': 0.051672, 'P1': 4.9, 'P2': 0.186, 'SASA': 1.931, 'V': 93.5},
    'M': {'H1': 0.64, 'H2': -1.3, 'NCI': 0.002683, 'P1': 5.7, 'P2': 0.221, 'SASA': 2.034, 'V': 94.1},
    'N': {'H1': -0.78, 'H2': 2.0, 'NCI': 0.005392, 'P1': 11.6, 'P2': 0.134, 'SASA': 1.655, 'V': 58.7},
    'P': {'H1': 0.12, 'H2': 0.0, 'NCI': 0.239531, 'P1': 8.0, 'P2': 0.131, 'SASA': 1.468, 'V': 41.9},
    'Q': {'H1': -0.85, 'H2': 0.2, 'NCI': 0.049211, 'P1': 10.5, 'P2': 0.180, 'SASA': 1.932, 'V': 80.7},
    'R': {'H1': -2.53, 'H2': 3.0, 'NCI': 0.043587, 'P1': 10.5, 'P2': 0.291, 'SASA': 2.560, 'V': 105.0},
    'S': {'H1': -0.18, 'H2': 0.3, 'NCI': 0.004627, 'P1': 9.2, 'P2': 0.062, 'SASA': 1.298, 'V': 29.3},
    'T': {'H1': -0.05, 'H2': -0.4, 'NCI': 0.003352, 'P1': 8.6, 'P2': 0.108, 'SASA': 1.525, 'V': 51.3},
    'V': {'H1': 1.08, 'H2': -1.5, 'NCI': 0.057004, 'P1': 5.9, 'P2': 0.140, 'SASA': 1.645, 'V': 71.5},
    'W': {'H1': 0.81, 'H2': -3.4, 'NCI': 0.037977, 'P1': 5.4, 'P2': 0.409, 'SASA': 2.663, 'V': 145.5},
    'Y': {'H1': 0.26, 'H2': -2.3, 'NCI': 117.3, 'P1': 6.2, 'P2': 0.298, 'SASA': 2.368, 'V': 0.023599},
}

def avg_sd(values):
    avg = sum(values) / len(values)
    sd = math.sqrt(sum((v - avg) ** 2 for v in values) / len(values))
    return avg, sd

def normalize_pcp():
    PCPVS = {pcp: [v[pcp] for v in AAPCPVS.values()] for pcp in PCPNS}
    avg_sd_dict = {pcp: avg_sd(vals) for pcp, vals in PCPVS.items()}
    norm = {}
    for aa, props in AAPCPVS.items():
        norm[aa] = {pcp: (props[pcp] - avg_sd_dict[pcp][0]) / avg_sd_dict[pcp][1] for pcp in PCPNS}
    return norm

NORMALIZED = normalize_pcp()

def ac_feature(seq, lag=30):
    seq = [aa for aa in seq.upper() if aa in NORMALIZED]
    L = len(seq)
    if L < lag + 1:
        return [0.0] * (len(PCPNS) * lag)
    vec = []
    for pcp in PCPNS:
        values = [NORMALIZED[aa][pcp] for aa in seq]
        mean = sum(values) / L
        values = [v - mean for v in values]
        for d in range(1, lag + 1):
            ac = sum(values[i] * values[i + d] for i in range(L - d)) / (L - d)
            vec.append(ac)
    return vec

def extract_573_features(seq):
    return aac_feature(seq) + ct_feature(seq) + ac_feature(seq)

def fasta_to_feature_dict(fasta_path, output_h5_path):
    feature_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        pid = record.id
        seq = str(record.seq)
        try:
            vec = extract_573_features(seq)
            if len(vec) == 573:
                feature_dict[pid] = np.array(vec, dtype=np.float64)
        except Exception as e:
            print(f"Error on {pid}: {e}")
    with h5py.File(output_h5_path, "w") as f:
        for k, v in feature_dict.items():
            f.create_dataset(k, data=v)
    print(f"Saved {len(feature_dict)} proteins to {output_h5_path}")

if __name__ == "__main__":
    fasta_file = "combin3_dedup.fasta"
    output_h5 = "combin3_dedup.h5"
    fasta_to_feature_dict(fasta_file, output_h5)


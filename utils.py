import numpy as np
from scipy.spatial.distance import pdist, squareform
from Bio.PDB import *
from matplotlib import pyplot as plt


# PDB conventions: https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/COMPOUND_NAME
nonstandard_aa_substitutions = {
    "MSE": "M",
    "ADP": "A",
    "HG": "X",
    "ACE": "X",
    "ATP": "X",
    "CA": "X",
    "CLR": "X",
    "HOH": "X",
    "MG": "X",
    "NH2": "X",
    "NI": "X",
    "PNM": "X",
    "PO4": "X",
    "SO4": "X",
    "ZN": "X",
    "CL": "X",
    "ANP": "X",
    "GOL": "X",
    "CXT": "X",
    "CDL": "X",
    "PC1": "X",
    "LDM": "X",
    "BOG": "X",
    "GDM": "X",
    "FOS": "X",
    "V4O": "X",
    "FUC": "X",
    "NGA": "X",
    "B46": "X",
    "ACT": "X",
    "EDO": "X",
    "AF3": "X",
    "MK7": "X"
}


def to_one_letter_seq(chain):
    seq = ""
    for residue in chain.get_residues():
        three_letter_name = residue.get_resname()
        try:
            seq += Polypeptide.three_to_one(three_letter_name)
        except KeyError:
            if three_letter_name in nonstandard_aa_substitutions:
                seq += nonstandard_aa_substitutions[three_letter_name]
            else:
                print(f"Non-standard amino acid {three_letter_name}, placed an X")
                seq += "X"
    return seq


def indices_in_ref_and_query(alignment):
    i_ref = alignment.r_pos
    i_query = alignment.q_pos
    idxs_ref = []
    idxs_query = []
    for pair in alignment.cigar:
        if pair[1] == "M":
            j_ref = i_ref + pair[0]
            idxs_ref += list(range(i_ref, j_ref))
            i_ref = j_ref
            j_query = i_query + pair[0]
            idxs_query += list(range(i_query, j_query))
            i_query = j_query
        elif pair[1] == "D":
            i_ref += pair[0]

        elif pair[1] == "I":
            i_query += pair[0]

    return idxs_ref, idxs_query


def calc_CA_dist_matrix(chain, idx_subset):
    """Returns a matrix of C-alpha distances in a (subset of a) chain"""
    idx_subset = set(idx_subset)
    residue_coords = [residue["CA"].coord for i, residue in enumerate(chain.get_residues()) if i in idx_subset]

    return squareform(pdist(residue_coords))


def calc_min_dist_matrix(chain, idx_subset):
    """Returns a matrix of minimum distances between residues in a (subset of a) chain"""
    idx_subset = set(idx_subset)

    return squareform(
        np.array([min([atom_in_res_i - atom_in_res_j for atom_in_res_i in res_i for atom_in_res_j in res_j])
                 for i, res_i in enumerate(chain.get_residues()) if i in idx_subset
                 for j, res_j in enumerate(chain.get_residues()) if j > i and j in idx_subset])
    )


def top_n_contact_mask(scores_matrix, sequence_proximity_mask, n, return_flattened_scores=False):
    sq = squareform(scores_matrix, checks=False)
    sq_proximity = np.logical_not(squareform(sequence_proximity_mask, checks=False))
    sq[sq_proximity] = -np.inf
    if return_flattened_scores:
        return sq
    contact_mask = np.zeros(len(sq), dtype=bool)
    argsrt = np.argsort(sq)[::-1]
    contact_mask[argsrt[:n]] = True
    
    return squareform(contact_mask)


def contact_matrix_comparison(dist_matrix,
                              scores_matrix,
                              max_eucl_dist=8,
                              n_pred=None,
                              min_sequence_dist=5,
                              contact_mask=None,
                              return_data_for_roc=False):
    assert dist_matrix.shape == scores_matrix.shape
    n_residues = dist_matrix.shape[0]
    sequence_proximity_mask = np.abs(
        np.arange(n_residues) - np.arange(n_residues)[:, None]
    ) >= min_sequence_dist
    
    if contact_mask is None:
        eucl_dist_mask = dist_matrix < max_eucl_dist
        contact_mask = np.triu(np.logical_and(eucl_dist_mask, sequence_proximity_mask))        

    if n_pred is None:
        n_contacts = int(np.sum(contact_mask))
        n_pred = n_contacts

    contact_mask_pred = top_n_contact_mask(scores_matrix,
                                           sequence_proximity_mask,
                                           n=n_pred,
                                           return_flattened_scores=return_data_for_roc)
    if return_data_for_roc:
        contact_scores = contact_mask_pred
        return contact_scores, squareform(contact_mask, checks=False)

    contact_mask_pred = np.tril(contact_mask_pred)

    contact_mask_true_pos = np.logical_and(contact_mask_pred.T, contact_mask).T
    contact_mask_false_pos = np.logical_xor(contact_mask_true_pos, contact_mask_pred)
    contact_mask_false_neg = np.logical_and(np.logical_not(contact_mask_pred).T, contact_mask)

    full_matrix = contact_mask_true_pos.T * 1. + contact_mask_false_neg * 3 + contact_mask_true_pos * 1. + contact_mask_false_pos * 4.
    full_matrix[full_matrix == 0] = np.nan
    
    return full_matrix


def plot_contact_matrices(dist_matrix,
                          scores_matrix,
                          max_eucl_dist=8,
                          n_pred=None,
                          min_sequence_dist=5,
                          title=None):
    full_matrix = contact_matrix_comparison(dist_matrix,
                                            scores_matrix,
                                            max_eucl_dist=max_eucl_dist,
                                            n_pred=n_pred,
                                            min_sequence_dist=min_sequence_dist)

    plt.figure(figsize=(10, 10))
    plt.matshow(full_matrix, fignum=1, cmap="viridis")

    plt.colorbar()
    if title is not None:
        plt.title(title)

    plt.show()


def read_bmDCA_parameters(path):
    """Read in a bmDCA-generated text file containing learnt Potts model parameters, and return
    the field and coupling tensors."""
    J2_idxs = []
    J2_vals = []
    h_idxs = []
    h_vals = []
    with open(path, "r") as f:
        for line in f:
            l = line.rstrip("\n").split(" ")
            val = float(l[-1])
            idxs = list(map(int, l[1:-1]))
            if l[0] == "J":
                J2_idxs.append(idxs)
                J2_vals.append(val)
            elif l[0] == "h":
                h_idxs.append(idxs)
                h_vals.append(val)

    J2_idxs = np.asarray(J2_idxs)
    J2_vals = np.asarray(J2_vals, dtype=np.float64)
    h_idxs = np.asarray(h_idxs)
    h_vals = np.asarray(h_vals, dtype=np.float64)

    n_sites = np.max(h_idxs[:, 0]) + 1
    n_states = np.max(h_idxs[:, 1]) + 1

    h = np.zeros((n_sites, n_states), dtype=np.float64)
    J2 = np.zeros((n_sites, n_sites, n_states, n_states), dtype=np.float64)

    h[tuple(h_idxs.T)] = h_vals
    J2[tuple(J2_idxs.T)] = J2_vals

    # Symmetrize J2
    J2 += np.transpose(J2, (1, 0, 3, 2))

    return h, J2


def zero_sum_gauge_frob_scores(J2, apc=True):
    """Compute a score for contacts between sites by first passing to a zero-sum gauge, then computing a Frobenius
    norm, and finally applying the average product correction (APC) if desired."""
    # Pass to zero-sum gauge
    J2_zs = J2
    J2_zs -= np.mean(J2, axis=2, keepdims=True)
    J2_zs -= np.mean(J2, axis=3, keepdims=True)
    J2_zs += np.mean(J2, axis=(2, 3), keepdims=True)

    # Frobenius norm
    S_frob = np.linalg.norm(J2_zs, axis=(2, 3), ord='fro')

    # Average-product correction
    S = S_frob
    if apc:
        S -= (np.mean(S_frob, axis=1, keepdims=True) * np.mean(S_frob, axis=0, keepdims=True)) / np.mean(S_frob)

    return S

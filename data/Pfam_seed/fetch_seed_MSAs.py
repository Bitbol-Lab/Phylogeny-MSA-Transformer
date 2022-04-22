"""Script to fetch Pfam seed MSAs and align them to their corresponding Pfam HMMs."""

import pathlib
import os
from subprocess import run
import requests

import tqdm

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from prody.database.pfam import fetchPfamMSA


pfam_families = [
    "PF00004",
    "PF00005",
    "PF00041",
    "PF00072",
    "PF00076",
    "PF00096",
    "PF00153",
    "PF00271",
    "PF00397",
    "PF00512",
    "PF00595",
    "PF01535",
    "PF02518",
    "PF07679",
    "PF13354"
]

fasta_folder = pathlib.Path("msa")
hmm_folder = pathlib.Path("hmm")

fmt = "fasta"
alignment = "seed"
gaps = None


def _sequences_to_str_array(filename, fmt="fasta"):
    return np.asarray([list(str(record.seq)) for record in SeqIO.parse(filename, fmt)])


if __name__ == "__main__":
    for pfam_family in tqdm.tqdm(pfam_families):
        # 1) Download seed MSAs and PFAM HMMs
        fetchPfamMSA(pfam_family,
                     format=fmt,
                     alignment=alignment,
                     order="tree",
                     inserts="lower",
                     gaps=gaps,
                     folder=fasta_folder,
                     outname=pfam_family,
                     timeout=1000)
        hmm_url = f"https://pfam.xfam.org/family/{pfam_family}/hmm"
        hmm_filename = hmm_folder / f"{pfam_family}.hmm"
        r = requests.get(hmm_url, allow_redirects=True)
        with open(hmm_filename, "wb") as f:
            f.write(r.content)

        # 2) Align seed MSAs to downloaded HMMs with hmmalign, producing Stockholm file
        msa_filename = fasta_folder / f"{pfam_family}_{alignment}.{fmt}"
        aligned_msa_filename_stockholm = fasta_folder / f"{pfam_family}_{alignment}_hmmalign.stockholm"
        run(["hmmalign", "--amino", "-o", aligned_msa_filename_stockholm, hmm_filename, msa_filename])

        # 3) Convert Stockholm aligned MSA to FASTA format
        aligned_msa_filename_fasta = fasta_folder / f"{pfam_family}_{alignment}_hmmalign.fasta"
        parsed = list(SeqIO.parse(aligned_msa_filename_stockholm, "stockholm"))
        with open(aligned_msa_filename_fasta, "w") as output_handle:
            SeqIO.write(parsed, output_handle, "fasta")

        # 4) Keep only match and deletion states
        msa_arr = _sequences_to_str_array(aligned_msa_filename_fasta)
        valid_idxs = []
        for col in range(msa_arr.shape[1]):
            unique_set = set(np.unique(msa_arr[:, col]))
            if not set(string.ascii_lowercase).intersection(unique_set):
                valid_idxs.append(col)

        no_inserts_filename = fasta_folder / f"{pfam_family}_{alignment}_hmmalign_no_inserts.fasta"
        parsed = list(SeqIO.parse(aligned_msa_filename_fasta, fmt))
        with open(no_inserts_filename, "w") as output_handle:
            for i, record in enumerate(parsed):
                record.seq = Seq("".join(msa_arr[i, valid_idxs]))
            SeqIO.write(parsed, output_handle, fmt)

        os.remove(aligned_msa_filename_fasta)
        os.remove(aligned_msa_filename_stockholm)

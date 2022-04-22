# Phylogeny-MSA-Transformer
Supporting repository for "Protein language models trained on multiple sequence alignments learn phylogenetic relationships" (preprint: https://doi.org/10.1101/2022.03.29.486219)

## Getting started

Clone this repository on your local machine by running
```bash
git clone git@github.com:Bitbol-Lab/Phylogeny-MSA-Transformer.git
```
and move inside the root folder.
We recommend creating and activating a dedicated ``conda`` or ``virtualenv`` Python virtual environment.
Then, install the required libraries:
```python
python -m pip install -U -r requirements.txt
```

## Requirements
In order to run the notebooks, the following python packages are required:

- tqdm
- jupyter
- matplotlib
- statsmodels
- biopython
- esm==0.4.0

``prody`` and ``HMMER`` are required to run the Python script ``data/Pfam_Seed/fetch_seed_MSAs.py``, if you wish to create new
Pfam full MSAs instead of using the ones provided.

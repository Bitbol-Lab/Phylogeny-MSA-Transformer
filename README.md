# Phylogeny-MSA-Transformer

[![DOI](https://zenodo.org/badge/483996183.svg)](https://zenodo.org/badge/latestdoi/483996183)

Supporting repository for ["Protein language models trained on multiple sequence alignments learn phylogenetic relationships" (Lupo, Sgarbossa, and Bitbol, 2022)](https://www.nature.com/articles/s41467-022-34032-y). The MSA Transformer model used here was introduced in [(Rao el al, 2021)](https://proceedings.mlr.press/v139/rao21a.html).

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
- swalign
- esm==0.4.0

``prody`` and ``HMMER`` are required to run the Python script ``data/Pfam_Seed/fetch_seed_MSAs.py``, if you wish to create new
Pfam full MSAs instead of using the ones provided.

### Citation

Our work can be cited using the following BibTeX entry:

```bibtex
@article{lupo2022protein,
  title={Protein language models trained on multiple sequence alignments learn phylogenetic relationships},
  author={Lupo, Umberto and Sgarbossa, Damiano and Bitbol, Anne-Florence},
  year={2022},
  volume={13},
  number={6298},
  journal={Nat. Commun.},
  doi={10.1038/s41467-022-34032-y}
}
```

#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#
#  All functions which prepare files for the meme analysis

import pandas as pd

from helpers.common import get_de_genes, extract_utr_sequence
from helpers.constants import *


def prepare_fasta(utr, filename=None):
    genes = get_de_genes()
    t2g = pd.read_csv(FILE_MAPPING)
    t2g = dict(zip(t2g["Transcript stable ID"], t2g["Gene stable ID"]))
    seqs = extract_utr_sequence(utr)
    seqs = {t2g[k]: v for k, v in seqs.items() if t2g[k] in genes}
    if filename is None:
        filename = f"{utr}utr.fasta"
    with open(filename, "w") as f:

        for k, v in seqs.items():
            if len(v) >= 8:  # This is condition for meme
                print(f">{k}", file=f)
                print(v, file=f)

    return filename


def run():
    prepare_fasta(5)

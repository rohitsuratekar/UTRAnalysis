#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#
#  UTR Analysis and related statistics
#
#  Common functions

import numpy as np
import pandas as pd
from helpers.constants import *
from collections import defaultdict


def get_de_genes() -> np.ndarray:
    df = pd.read_csv(FILE_DIFF)
    df = df[df["padj"] <= P_VALUE_CUTOFF].reset_index(drop=True)
    df = df[(df["log2FoldChange"] <= LOG_2_FOLD_CUTOFF * -1) | (
            df["log2FoldChange"] >= LOG_2_FOLD_CUTOFF)]
    return df["gene_id"].values


def extract_utr_sequence(utr: int) -> dict:
    if utr not in [5, 3]:
        raise Exception(f"UTR should be either 5 or 3. You have provided '"
                        f"{utr}'")
    filename = FILE_5UTR if utr == 5 else FILE_3UTR

    with open(filename) as f:
        lines = f.readlines()
        seq = defaultdict(list)
        current_id = None
        for line in lines:
            if line.strip().startswith(">"):
                current_id = line.strip().split("|")[2]
                continue
            seq[current_id].append(line.strip())

    seq = {k: "".join(v) for k, v in seq.items() if
           v[0] != "Sequence unavailable"}
    return seq

#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#
#  UTR Analysis and related statistics
#
#  Common functions

from collections import defaultdict
import pandas as pd
from helpers.constants import *


def get_genes(*, direction, condition, atlas):
    condition = condition.strip().replace(" ", "_")
    filename = f"mito/{direction}_{atlas}_{condition}.csv"
    df = pd.read_csv(filename)
    if "FDR" in df.columns:
        df = df[df["FDR"] <= 0.05]
    if "log2FC" in df.columns:
        if direction == "up":
            df = df[df["log2FC"] > 0]
        else:
            df = df[df["log2FC"] < 0]
    return df["gene_id"].to_numpy()


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

#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#
#  All functions which prepare files for the meme analysis

import os
import subprocess

import pandas as pd

from helpers.common import extract_utr_sequence, get_genes
from helpers.constants import *


def prepare_fasta(utr, genes, filename=None):
    t2g = pd.read_csv(FILE_MAPPING)
    t2g = dict(zip(t2g["Transcript stable ID"], t2g["Gene stable ID"]))
    seqs = extract_utr_sequence(utr)
    if genes is None:
        seqs = {t2g[k]: v for k, v in seqs.items()}
    else:
        seqs = {t2g[k]: v for k, v in seqs.items() if t2g[k] in genes}
    if filename is None:
        filename = f"{utr}utr.fasta"
    with open(filename, "w") as f:
        for k, v in seqs.items():
            if len(v) >= 6:
                print(f">{k}", file=f)
                print(v, file=f)

    return filename


def generate_all():
    folder = "fasta"
    if not os.path.isdir(folder):
        os.mkdir(folder)
    direction = "up"
    atlas = "mitocarta"
    conditions = ["hsf wt", "hsf mia40", "mars wt", "mars mia40"]
    for condition in conditions:
        condition = "_".join(condition.strip().split(" "))
        filename = f"{folder}/{direction}_{atlas}_{condition}_3utr.fasta"
        genes = get_genes(direction=direction,
                          atlas=atlas,
                          condition=condition)
        prepare_fasta(3, genes, filename=filename)


def run_streme():
    meme_path = ":/home/dex/Softwares/meme/bin"
    meme_path += ":/home/dex/Softwares/meme/libexec/meme-5.2.0"
    targets = [
        "up_mitocarta_mars_mia40_3utr.fasta",
        "up_mitocarta_mars_wt_3utr.fasta",
        "up_mitocarta_hsf_mia40_3utr.fasta",
        "up_mitocarta_hsf_wt_3utr.fasta"
    ]
    control = "./fasta/control.fasta"
    env = os.environ.copy()
    env['PATH'] += meme_path
    for i in range(len(targets)):
        name = targets[i].replace("_3utr.fasta", '')
        output_path = f"./streme/{name}"
        primary = f"./fasta/{targets[i]}"

        opts_random = [
            "streme",
            "--o", output_path,
            "--objfun", "de",
            "--pvt", "0.05",
            "--patience", "1",
            "--minw", "6",
            "--maxw", "12",
        ]
        opts_all = [x for x in opts_random]
        opts_all.extend([
            "--n", control
        ])
        opts_all.extend(["--p", primary])
        opts_random.extend(["--p", primary])
        # For random
        subprocess.run(opts_random, env=env)
        subprocess.run(opts_all, env=env)


def run():
    run_streme()

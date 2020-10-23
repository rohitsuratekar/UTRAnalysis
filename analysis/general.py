#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from scipy.stats import gaussian_kde
from helpers.constants import *
from helpers.common import get_de_genes, extract_utr_sequence

p = Palette()


def filter_utrs(genes, seq) -> list:
    df = pd.read_csv(FILE_MAPPING)
    df = df.sort_values(by="Transcript length (including UTRs and CDS)",
                        ascending=False)
    seq = pd.DataFrame.from_dict(seq, orient="index")
    seq = seq.rename(columns={0: "sequence"})
    df = df.set_index("Transcript stable ID")
    df = pd.concat([df, seq], axis=1)
    del seq
    df = df[df["Gene stable ID"].isin(genes)]
    df = df[df['sequence'].notna()]
    df = df.sort_values(by=["Transcript length (including UTRs "
                            "and CDS)"], ascending=False)
    df = df.drop_duplicates(subset=["Gene stable ID"])
    df = dict(zip(df["Gene name"], df["sequence"]))
    df = [(k, v) for k, v in df.items()]
    return df


def plot_utr_length():
    x_lim = 2000
    genes = get_de_genes()
    for i in [5, 3]:
        raw = extract_utr_sequence(i)
        all_utr = [len(x) for x in raw.values()]
        utr = filter_utrs(genes, raw)
        lens = [len(x[1]) for x in utr]
        density_utr = gaussian_kde(all_utr)
        density_len = gaussian_kde(lens)
        xs = np.linspace(0, x_lim, 200)
        line_style = "--" if i == 3 else None
        plt.plot(xs, density_utr(xs), color=p.blue(),
                 lw=2, label=f"{i}' UTR Global", ls=line_style,
                 zorder=3)
        plt.plot(xs, density_len(xs), color=p.red(), ls=line_style,
                 lw=2, label=f"{i}' UTR DE",
                 zorder=3)
    plt.xlabel("UTR length")
    plt.ylabel("Frequency (density)")
    plt.legend(loc=0)
    plt.grid(zorder=0, ls=":")
    ax = plt.gca()
    ax.set_facecolor(p.gray(shade=15))
    plt.annotate(f"Lengths above {x_lim} are omitted",
                 (1, 1.05), ha="right", xycoords="axes fraction", va="top",
                 fontstyle="italic", color=p.gray(shade=70))
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)

    plt.show()


def run():
    plot_utr_length()

#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from scipy.stats import gaussian_kde
from scipy.stats import ttest_ind

from helpers.common import extract_utr_sequence, get_genes
from helpers.constants import *

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


def plot_utr_length_distribution():
    x_lim = 2000
    colors = [p.cyan, p.magenta]
    conditions = ["hsf mia40", "mars mia40"]
    atlas = "mitocarta"
    direction = "down"
    samples = []
    if atlas == "mitocarta":
        colors = [p.blue, p.red]
    total = []
    for i, condition in enumerate(conditions):
        genes = get_genes(direction=direction,
                          condition=condition,
                          atlas=atlas)
        raw = extract_utr_sequence(3)
        utr = filter_utrs(genes, raw)
        values = [len(x[1]) for x in utr]
        samples.append([x for x in values])
        mean_value = sum(values) / len(values)
        total.append(len(values))
        print(f"{condition.upper()} Total : {len(values)}")
        print(f"{condition.upper()} Average : {round(mean_value, 2)}")

        density_utr = gaussian_kde(values)
        xs = np.linspace(0, x_lim, 200)
        plt.plot(xs, density_utr(xs), color=colors[i](),
                 lw=2, label=f"{condition.upper()}",
                 zorder=3)
        plt.axvline(mean_value, color=colors[i](), ls="--")
        txt_loc = plt.gca().get_ylim()[1] / 2
        align = "right"
        if i == 0:
            align = "left"
        plt.annotate(f"{round(mean_value, 2)}",
                     xy=(mean_value, txt_loc),
                     rotation=90, ha=align, va="center",
                     bbox=dict(fc=p.white(), ec=colors[i]()))

    # Welch's unequal variances t-test
    tt = ttest_ind(a=samples[0], b=samples[1], equal_var=False)

    plt.xlabel("UTR length")
    plt.ylabel("Frequency (density)")
    plt.legend(loc=0)
    plt.grid(zorder=0, ls=":")
    ax = plt.gca()
    ax.set_facecolor(p.gray(shade=15))
    plt.annotate(f"p-value : {round(tt[1], 6)}"
                 f"\n\nLengths above {x_lim}\nare not shown"
                 f"\n\n n = {total[0]} ({conditions[0]}),"
                 f"\n {total[1]} ({conditions[1]})",
                 (0.96, 0.81), ha="right", xycoords="axes fraction", va="top",
                 fontstyle="italic", color=p.gray(shade=70))

    plt.title(f"3' UTR Length Distribution of {direction} regulated genes "
              f"({atlas})")
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)

    plt.show()


def run():
    plot_utr_length_distribution()

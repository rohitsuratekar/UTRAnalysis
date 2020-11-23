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


def common_gene_analysis():
    x_lim = 2000
    direction = "up"
    atlas = "general"
    wt = get_genes(direction=direction,
                   condition="mars wt",
                   atlas=atlas)
    mia40 = get_genes(direction=direction,
                      condition="mars mia40",
                      atlas=atlas)

    common = list(set(wt).intersection(set(mia40)))
    wt_only = list(set(wt).difference(set(mia40)))
    mia40_only = list(set(mia40).difference(set(wt)))

    conditions = [common, wt_only, mia40_only]
    labels = [f"common (n={len(common)})",
              f"only wt (n={len(wt_only)})",
              f"only mia40 (n={len(mia40_only)})"]
    colors = [p.gray, p.blue, p.red]
    samples = []

    for i in range(len(conditions)):
        raw = extract_utr_sequence(3)
        utr = filter_utrs(conditions[i], raw)
        values = [len(x[1]) for x in utr]
        samples.append([x for x in values])
        density_utr = gaussian_kde(values)
        xs = np.linspace(0, x_lim, 200)
        plt.plot(xs, density_utr(xs), color=colors[i](),
                 lw=2, label=f"{labels[i]}",
                 zorder=3)

    tt1 = ttest_ind(a=samples[0], b=samples[1], equal_var=False)
    tt2 = ttest_ind(a=samples[0], b=samples[2], equal_var=False)

    plt.xlabel("UTR length")
    plt.ylabel("Frequency (density)")
    plt.legend(loc=0)
    plt.grid(zorder=0, ls=":")
    ax = plt.gca()
    ax.set_facecolor(p.gray(shade=15))

    plt.annotate(f"p-value : "
                 f"\n common vs wt : {round(tt1[1], 30)}"
                 f"\n common vs mia40 : {round(tt2[1], 8)}"
                 f"\n\nLengths above {x_lim}\nare not shown",
                 (0.96, 0.78), ha="right", xycoords="axes fraction", va="top",
                 fontstyle="italic", color=p.gray(shade=70))

    plt.title(f"{direction}-regulated genes ({atlas})")
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)
    plt.show()


def run():
    common_gene_analysis()

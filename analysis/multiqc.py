#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#
#  All plots for MultiQC data

from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from SecretColors.cmaps import BrewerMap
from matplotlib.patches import Patch

p = Palette()

color_list = [p.red(), p.red(shade=20), p.blue(shade=40)]
cm = BrewerMap(matplotlib).from_list(color_list, is_qualitative=True)

FILE_FASTQC = "multiqcdata/fastqc-status-check-heatmap.csv"
FILE_SORTMERNA = "multiqcdata/sortmerna-detailed-plot.csv"
FILE_STAR = "multiqcdata/star_alignment_plot.csv"


def plot_fastqc():
    df = pd.read_csv(FILE_FASTQC)
    analysis_type = df["Section Name"].values
    samples = df["Series 1 (y)"].values
    values = df["Series 1 (value)"].values
    analysis = defaultdict(list)
    for i, value in enumerate(analysis_type):
        analysis[value].append((samples[i], values[i]))

    labels = None
    x_labels = []
    full_data = []
    for i, key in enumerate(analysis.keys()):
        x_labels.append(key)
        if labels is None:
            labels = sorted([x[0] for x in analysis[key]])
        data = sorted(analysis[key], key=lambda x: x[0])
        data = [x[1] for x in data]
        full_data.append(data)

    fig = plt.figure(figsize=(14, 8))  # type: plt.Figure
    ax = fig.add_subplot()  # type: plt.Axes
    full_data = np.asarray(full_data).transpose()
    ax.imshow(full_data, cmap=cm, aspect='auto')
    ax.set_yticks(range(0, len(labels)))
    ax.set_yticklabels(labels)
    x_labels = ["\n".join(str(x).strip().rsplit(" ", 2)) for x in x_labels]
    ax.set_xticks(range(0, len(x_labels)))
    ax.set_xticklabels(x_labels)
    for i in range(0, len(x_labels) - 1):
        ax.axvline(i + 0.5, color=p.white(), alpha=0.2)
    for i in range(0, len(labels) - 1):
        ax.axhline(i + 0.5, color=p.white(), alpha=0.2)
        if i % 3 == 0:
            ax.axhline(i - 0.5, color=p.white(), alpha=0.8, ls="--")
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300, transparent=True)
    plt.show()


def plot_sortmerna():
    c1 = "silva-euk-18s-id95_count"
    c2 = "silva-euk-28s-id98_count"
    df = pd.read_csv(FILE_SORTMERNA)
    df = df.sort_values(by="Category", ascending=False)
    df[c1] = df[c1] / 1000000
    df[c2] = df[c2] / 1000000
    labels = df["Category"].values
    c1_data = df[c1].values
    c2_data = df[c2].values

    fig = plt.figure(figsize=(11, 8))

    for i in range(0, len(labels)):
        plt.barh(i, c1_data[i], color=p.aqua(shade=40))
        plt.barh(i, c2_data[i], left=c1_data[i], color=p.aqua(shade=60))

    plt.yticks(range(0, len(labels)), labels)
    plt.xlabel("Reads Filtered (in Millions)")

    handles = [
        Patch(facecolor=p.aqua(shade=40), label=c1),
        Patch(facecolor=p.aqua(shade=60), label=c2)
    ]

    plt.legend(handles=handles, loc=0)
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300, transparent=True)
    plt.show()


def plot_star():
    df = pd.read_csv(FILE_STAR)
    df = df.sort_values(by="Category", ascending=False)
    c1 = "Uniquely mapped"
    c2 = "Mapped to multiple loci"
    c3 = "Mapped to too many loci"
    c4 = "Unmapped: too short"
    c5 = "Unmapped: other"
    labels = df["Category"].values
    df = df.set_index("Category")
    df = df.apply(lambda x: x / 1000000)
    fig = plt.figure(figsize=(12, 8))
    for i in range(len(labels)):
        d1 = df[c1].values[i]
        d2 = df[c2].values[i]
        d3 = df[c3].values[i]
        d4 = df[c4].values[i]
        d5 = df[c5].values[i]
        plt.barh(i, d1, color=p.cyan(shade=30))
        plt.barh(i, d2, left=d1, color=p.cyan(shade=40))
        plt.barh(i, d3, left=d1 + d2, color=p.cyan(shade=50))
        plt.barh(i, d4, left=d1 + d2 + d3, color=p.cyan(shade=60))
        plt.barh(i, d5, left=d1 + d2 + d3 + d4, color=p.cyan(shade=70))

    plt.yticks(range(len(labels)), labels)
    handles = [
        Patch(facecolor=p.cyan(shade=30), label=c1),
        Patch(facecolor=p.cyan(shade=40), label=c2),
        Patch(facecolor=p.cyan(shade=50), label=c3),
        Patch(facecolor=p.cyan(shade=60), label=c4),
        Patch(facecolor=p.cyan(shade=70), label=c5),
    ]

    plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5),
               fancybox=False, shadow=False)
    plt.grid(axis="x", color=p.black(), alpha=0.5, ls="--")
    plt.xlabel("Reads (in Millions)")
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300, transparent=True)
    plt.show()


def run():
    plot_sortmerna()

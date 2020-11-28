#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SecretColors import Palette
from matplotlib.patches import Patch
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


def plot_utr_length_distribution(conditions, direction, atlas,
                                 filename="plot.png",
                                 threshold=0, use_limit=None
                                 ):
    plt.clf()
    x_lim = 2000
    colors = [p.cyan, p.magenta]

    # Adjustments
    if direction == "down":
        colors = [p.green, p.orange]
    samples = []
    total = []
    for i, condition in enumerate(conditions):
        genes = get_genes(direction=direction,
                          condition=condition,
                          atlas=atlas,
                          threshold=threshold,
                          use_limit=use_limit)
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
    round_fig = int(np.log10(tt[1])) * -1 + 1
    plt.annotate(f"p-value : {round(tt[1], round_fig)}"
                 f"\n\nLengths above {x_lim}\nare not shown"
                 f"\n\n n = {total[0]} ({conditions[0]}),"
                 f"\n {total[1]} ({conditions[1]})",
                 (0.96, 0.81), ha="right", xycoords="axes fraction", va="top",
                 fontstyle="italic", color=p.gray(shade=70))

    if use_limit is not None:
        tsh = f"{use_limit[0]} <= Log2FC <= {use_limit[1]}"
    else:
        if direction == "up":
            tsh = f"Log2FC >= {threshold}"
        else:
            tsh = f"Log2FC <= {threshold}"
    plt.title(f"{direction} regulated genes "
              f"[{atlas}, {tsh}]")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    # plt.show()


def common_gene_analysis(conditions, *,
                         direction, atlas, threshold, use_limit,
                         filename="plot.png"):
    plt.clf()
    x_lim = 2000
    wt = get_genes(direction=direction,
                   condition=conditions[0],
                   atlas=atlas,
                   use_limit=use_limit,
                   threshold=threshold)
    mia40 = get_genes(direction=direction,
                      condition=conditions[1],
                      atlas=atlas,
                      use_limit=use_limit,
                      threshold=threshold)

    common = list(set(wt).intersection(set(mia40)))
    wt_only = list(set(wt).difference(set(mia40)))
    mia40_only = list(set(mia40).difference(set(wt)))

    new_conds = [common, wt_only, mia40_only]
    labels = [f"common (n={len(common)})",
              f"{conditions[0]} (n={len(wt_only)})",
              f"{conditions[1]} (n={len(mia40_only)})"]
    colors = [p.gray, p.blue, p.red]
    samples = []

    for i in range(len(new_conds)):
        raw = extract_utr_sequence(3)
        utr = filter_utrs(new_conds[i], raw)
        values = [len(x[1]) for x in utr]
        samples.append([x for x in values])
        try:
            density_utr = gaussian_kde(values)
        except ValueError:
            print("Skiped.....")
            print(filename)
            return
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

    round_fig1 = int(np.log10(tt1[1])) * -1 + 1
    round_fig2 = int(np.log10(tt2[1])) * -1 + 1
    plt.annotate(f"p-value : "
                 f"\n common vs {conditions[0]} : {round(tt1[1], round_fig1)}"
                 f"\n common vs {conditions[1]} : {round(tt2[1], round_fig2)}"
                 f"\n\nLengths above {x_lim}\nare not shown",
                 (0.96, 0.78), ha="right", xycoords="axes fraction", va="top",
                 fontstyle="italic", color=p.gray(shade=70))

    if use_limit is not None:
        tsh = f"{use_limit[0]} <= Log2FC <= {use_limit[1]}"
    else:
        if direction == "up":
            tsh = f"Log2FC >= {threshold}"
        else:
            tsh = f"Log2FC <= {threshold}"
    plt.title(f"{direction} regulated genes "
              f"[{atlas}, {tsh}]")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    # plt.show()


def violin_plots():
    direction = "up"
    atlas = "general"

    conditions = ["mars wt", "mars mia40", "hsf wt ", "hsf mia40"]
    colors = [p.blue(shade=40), p.red(shade=40),
              p.blue(shade=40), p.red(shade=40)]
    labels = ["mars", "hsf"]
    samples = []
    counter = 0
    for i, condition in enumerate(conditions):
        genes = get_genes(direction=direction,
                          condition=condition,
                          atlas=atlas)
        raw = extract_utr_sequence(3)
        utr = filter_utrs(genes, raw)
        values = [len(x[1]) for x in utr]
        samples.append([x for x in values])
        vl = plt.violinplot(values,
                            vert=False,
                            positions=[counter],
                            showextrema=False)
        for v in vl:
            if v == "bodies":
                for b in vl[v]:
                    b.set_color(colors[i])
                    b.set_alpha(1)
                    m = np.mean(b.get_paths()[0].vertices[:, 1])
                    if i % 2 == 1:
                        b.get_paths()[0].vertices[:, 1] = np.clip(
                            b.get_paths()[0].vertices[:, 1], -np.inf, m)
                    else:
                        b.get_paths()[0].vertices[:, 1] = np.clip(
                            b.get_paths()[0].vertices[:, 1], m, np.inf)

        if i % 2 != 0:
            counter += 1

    # Welch's unequal variances t-test
    tt1 = ttest_ind(a=samples[0], b=samples[1], equal_var=False)
    tt2 = ttest_ind(a=samples[2], b=samples[3], equal_var=False)

    plt.xlabel("UTR Length")
    plt.gca().set_facecolor(p.gray(shade=10))
    plt.title(f"{direction}-regulated genes ({atlas})")
    plt.yticks(range(len(labels)), labels)
    handles = [
        Patch(facecolor=p.blue(shade=40), label="wt"),
        Patch(facecolor=p.red(shade=40), label="mia40"),
    ]
    plt.annotate(f"p-value : "
                 f"\n MARS (wt vs mia40) : {tt1[1]}"
                 f"\n  HSF (wt vs mia40) : {tt2[1]}",
                 (0.96, 0.5), ha="right",
                 xycoords="axes fraction",
                 va="center",
                 fontstyle="italic", color=p.gray(shade=70))

    plt.legend(handles=handles, loc=0)
    plt.ylim(-0.5, 1.5)
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)
    plt.show()


def generate_combinations():
    atlas = "mitocarta"
    thresholds = [
        0.5, -0.5, 1, -1,
        [-1, -0.5], [0.5, 1]
    ]
    cond_pairs = [
        ("hsf wt", "mars wt", "wt"),
        ("hsf mia40", "mars mia40", "mia40"),
        ("hsf wt", "hsf mia40", "hsf"),
        ("mars wt", "mars mia40", "mars"),
    ]

    fig_pair = []
    com_pair = []
    current = []
    current_com = []
    for pair in cond_pairs:
        for threshold in thresholds:
            use_limit = None
            if isinstance(threshold, list):
                t = threshold[0]
                use_limit = threshold
                th_str = "_".join([str(x) for x in threshold])
            else:
                t = threshold
                th_str = f"{threshold}"
            if t < 0:
                direction = "down"
            else:
                direction = "up"

            th_str = th_str.replace("-", "")
            filename = f"figs/{atlas[:3]}_{direction}_{pair[2]}_{th_str}"
            filename = filename.replace("-", "")
            filename_com = f"{filename}_cog.png"
            filename += ".png"
            plot_utr_length_distribution(direction=direction,
                                         conditions=(pair[0], pair[1]),
                                         atlas=atlas,
                                         threshold=threshold,
                                         use_limit=use_limit,
                                         filename=filename)
            common_gene_analysis(direction=direction,
                                 conditions=(pair[0], pair[1]),
                                 atlas=atlas,
                                 threshold=threshold,
                                 use_limit=use_limit,
                                 filename=filename_com)

            current.append(filename)
            current_com.append(filename_com)
            if len(current) == 2:
                fig_pair.append((th_str, current))
                com_pair.append((th_str, current_com))
                current = []
                current_com = []

    with open("figs.tex", "w") as f:
        for d in fig_pair:
            print(get_latex_fig(d, "template"), file=f)
        for d in com_pair:
            print(get_latex_fig(d, "template2", prefix="cog_"), file=f)


def get_latex_fig(data, template, prefix=""):
    log2 = data[0]
    files = sorted(data[1])
    with open(template) as f:
        tp = f.readlines()
    lines = "".join(tp)
    lines = lines.replace("$FIG1", files[0])
    lines = lines.replace("$FIG2", files[1])
    label = str(files[0]).replace("figs/gen_", "")
    label = label.replace("up_", "").replace("down_", "").replace(".png", "")
    label = prefix + label
    lines = lines.replace("$LAB", label)
    lines = lines.replace("$COND", label.split("_")[0])
    if "_" in log2:
        log_txt = f"Log2FC between {' and '.join(log2.split('_'))}"
    else:
        log_txt = f"Log2FC threshold of {log2}"

    lines = lines.replace("$LOG2FC", log_txt)
    return lines


def run():
    generate_combinations()

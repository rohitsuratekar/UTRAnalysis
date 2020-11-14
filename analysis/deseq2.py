#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  All analysis related to DESeq2 output

import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pygraphviz as pgv
from SecretColors import Palette
from SecretColors.utils import text_color
from gprofiler import GProfiler
from matplotlib.patches import Patch

matplotlib.rc("font", family="IBM Plex Sans")

p = Palette()

SOURCE_NODES = ["GO:0005575", "KEGG:00000", "GO:0008150", "GO:0003674"]


def go_enrichment_analysis(filename: str):
    df = pd.read_csv(filename)
    df = df[df["padj"] < 0.05]
    df = df.sort_values(by="log2FoldChange", ascending=False)
    up_genes = df[df["log2FoldChange"] > 0]["gene_id"].to_numpy()
    # Need to switch the order so that more negative will come on top as
    # those are more significant
    df = df.sort_values(by="log2FoldChange", ascending=True)
    down_genes = df[df["log2FoldChange"] < 0]["gene_id"].to_numpy()
    gp = GProfiler(return_dataframe=True)
    prof = gp.profile(organism="drerio",
                      ordered=True,
                      no_evidences=True,
                      query={
                          "up": list(up_genes),
                          "down": list(down_genes)
                      },
                      domain_scope="known",
                      user_threshold=0.05)
    prof.to_csv("out.csv", index=False)


# 'source', 'native', 'name', 'p_value', 'significant', 'description',
#        'term_size', 'query_size', 'intersection_size', 'effective_domain_size',
#        'precision', 'recall', 'query', 'parents'

def draw_network(source):
    df = pd.read_csv("out.csv", encoding="utf8")
    df = df[df["p_value"] < 0.05]
    df = df.sort_values(by="p_value")
    df = df[df["source"] == source]
    df["p_value"] = -1 * np.log(df["p_value"])
    # df["p_value"] = np.log(df["p_value"])
    max_p = df["p_value"].max()
    # df["p_value"] = 20 + df["p_value"] * 50 / max_p
    df["p_value"] = df["p_value"] / max_p

    colors = dict(zip(df['native'], df["p_value"]))

    name_maps = {}
    # Add all children in given order
    counter = 1

    for c in df["native"].values:
        if c not in name_maps:
            name_maps[c] = counter
            counter += 1

    graph = pgv.AGraph(directed=True)
    for (i, d) in df.iterrows():
        parents = d["parents"]
        parents = [re.sub("[^A-Za-z0-9:]+", "", x) for x in parents.split(",")]
        parents = [x for x in parents if len(x) > 0]
        child = d["native"]

        for k in parents:
            if k not in name_maps.keys():
                name_maps[k] = counter
                counter += 1
            graph.add_edge(name_maps[k], name_maps[child])

    id_mapping = {v: k for k, v in name_maps.items()}
    graph.graph_attr['dpi'] = 300
    graph.graph_attr['bgcolor'] = "transparent"
    graph.graph_attr['overlap'] = False
    for n in graph.nodes():
        n.attr['style'] = "filled"
        n.attr['color'] = p.gray(shade=20)
        name = id_mapping[int(n.name)]
        if name in colors:
            clr = 15 + colors[name] * 70
            clr = p.red(shade=clr)
            n.attr['color'] = clr
            n.attr['fontcolor'] = text_color(clr)
            n.attr['fontsize'] = 18

    graph.layout()
    graph.draw("plot.png")


def map_terms(no_of_terms, source, condition, axis_offset=50):
    # plt.figure(figsize=(9, 5))
    plt.figure(figsize=(7, 5))
    up_color = p.blue
    down_color = p.red
    df = pd.read_csv("out.csv")
    # Remove the base node
    df = df[~df["native"].isin(SOURCE_NODES)]
    df = df.sort_values(by="p_value")
    df["p_value"] = np.log(df["p_value"])
    df = df[df["source"] == source]
    upr = df[df["query"] == "up"]
    up_mapping = tuple(zip(upr["name"], upr["p_value"]))
    up_mapping = sorted(up_mapping, key=lambda x: x[1])
    down = df[df["query"] == "down"]
    down_mapping = tuple(zip(down["name"], down["p_value"]))
    down_mapping = sorted(down_mapping, key=lambda x: x[1])

    max_no = max(abs(up_mapping[0][1]), abs(down_mapping[0][1])) + axis_offset

    for i in range(no_of_terms):
        up_no = -1 * up_mapping[i][1]
        down_no = down_mapping[i][1]
        plt.barh(i, max_no, color=up_color(shade=15))
        plt.annotate(up_mapping[i][0], xy=(5, i),
                     va="center", ha="left")
        plt.barh(i, up_no, color=up_color())
        plt.barh(i, -1 * max_no, color=down_color(shade=15))
        plt.barh(i, down_no, color=down_color(shade=40))
        plt.annotate(down_mapping[i][0], xy=(-5, i),
                     va="center", ha="right")

    plt.xlim([-1 * max_no, max_no])
    plt.yticks([])
    plt.xlabel("Log$_{10}$ (p-value) [for Up: -1xValue]")

    handles = [
        Patch(facecolor=down_color(shade=40), label="Down"),
        Patch(facecolor=up_color(), label="Up"),

    ]
    plt.legend(handles=handles, loc="upper right",
               ncol=2,
               frameon=False,
               bbox_to_anchor=(1.01, 1.08))
    plt.annotate(f"{source} terms enrichment in {condition}", xy=(0, 1.02),
                 va="bottom",
                 xycoords="axes fraction")
    plt.gca().invert_yaxis()
    plt.savefig("plot.png", dpi=300, transparent=True)
    plt.show()


def plot_network(source):
    df = pd.read_csv("out.csv")
    df = df[df["p_value"] < 0.05]
    df = df[df["source"] == source]
    df = df.sort_values(by="p_value")
    df = df[~df["native"].isin(SOURCE_NODES)]
    up = df[df["query"] == "up"]
    down = df[df["query"] == "down"]
    counter = 1
    names = {}
    for it in df["native"].values:
        if it not in names:
            names[it] = counter
            counter += 1

    graph = pgv.AGraph(directed=True)
    for _, d in df.iterrows():
        parents = d["parents"]
        parents = [re.sub("[^A-Za-z0-9:]+", "", x) for x in parents.split(",")]
        for k in parents:
            if k not in names:
                names[k] = counter
                counter += 1
            graph.add_edge(names[k], names[d["native"]])

    for n in names:
        cond_up = n in up['native'].values
        cond_down = n in down['native'].values
        nd = graph.get_node(names[n])
        if any([cond_up, cond_down]):
            nd.attr['style'] = "filled"
        if all([cond_up, cond_down]):
            nd.attr['color'] = p.gray()
        elif cond_up:
            nd.attr['color'] = p.blue(shade=40)
        elif cond_down:
            nd.attr['color'] = p.red(shade=40)

    graph.graph_attr['overlap'] = False
    graph.graph_attr['K'] = 10
    graph.layout("sfdp")
    graph.draw("plot.png")


def run():
    source = "KEGG"
    # go_enrichment_analysis("deseq2/salmon_hsf_vs_whole.csv")
    map_terms(10, source, "HSF/W (WT)", axis_offset=10)

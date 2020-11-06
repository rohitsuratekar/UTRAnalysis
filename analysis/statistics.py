#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
from matplotlib.collections import LineCollection
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from SecretColors import Palette
from matplotlib.patches import Patch
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster import hierarchy

matplotlib.rc('font', family='IBM Plex Sans')

p = Palette()
BASE_FOLDER = "salmon_counts"
CONDITIONS = ["wt_hsf", "wt_mars", "wt_whole"]
COND_COLORS = [p.cyan(shade=50),
               p.magenta(shade=50),
               p.yellow(shade=30),
               p.gray(shade=50),
               p.blue(shade=30),
               p.red(shade=30)]
REPLICATES = 3


def get_file_names():
    all_samples = []
    c = 0
    for name in CONDITIONS:
        for i in range(REPLICATES):
            all_samples.append(f"{BASE_FOLDER}/{name}_{i + 1}_quant.sf")
        c += 1
    return all_samples


def extract_expression(filenames: list):
    col = "TPM"
    name = "Name"
    all_dfs = []
    for f in filenames:
        df = pd.read_csv(f, sep="\t")
        df = df.set_index(name)
        df = df[[col]]
        df = df.rename(columns={col: f})
        all_dfs.append(df)
    df = pd.concat(all_dfs, axis=1, join="inner")
    return df.to_numpy()


def perform_pca(tpms):
    pca = PCA()
    pca.fit(tpms)
    var = pca.components_
    lc = []
    colors = []
    for con, col in zip(CONDITIONS, COND_COLORS):
        lc.append(Patch(label=con, color=col))
        for i in range(REPLICATES):
            colors.append(col)

    percentages = [round(x, 2) for x in pca.explained_variance_ratio_ * 100]
    plt.legend(handles=lc, loc=0)
    plt.scatter(var[0, :], var[1, :], color=colors, s=100)
    plt.xlabel(f"PC1: {percentages[0]} %")
    plt.ylabel(f"PC2: {percentages[1]} %")
    plt.gca().set_facecolor(p.gray(shade=10))
    plt.grid(axis="both", ls=":", alpha=0.7)
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)
    plt.show()


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    labels = []
    for c in CONDITIONS:
        for i in range(REPLICATES):
            labels.append(f"{c}_{i + 1}")
    labels = np.asarray(labels)
    # Plot the corresponding dendrogram
    hierarchy.set_link_color_palette(COND_COLORS)
    dendrogram(linkage_matrix, **kwargs, labels=labels)


def hierarchical_clustering(tpms):
    clustering = AgglomerativeClustering(
        distance_threshold=0, n_clusters=None,
        affinity="manhattan",
        linkage="complete")
    model = clustering.fit(tpms.transpose())
    # fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)  # type: plt.Axes
    plot_dendrogram(model, truncate_mode=None, ax=ax,
                    orientation="right")
    plt.xlabel("Distance")
    for c in ax.get_children():
        if isinstance(c, LineCollection):
            c.set_linewidth(2)
    ax.set_facecolor(p.gray(shade=10))
    plt.tight_layout()
    plt.savefig("plot.png", dpi=300)
    plt.show()


def plot_ma():
    pass


def run():
    names = get_file_names()
    tpms = extract_expression(names)
    hierarchical_clustering(tpms)
    # perform_pca(tpms)

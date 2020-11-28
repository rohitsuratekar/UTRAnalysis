#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics

import pygraphviz as pgv


def draw_method():
    graph = pgv.AGraph(directed=True)

    node_mapping = {
        "a": "Gene List",
        "b": "DE gene below FDR 0.05",
        "c": "3'UTR sequences from BioMart",
        "d": "Select longest transcript for any gene"
    }

    graph.add_edge("a", "b")
    graph.add_edge("c", "d")

    for node in graph.nodes():
        node.attr['label'] = node_mapping[node]

    graph.layout("dot")
    graph.draw("plot.png")


def run():
    pass

#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#  Homer de-novo motif predictions

import os
import subprocess


def run_homer():
    folder = "homer"
    if not os.path.isdir(folder):
        os.mkdir(folder)

    target_suffix = [
        "_mitocarta_mars_mia40_3utr.fasta",
        "_mitocarta_mars_wt_3utr.fasta",
        "_mitocarta_hsf_mia40_3utr.fasta",
        "_mitocarta_hsf_wt_3utr.fasta"
    ]
    control = "./fasta/control.fasta"
    targets = []
    prefix = ["up", "down"]
    for t in target_suffix:
        for p in prefix:
            targets.append(f"{p}{t}")

    homer_path = ":/home/dex/Softwares/homer/bin"
    env = os.environ.copy()
    env['PATH'] += homer_path

    for target in targets:
        output = target.replace("_3utr.fasta", "")
        output = f"{folder}/{output}"
        target = f"./fasta/{target}"
        random_opts = [
            "findMotifs.pl",
            target,
            "fasta",
            f"{output}_control",
            '-rna',
            '-homer2',
            "-p",
            "7",
        ]
        all_opts = [x for x in random_opts]
        all_opts[3] = f"{output}_all"
        all_opts.extend(["-fasta", control])
        subprocess.run(random_opts, env=env)
        subprocess.run(all_opts, env=env)


def run():
    run_homer()

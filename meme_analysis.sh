#!/usr/bin/env sh
# Copyright (c) 2020.
# Author: Rohit Suratekar, IIMCB
#
# UTR Analysis and related statistics

MEME_PATH="/home/dex/Softwares/meme/bin/meme"
FASTA="5utr.fasta"

export PATH=$HOME/Softwares/meme/bin:$HOME/Softwares/meme/libexec/meme-5.1.1:$PATH

$MEME_PATH -dna $FASTA

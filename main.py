#  Copyright (c) 2020.
#  Author: Rohit Suratekar, IIMCB
#
#  UTR Analysis and related statistics
#
#  utr5.fasta, utr3.fasta, trans_to_gene.csv files (in data folder)
#  are downloaded from the  ENSEMBL's BioMart portal on 20 Oct 2020.
#  https://www.ensembl.org/biomart
#  Ensembl Genes 101 > Zebrafish Genes (GRCz11) > Export


from analysis.deseq2 import run

run()

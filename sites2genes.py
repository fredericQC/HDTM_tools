#!/usr/bin/env python
# coding: utf-8

"""
Usage : python insertions_to_genes.py CDS.bed intersect.bed
Input : CDS.bed intersect.bed
Output: genes_insertions.tsv
"""

import sys
import numpy as np
import pandas as pd

length_pc_start = 5.0
length_pc_stop = 85.0

# Read files
gene_df = pd.read_csv(
    sys.argv[1],
    dtype={
        "chr_name": str,
        "start": "Int64",
        "end": "Int64",
        "name": str,
        "score": "Int64",
        "strand": str
    },
    sep="\t",
    header=None,
    names=["chr_name", "start", "end", "name", "score", "strand"],
)

intersect_df = pd.read_csv(
    sys.argv[2],
    dtype={
        "chr_name_gene": str,
        "start_gene": "Int64",
        "end_gene": "Int64",
        "name_gene": str,
        "score_gene": "Int64",
        "strand_gene": str,
		"chr_name_insertion": str,
		"start_insertion": "Int64",
		"end_insertion": "Int64",
		"score_insertion": "Int64"
    },
    sep="\t",
    header=None,
	usecols=range(10),
    names=["chr_name_gene", "start_gene", "end_gene", "name_gene", "score_gene", "strand_gene", "chr_name_insertion", "start_insertion", "end_insertion", "score_insertion"],
)

# filter insertions outside of genes (due to insertion duplication and excluding length_pc_start length_pc_end)
intersect_df.replace({'strand_gene': {"+": "1", "-": "-1"}},inplace=True)
gene_df.replace({'strand': {"+": "1", "-": "-1"}},inplace=True)

intersect_df = intersect_df.loc[
	(
	(intersect_df['strand_gene']=="1") &
	((intersect_df['start_insertion'] - 5) >= (intersect_df['start_gene']+((intersect_df['end_gene']-intersect_df['start_gene']) * length_pc_start/100.0))) &
	((intersect_df['start_insertion'] + 5) <= (intersect_df['end_gene'] - 1 - ((intersect_df['end_gene']-intersect_df['start_gene']) * (100-length_pc_stop)/100.0)))
	)
	|
	(
	(intersect_df['strand_gene']=="-1") &
	((intersect_df['start_insertion']+5) <= (intersect_df['end_gene'] - 1 - ((intersect_df['end_gene']-intersect_df['start_gene'])*length_pc_start/100.0))) &
	((intersect_df['start_insertion']-5) >= (intersect_df['start_gene'] + ((intersect_df['end_gene']-intersect_df['start_gene'])*(100-length_pc_stop)/100.0)))
	)
]

# make columns required for biotradis toolkit (almost, usually no chr_name)
gene_df["locus_tag"] = gene_df["name"]
gene_df["gene_name"] = gene_df["name"]
gene_df["ncrna"] = 0
gene_df["read_count"] = 0
gene_df["ins_index"] = 0.0
gene_df["gene_length"] = gene_df["end"] - gene_df["start"]
gene_df["ins_count"] = 0
gene_df["fcn"] = "NA"
gene_df = gene_df.set_index('name')

# Compute insertion index and read_count
intersect_df = intersect_df.groupby(['name_gene'],sort=False).agg(['sum','mean','count'])
for index, gene in intersect_df.iterrows():
	gene_df.loc[gene.name,"read_count"] = gene[('score_insertion', 'sum')]
	gene_df.loc[gene.name,"ins_count"] = gene[('score_insertion', 'count')]
gene_df["ins_index"] = gene_df["ins_count"]/gene_df["gene_length"]

gene_df[["chr_name","locus_tag", "gene_name", "ncrna", "start", "end", "strand", "read_count", "ins_index", "gene_length", "ins_count", "fcn"]].to_csv(
    sys.argv[3], header=True, index=False, sep="\t"
)

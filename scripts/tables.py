#!/usr/bin/env python
# Encoding: utf8

'''

'''
from Bio import SeqIO
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
import os

def calculate_statistics(genomes_lengths, genes_file):
    #print(genes_file)
    # read genes_table for the genome
    # print(genes_file)
    # print(str(genes_file))
    # print(os.path.abspath(genes_file))
    table = pd.read_csv(genes_file, header=0, sep="\t")

    # separate according to types and count
    cds = table[table["type"] == "CDS"]
    trna = table[table["type"] == "tRNA"]
    n_cds = cds.shape[0]
    n_trna = trna.shape[0]

    # get genome_id and caller
    genome_id = table["genome"][0]
    caller    = table["caller"].value_counts().idxmax()

    # compute coding length, only cds. It colapses overlaping ORFs in the same
    # interval
    start = cds["start"].tolist()
    end   = cds["end"].tolist()

    intervals = [[s,e] for s,e in zip(start, end)]

    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)

    # coding length
    coding_bases = float(0)
    for interval in merged:
        coding_bases += 1 + (interval[1] - interval[0])

    tow = [
            genome_id,
            genomes_lengths[genome_id],
            str(n_cds),
            "{:.3f}".format(cds.mean(axis=0)["aa_length"]),
            "{:.3f}".format(coding_bases/float(genomes_lengths[genome_id])),
            caller,
            str(n_trna)
          ]

    return tow




def coding_statistics(tables_genes, lengths_file, output_file):
    lines = [line.strip().split("\t") for line in open(lengths_file).readlines()]
    lengths = {line[0]:line[1] for line in lines}

    # func = partial(calculate_statistics, lengths)
    # pool = multiprocessing.Pool(processes=3)
    # lines_calculated = pool.map(func, tables_genes)
    # pool.close()
    # pool.join()

    # with open(output_file, "w") as fout:
    #     for line in lines_calculated:
    #         fout.write("{}\n".format("\t".join(line)))


    with open(output_file, "w") as fout:
        for table in tables_genes:
            print(table)
            tow = calculate_statistics(lengths, table)
            fout.write("{}\n".format("\t".join(tow)))

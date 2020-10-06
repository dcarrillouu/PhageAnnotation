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
    print(genes_file)
    # read genes_table for the genome
    table = pd.read_csv(genes_file, header=0, sep="\t")

    # separate according to types and count
    cds = table[table["type"] == "CDS"]
    trna = table[table["type"] == "tRNA"]
    n_cds = cds.shape[0]
    n_trna = trna.shape[0]

    # get genome_id and caller
    genome_id = cds["genome"][0]
    caller    = cds["caller"][0]

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

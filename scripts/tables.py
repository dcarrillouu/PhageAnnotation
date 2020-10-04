#!/usr/bin/env python
# Encoding: utf8

'''

'''

import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
import os

def genes2table(fasta_cds, trna_file, outable, genome, caller):
    '''

    '''
    # Initialize output list
    final_table = list()

    # Process fasta cds
    records = SeqIO.parse(fasta_cds, "fasta")


    for record in records:
        splited = record.id.split("|")

        start = splited[2].split("..")[0]
        end   = splited[2].split("..")[1]

        towt  = [
                 genome,
                 start,
                 end,
                 splited[4],
                 splited[3],
                 caller,
                 "CDS",
                 record.id,
                 splited[5]
                 ]

        final_table.append(towt)


    # Process trna
    # don't strip() so the 'notes' column is preserved when empty
    lines = [line.split("\t") for line in open(trna_file).readlines()[3:]]
    for line in lines:
        print("#########################", line)
        # Bitscore filter >35
        if float(line[-2]) >= 35:

            start = str(min(int(line[2]), int(line[3])))
            end   = str(max(int(line[2]), int(line[3])))

            if int(line[2]) < int(line[3]):
                strand = "+"
            else:
                strand = "-"


            towt = [
                    line[0],
                    start,
                    end,
                    strand,
                    "",
                    "tRNAscan-SE",
                    "tRNA",
                    "tRNA_{}".format(line[1]),
                    "tRNA_{}".format(line[1])
                    ]
            final_table.append(towt)

    # Format final table
    # Header
    header = ["genome", "start", "end", "strand", "aa_length", "caller", "type", "gene_header", "gene_id"]

    # sort by start position
    final_table_sorted = sorted(final_table, key=lambda line: int(line[1]))

    # write file
    with open(outable, "w") as fout:
        fout.write("\t".join(header) + '\n')
        for line in final_table_sorted:
            fout.write("{}\n".format("\t".join(line)))

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





# Notice lengths_file may change in the future. By now I compute the genome
# length with BioPython and save to file.
def coding_statistics(tables_genes, lengths_file, output_file):
    # '''
    # coding density (%)
    # gene length (mean)
    # number of ORFs
    # number of tRNAs
    # ''''
    # read lengths file and save to dict
    print(tables_genes)
    lines = [line.strip().split("\t") for line in open(lengths_file).readlines()]
    lengths = {line[0]:line[1] for line in lines}

    func = partial(calculate_statistics, lengths)
    pool = multiprocessing.Pool(processes=3)
    lines_calculated = pool.map(func, tables_genes)
    pool.close()
    pool.join()

    with open(output_file, "w") as fout:
        for line in lines_calculated:
            fout.write("{}\n".format("\t".join(line)))

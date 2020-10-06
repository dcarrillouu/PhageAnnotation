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



def genes2table(fasta_cds, trna_file, outable, genome, caller):
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
        # Bitscore filter >35
        if float(line[-2]) >= 35:

            start = str(min(int(line[2]), int(line[3])))
            end   = str(max(int(line[2]), int(line[3])))

            if int(line[2]) < int(line[3]):
                strand = "+"
            else:
                strand = "-"


            towt = [
                    line[0].strip(), # 'genome id' column of tRNAscan has an extra space
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

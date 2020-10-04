#!/usr/bin/env python
# Encoding: utf8

'''

'''
from Bio import SeqIO
import os



def fix_headers(infile, outfile, mode):
    print("####################################")
    print(infile)
    ("####################################")

    records = SeqIO.parse(infile, "fasta")

    with open(outfile, "w") as fout:

        for record in records:
            genome = "_".join(record.id.split("_")[:-1])

            descrp = record.description.split(" # ")
            start  = descrp[1]
            end    = descrp[2]

            aa_length = int((int(end) - int(start) + 1) / 3)

            if descrp[3] == "1":
                strand = "+"
            else:
                strand = "-"

            gene_id = descrp[-1].split("ID=")[1].split(";")[0]

            header_out = f'{genome}|prodigal-{mode}|{start}..{end}|{aa_length}|{strand}|{gene_id}'

            record.id = header_out
            record.description = ""

            SeqIO.write(record, fout, "fasta")

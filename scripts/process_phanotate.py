#!/usr/bin/env python
# Encoding: utf8

'''

'''

import argparse
from pathlib import Path
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def extract_ORF_sequences(genome_fasta, table):
    # get the id of the file/genome
    #id = os.path.basename(table).replace("_phnot.tbl","")

    # read fasta genome file
    genome_record = SeqIO.read(genome_fasta, "fasta")

    # read PHANOTATE table
    tbl = [line.strip().split("\t") for line in open(table) if not line.startswith("#")]

    # open files to write the sequences
    out_nt  = open(table.replace(".tab", ".fna"), "w")
    out_aa  = open(table.replace(".tab", ".faa"), "w")
    #out_tbl = open(os.path.join(outdir, id + "_phnot.tblgenes"), "w")

    # go through the table and write the sequences
    # subtract one to the start position
    cont = 0
    for line in tbl:
        start = min(int(line[0]), int(line[1])) - 1
        end   = max(int(line[0]), int(line[1]))
        strand = line[2]


        if strand == "+":
            seq = genome_record.seq[start:end]
        else:
            seq = genome_record.seq[start:end].reverse_complement()

        aa_length = int(len(seq)/3)
        header_out = f'{genome_record.id}|phanotate|{start+1}..{end}|{aa_length}|{strand}|gene_{cont}'

        nt_record = SeqRecord(seq, id=header_out, description='')
        aa_record = SeqRecord(seq.translate(), id=header_out, description='')

        SeqIO.write(nt_record, out_nt, "fasta")
        SeqIO.write(aa_record, out_aa, "fasta")

        # write to the tablegenes file
        #out_tbl.write(f'{line[3]}\t{start+1}\t{end}\t{aa_length}\t{strand}\tgene_{cont}\t{header_out}\tphanotate\n')

        cont += 1

    out_nt.close()
    out_aa.close()
    #out_tbl.close()

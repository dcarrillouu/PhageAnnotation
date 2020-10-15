'''
vContact2 needs as input the proteins and a table describing which genome they
belong to, like this:

protein_id,contig_id,keywords
ref|NP_039777.1|,Sulfolobus spindle-shaped virus 1,ORF B-251
ref|NP_039778.1|,Sulfolobus spindle-shaped virus 1,ORF D-335
ref|NP_039779.1|,Sulfolobus spindle-shaped virus 1,ORF E-54

This script takes as input a fasta with all the proteins that were selected as
best-coding (prodigal-11, TAG or TGA) and generates such table.

Notice I can do this because the genome_id information is in the header of the
protein, which is like this:

>crAssphage|prodigal-11|191..346|52|+|1_1

input fasta file was generated as follows:

'''

import argparse
from pathlib import Path
import os
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_proteins_file',
                               dest='in_faa',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-o', '--output_table',
                               dest='out_table',
                               type=lambda p: Path(p).resolve(),
                               required=True,
                               help='Extension must be csv or vContact2 will fail'
                               )

    return parser.parse_args()

def main():
    args = parse_args()

    with open(args.out_table, "w") as fout:
        fout.write("protein_id,contig_id,keywords\n")

        records = SeqIO.parse(args.in_faa, "fasta")
        for record in records:
            tow = list()
            genome = record.description.split("|")[0]
            tow = f'{record.description},{genome},not_provided\n'
            fout.write(tow)







if __name__ == '__main__':
    main()

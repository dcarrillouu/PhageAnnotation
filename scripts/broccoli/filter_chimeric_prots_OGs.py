'''
Takes the OG raw fasta files and removes chimeric proteins. As input, it needs
the "chimeric_proteins.txt" file from Broccoli dir_step3 and the folder with all
the raw OG_XX.faa
'''

import argparse
import glob
import os
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-f', '--OG_faa_dir',
                               dest='OG_faa_dir',
                               required=True,
                               help='directory with genome tables annotated with OG and yutin.'
                               'Generated with "yutin_to_genome-tables.py", extension ".genes-OG-profiles"'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help='dir to put the filtered OG_XX_filtered.faa'
                               )

    requiredArgs.add_argument('-c', '--chimeric_file',
                               dest='chimeric_file',
                               required=True,
                               help='file "chimeric_proteins.txt"'
                               )


    return parser.parse_args()


def main():

    args = parse_args()

    # read chimeric file information and store id for the chimeric proteins
    chimeric = [line.split("\t")[1] for line in open(args.chimeric_file).readlines()[1:]]
    print(chimeric)

    OG_faas  = glob.glob(f"{args.OG_faa_dir}/*faa")
    for faa in OG_faas:
        outfile = f"{args.out_dir}/" + os.path.basename(faa).split(".faa")[0] + "_chim.faa"

        tow = list()
        records = SeqIO.parse(faa, "fasta")
        for record in records:
            if record.id not in chimeric:
                tow.append(record)
            else:
                print(record.id)

        SeqIO.write(tow, outfile, "fasta")



















if __name__ == "__main__":
    main()

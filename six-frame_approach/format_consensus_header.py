'''
Script to format the header name of the consensus sequence after running hhconsensus
'''

import argparse
import glob
import os
from Bio import SeqIO, SearchIO



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_dir',
                               dest='input_dir',
                               required=True,
                               help='directory with the mafft aligned OGs containing '
                               'a consensus sequence.'
                               )


    return parser.parse_args()


def main():

    args = parse_args()

    # get the mafft msa files
    msa_files = glob.glob(f"{args.input_dir}/*.mafft-einsi-cons")

    for file in msa_files:
        # get the OG_id from the file name
        og_id = os.path.basename(file).split(".")[0]
        print(og_id)
        # read file as lines. Discard first "#" line introduced by hhconsensus
        lines = [line for line in open(file).readlines() if not line.startswith("#")]

        tow = list()
        for line in lines:
            if ">" in line and "_consensus" in line:
                tow.append(f">{og_id}_consensus\n")
            else:
                tow.append(line)

        print(tow[0])
















if __name__ == "__main__":
    main()

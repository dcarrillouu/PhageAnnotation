'''
Collapse all the genomes from Bas at the family level.
'''


import argparse
import os
import glob
from functools import partial
import multiprocessing


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_dir',
                               dest='in_dir',
                               required=True,
                               help='Directory with the genomes fasta files'
                               )
    requiredArgs.add_argument('-out', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help='Directory to store the results'
                               )

    requiredArgs.add_argument('-t', '--taxa_file',
                               dest='taxa_file',
                               required=True,
                               help='Tabular file with taxonomic assignments for the genomes'
                               )



    return parser.parse_args()



def main():

    args = parse_args()

    # read taxa file and store in 2 dicts: k=taxa , v=[genomes] ; k=genome, v=[fam, subfam, genus]
    subfamilies = dict()

    lines = [line.strip().split("\t") for line in open(args.taxa_file, "r").readlines()[1:]] # discard header line
    # notice the column slices are very specific to this file. It will vary with other taxa file, like Yutin
    taxa = {line[1]:line[17:20] for line in lines}

    for line in lines:
        if line[18] not in subfamilies:
            subfamilies[line[18]] = [line[1]]
        else:
            subfamilies[line[18]].append(line[1])



    for subfamily, genomes in subfamilies.items():
        to_cat = [f"{args.in_dir}/{genome}.fasta" for genome in genomes]
        os.system("cat {} > {}/{}.fasta".format(" ".join(to_cat), args.out_dir, subfamily))


if __name__ == "__main__":
    main()

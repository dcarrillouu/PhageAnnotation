'''

'''

import argparse
import multiprocessing
from functools import partial
import pandas as pd
import numpy as np
import glob
import os


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--tables_directory',
                               dest='tables_dir',
                               required=True,
                               help='directory with genome tables annotated with OG.'
                               'Generated with "broccoli_to_gggenes-tables.py", extension ".genes-OG"'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-og', '--OG_ID',
                               dest='og_id',
                               required=True,
                               help='plot the context por this OG_ID'
                               )

    requiredArgs.add_argument('-n', '--n_context',
                               dest='n_context',
                               required=False,
                               type=int,
                               default=2,
                               help='number of context genes to plot'
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    # get genomes tables with genes and OG information
    tables = glob.glob(f"{args.tables_dir}/*.genes-OG")

    tow = list()
    n_context = int(args.n_context)

    # iterate trhough the tables
    for table in tables:
        # iterate through the lines of the tables, and store the index in case of a match with the OG
        lines = [line.split("\t") for line in open(table).readlines()[1:]]

        for i in range(0, len(lines)):
            if lines[i][-1].replace("\n", "") == args.og_id:

                # check if there is enough context to plot upstream...
                if (i - n_context) < 0:
                    for j in range(0, i+n_context+1):

                        tow.append(lines[j])
                # ... and downstream
                elif (i + n_context) > len(lines):
                    for j in range(i-n_context, len(lines)):
                        tow.append(lines[j])
                # if there is enough context upstream and downstream
                else:
                    for j in range(i-n_context, i+n_context):

                        tow.append(lines[j])

    header = ["genome", "start", "end", "strand", "complete_id", "OG"]
    tow.insert(0, header)

    with open(f"{args.out_dir}/{args.og_id}_n{str(n_context)}.txt", "w") as fout:
        for line in tow:
            fout.write("\t".join(line))

if __name__ == "__main__":
    main()

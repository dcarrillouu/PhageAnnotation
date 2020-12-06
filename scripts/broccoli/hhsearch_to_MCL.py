'''
 This script takes the hhsearch blasttab files result of comparing the OG profiles
 between them, and returns a unique file with three columns: query OG, hit OG and evalue.
 Notice that the name of the query OG in the files is the name of the first sequence
 in its alignment, so I have to pick the filename up, that contains the name of the OG
'''

import glob
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--OG_blasttab_dir',
                               dest='tab_dir',
                               required=True,
                               help='directory with the tables resulting from running hhsearch, '
                               'one profile vs the others'
                               )
    requiredArgs.add_argument('-o', '--output_file',
                               dest='out_file',
                               required=True,
                               help='file to store the 3 columns outfile usable by MCL'
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    tow = list()

    tables = glob.glob(f"{args.tab_dir}/*.tab")
    for table in tables:
        # remove all the tags from the file to retain only the OG_id
        query_OG = os.path.basename(table).split(".tab")[0]

        # read file
        # discard first line since it is the hit to itself
        lines = [line.strip().split("\t") for line in open(table).readlines()[1:]]
        for line in lines:
            if float(line[-2]) < 0.01:
                tow.append([query_OG, line[1], line[-2]])

    with open(args.out_file, "w") as fout:
        for hit in tow:
            fout.write("\t".join(hit) + "\n")






if __name__ == "__main__":
    main()

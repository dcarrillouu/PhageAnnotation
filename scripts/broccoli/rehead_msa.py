'''
This script changes the identifier (header) of the sequences in a MSA so it does
not give any problem when building and annotating a tree.

It returns a re-head msa and a tabular file with the relationship rehead - original_header.
Edit: add the taxa info
'''

import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--in_msa_file',
                               dest='in_file',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--out_msa_file',
                               dest='out_file',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-t', '--taxa_file',
                               dest='taxa_file',
                               required=True,
                               help=''
                               )



    return parser.parse_args()




def main():

    args = parse_args()

    # read taxa_file and store in a dictionary, key=genome_id, value=taxa_ranks
    lines = [line.strip().split("\t") for line in open(args.taxa_file, "r").readlines()[1:]] # discard header line
    # notice the column slices are very specific to this file. It will vary with other taxa file, like Yutin
    taxa_info = {line[1]:line[17:20] for line in lines}


    msa_out = open(args.out_file, "w")
    rel_out = open(os.path.dirname(args.out_file) + "/rehead.txt", "w")
    rel_out.write("rehead\toriginial_id\ttaxa\n")

    cont = 1
    lines = open(args.in_file).readlines()
    # iterate alignment file lines
    for line in lines:
        if line.startswith(">"):
            genome_id = line.split("|")[0].replace(">", "")
            taxa = "NA"
            if genome_id in taxa_info:
                taxa = taxa_info[genome_id][0]

            rel_out.write("a{}b\t{}\t{}\n".format(cont, line.strip().replace(">", ""), taxa))
            msa_out.write(">a{}b\n".format(cont))
            cont += 1
        else:
            msa_out.write(line)

















if __name__ == "__main__":
    main()

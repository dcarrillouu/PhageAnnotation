'''
It takes the 'orthologous_groups.txt' from Broccoli, step3, and creates a fasta
file per OG with all its proteins. It also requires a big fasta file with all
the proteins that were analyzed.

It requires 'seqtk' to be in $PATH
'''

import argparse
import os
import multiprocessing
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-og', '--og_file',
                               dest='og_file',
                               required=True,
                               help='"orthologous_groups.txt" from Broccoli step 3'
                               )
    requiredArgs.add_argument('-f', '--fasta_file',
                               dest='fasta_file',
                               required=True,
                               help='big fasta file with all the proteins. The script '
                               'will look in this file for the OG proteins, so be sure '
                               'they are there'
                               )

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def extract_sequences_OG(fasta_big_file, out_dir, og):
    og_id = og[0]
    proteins = og[1].split(" ")
    tmp_file = f"{out_dir}/{og_id}.tmp"
    out_file = f"{out_dir}/{og_id}.faa"
    with open(tmp_file, "w") as tmp_out:
        for protein in proteins:
            tmp_out.write(f"{protein}\n")

    os.system(f"seqtk subseq {fasta_big_file} {tmp_file} > {out_file}")
    os.system(f"rm -f {tmp_file}")




def main():

    args = parse_args()

    # create out folder
    os.system(f"mkdir -p {os.path.abspath(args.out_dir)}")

    # store OGs in a list, [0]=OG, [1]=[proteins]
    ogs = [line.strip().split("\t") for line in open(args.og_file).readlines()[1:]]

    # create partial function
    func = partial(extract_sequences_OG, args.fasta_file, args.out_dir)

    pool = multiprocessing.Pool()
    pool.map(func, ogs)
    pool.close()
    pool.join()




if __name__ == "__main__":
    main()

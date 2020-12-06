'''
This script reads the filtered OGs and moves to another folder those that have
only one copy per genome, aka "easy" OGs.
'''

import argparse
import glob
import os
from Bio import SeqIO, SearchIO



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-ogs', '--ogs_faa_dir',
                               dest='ogs_faa_dir',
                               required=True,
                               help='directory with the filtered OG_xx.faa'
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help='file to write the genomic coordinates of the hits'
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    #list the OG files
    faa_files = glob.glob(f"{args.ogs_faa_dir}/*faa")

    dupl_OGs = 0
    # iterate the OG faa files and check for the presence of duplicated genomes
    for faa_file in faa_files:
        genomes = list()
        records = SeqIO.parse(faa_file, "fasta")
        for record in records:
            # extract genome and append to the list
            genome = record.id.split("|")[0]
            genomes.append(genome)

        # get unique genomes
        uniq_genomes = list(set(genomes))
        # and check for the presence of more than 1 genome
        check = True
        for genome in uniq_genomes:
            if genomes.count(genome) > 1:
                check = False


        # copy to folder if check, ie. if there is not a genome > 1 time
        if check:
            os.system(f"cp {faa_file} {args.out_dir}/{os.path.basename(faa_file)}")
        else:
            dupl_OGs += 1

    print(dupl_OGs)


if __name__ == "__main__":
    main()

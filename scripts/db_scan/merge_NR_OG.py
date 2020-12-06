'''
merges faa files from NR hmmsearch hits with the OG faa files.
Run alignment, trimming and fasttree.
'''

import argparse
import glob
import os
from Bio import Entrez, SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-nr', '--nr_faa_dir',
                               dest='nr_faa_dir',
                               required=True,
                               help='directory with filtered tblout files'
                               )
    requiredArgs.add_argument('-og', '--og_faa_dir',
                               dest='og_faa_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=False,
                               help='full seq evalue threshold'
                               )

    requiredArgs.add_argument('-msa', '--msa_dir',
                               dest='msa_dir',
                               required=False,
                               help=''
                               )
    requiredArgs.add_argument('-t', '--tree_dir',
                               dest='tree_dir',
                               required=False,
                               help=''
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    nr_files = glob.glob(f"{args.nr_faa_dir}/*faa")

    for nr_file in nr_files:
        og_id = os.path.basename(nr_file).split("nr-")[1].replace(".faa", "")

        os.system(f"cat {nr_file} {args.og_faa_dir}/{og_id}.faa > {args.out_dir}/{og_id}-NR_merged.faa")

        os.system(f"sed -i -e 's/ /_/g' -e 's/,//g' -e 's/://g' -e 's/;//g' {args.out_dir}/{og_id}-NR_merged.faa")

        os.system(f"time muscle -in {args.out_dir}/{og_id}-NR_merged.faa -out {args.msa_dir}/{og_id}-NR_merged.msa")

        os.system(f"trimal -in {args.msa_dir}/{og_id}-NR_merged.msa -out {args.msa_dir}/{og_id}-NR_merged_trim01.msa")

        os.system(f"time fasttree {args.msa_dir}/{og_id}-NR_merged_trim01.msa > {args.tree_dir}/{og_id}-NR_merged_trim01.nwk")












if __name__ == "__main__":
    main()

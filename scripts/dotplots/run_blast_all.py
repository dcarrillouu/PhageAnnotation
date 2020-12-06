'''
Runs Blast analysis to check for colinearity.

v1.0:
        - uses only the 869 genomes from Bas study. Not all of them have
          taxonomic assignment
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



    return parser.parse_args()


def main():

    args = parse_args()

    files = glob.glob(args.in_dir + "/*.fasta")
    for query in files:
        query_id = os.path.basename(query).replace(".fasta", "")
        for ref in files:
            ref_id = os.path.basename(ref).replace(".fasta", "")
            ref_db = ref.replace(".fasta", "")
            print(f"{query_id}__vs__{ref_id}")

            # blastn
            #os.system(f"/home/danielc/software/ncbi-blast-2.10.1+/bin/blastn -query {query} -db {ref_db} -out {args.out_dir}/{query_id}__vs__{ref_id}.blastn -outfmt '6 qaccver saccver length pident qlen qstart qend slen sstart send evalue bitscore' -task blastn -num_threads 70")

            # tblastx
            os.system(f"time /home/danielc/software/ncbi-blast-2.10.1+/bin/tblastx -query {query} -db {ref_db} -out {args.out_dir}/{query_id}__vs__{ref_id}.tblastx -outfmt '6 qaccver saccver length pident qlen qstart qend slen sstart send evalue bitscore' -num_threads 70")









if __name__ == "__main__":
    main()

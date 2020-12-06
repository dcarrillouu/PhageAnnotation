'''
Script to parse results from hmmsearch. If -t is provided, results of a transeq
scan are supposed as the input.

'''

import argparse
import glob
import os
import multiprocessing
from functools import partial
from Bio import SeqIO, SearchIO
import statistics



def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--hmmsearch_file',
                               dest='hmmsearch_file',
                               required=True,
                               help='plain text file with the hmmsearch results. '
                               'It can be the results of scanning proteins or transeq genomes'
                               )

    parser.add_argument('-l', '--lengths_file',
                               dest='lengths_file',
                               required=False,
                               help='file with the genome lengths "coding_statistics.txt"'
                               )

    parser.add_argument('-transeq', '--transeq',
                               dest='transeq_check',
                               required=False,
                               default=False,
                               const=True,
                               nargs='?',
                               help='marks the input as he result of a transeq search'
                               )

    parser.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=False,
                               help=''
                               )

    return parser.parse_args()


def parse_hmmsearch_proteins(hmmsearch_file):
    '''
    Parses the result of an hmmsearch using proteins as ref
    '''

    # read file
    records = SearchIO.parse(hmmsearch_file, "hmmer3-text")

    # create a dictionary to store results. Store evalue, bitscore etc for each hit
    # k=protein, v=[[hit1], [hit2], [hit3] ...]
    proteins_results = dict()

    for record in records:
        print(record.id)
        for hit in record.hits:
            print("\t",hit.id)





def main():

    args = parse_args()

    print(args.transeq_check)
    parse_hmmsearch_proteins(args.hmmsearch_file)







if __name__ == "__main__":
    main()

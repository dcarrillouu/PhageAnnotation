'''
Script that reads in the directory where tblout and domtblout files from hmmsearch,
checks if the analysis is complete by looking at the file size, reads the file and
filter it according to certain thresholds. Results are put on tabular format in the
especified directory.

* For now I look at the tblout file
'''

import argparse
import glob
import os
from Bio import SearchIO

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hmmsearch_dir',
                               dest='hmmsearch_dir',
                               required=True,
                               help='directory with tblout and domtblout files'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-e', '--evalue',
                               dest='evalue',
                               required=False,
                               type=float,
                               default=1e-30,
                               help='full seq evalue threshold'
                               )

    requiredArgs.add_argument('-s', '--score',
                               dest='bitscore',
                               required=False,
                               type=int,
                               default=100,
                               help=''
                               )

    return parser.parse_args()



def check_complete(file):
    '''
    Checks if the size of the file is 0, indicating that the analysis started
    but it still running
    '''

    if os.path.getsize(file) == 0:
        return False
    else:
        return True


def parse_tblout(tblout, evalue, bitscore, ofile):

    records = SearchIO.parse(tblout, "hmmer3-tab")

    # store hits above thresholds
    hits_threshold = list()
    for record in records:
        hits = record.hits
        for hit in hits:
            if hit.evalue <= evalue and hit.bitscore >= bitscore:
                hits_threshold.append([record.id, hit.id, hit.description, str(hit.evalue), str(hit.bitscore)])

    # check if significant hits were detected, and write them
    if hits_threshold:
        with open(ofile, "w") as fout:
            for hit in hits_threshold:
                fout.write("\t".join(hit) + "\n")



def main():

    args = parse_args()

    # get tblout files
    tblout_files = glob.glob(f"{args.hmmsearch_dir}/*.tblout")


    for tblout in tblout_files:
        # check if the analysis finished
        if check_complete(tblout):
            ofile = f"{args.out_dir}/{os.path.basename(tblout)}-filtered"
            parse_tblout(tblout, args.evalue, args.bitscore, ofile)
        else:
            print(f"{tblout} still running")


if __name__ == "__main__":
    main()

import argparse
import os
import multiprocessing
import glob
from functools import partial



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-f', '--fas_tafile',
                               dest='faa_file',
                               required=True,
                               help='directory with the .faa proteins files'
                               )
    requiredArgs.add_argument('-hmms', '--hmms_folder',
                               dest='hmms_folder',
                               required=True,
                               help='hmm db built with hmmpress'
                               )

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def run_hmmsearch(faa_file, outdir, hmm_db):
    hmm_id   = os.path.basename(hmm_db).split(".hmm")[0]
    faa_id   = os.path.basename(faa_file).split(".faa")[0]
    #outfile   = f"{os.path.abspath(outdir)}/{file_id}.out"
    tblout    = f"{os.path.abspath(outdir)}/{faa_id}-{hmm_id}.tblout"
    domtblout = f"{os.path.abspath(outdir)}/{faa_id}-{hmm_id}.domtblout"

    os.system(f"time hmmsearch --noali --tblout {tblout} --domtblout {domtblout} --cpu 1 {hmm_db} {faa_file} > /dev/null")
    print(f"{hmm_id} completed")


def main():

    args = parse_args()

    # get all hmms to scan
    hmms = glob.glob(f"{args.hmms_folder}/*.hmm")

    # def partial function
    partial_func = partial(run_hmmsearch, args.faa_file, args.out_dir)

    # parallel
    pool = multiprocessing.Pool(processes=40)
    pool.map(partial_func, hmms)
    pool.close()
    pool.join()













if __name__ == "__main__":
    main()

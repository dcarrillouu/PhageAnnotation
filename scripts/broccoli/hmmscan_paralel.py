import argparse
import os
import multiprocessing
import glob
from functools import partial



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-f', '--fasta_directory',
                               dest='f_dir',
                               required=True,
                               help='directory with the .faa proteins files'
                               )
    requiredArgs.add_argument('-db', '--hmm_db',
                               dest='hmm_db',
                               required=True,
                               help='hmm db built with hmmpress'
                               )

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def run_hmmscan(db, outdir, faa_file):
    file_id   = os.path.basename(faa_file).split(".faa")[0]
    outfile   = f"{os.path.abspath(outdir)}/{file_id}.out"
    tblout    = f"{os.path.abspath(outdir)}/{file_id}.tblout"
    domtblout = f"{os.path.abspath(outdir)}/{file_id}.domtblout"
    print(file_id)
    os.system(f"hmmscan -o {outfile} --tblout {tblout} --domtblout {domtblout} --cpu 1 {db} {faa_file}")



def main():

    args = parse_args()

    # get all faa files to scan
    faa_files = glob.glob(f"{args.f_dir}/*.faa")

    # def partial function
    partial_func = partial(run_hmmscan, args.hmm_db, args.out_dir)

    # parallel
    pool = multiprocessing.Pool(processes=None)
    pool.map(partial_func, faa_files)
    pool.close()
    pool.join()













if __name__ == "__main__":
    main()

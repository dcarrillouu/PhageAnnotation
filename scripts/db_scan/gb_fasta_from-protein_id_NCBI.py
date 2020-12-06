'''
tblout-filtered files as input, takes the protein id from the second column,
downloads the genbank first, reads the genbank and converts to fasta.
'''

import argparse
import glob
import os
from Bio import Entrez, SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--tblout_fitered_dir',
                               dest='tblout_filt_dir',
                               required=True,
                               help='directory with filtered tblout files'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-n', '--n_seqs',
                               dest='n_seqs',
                               required=False,
                               type=int,
                               default=5,
                               help='full seq evalue threshold'
                               )

    return parser.parse_args()



def get_genbank_fasta(hits, file_id, outdir):

    # first download genbank
    protein_ids = [hit[1] for hit in hits]
    #print(protein_ids, file_id)

    handle = Entrez.efetch(db="protein", id=protein_ids, rettype='gb', retmode='text')
    data = handle.read()

    with open(f"{outdir}/{file_id}.gb", "w") as fout:
        fout.write(data)


    # then read genbank and covert to fasta. Count the number of records
    # and check if it lacks some
    records = SeqIO.parse(f"{outdir}/{file_id}.gb", "gb")
    to_count = [record for record in records]
    if len(to_count) == len(hits):
        with open(f"{outdir}/{file_id}.faa", "w") as fout:
            SeqIO.write(to_count, fout, "fasta")
    else:
        print(f"not all proteins were downloaded for {file_id}")



def main():

    args = parse_args()

    Entrez.email   = 'dcarrillo.bioinf@gmail.com'
    Entrez.api_key = 'b6db7fece605d37fcabd4b93749d2e46aa09'

    # get files
    files = glob.glob(f"{args.tblout_filt_dir}/*.tblout-filtered")

    # read files. If they have less than 100 seqs, proceed with gb download
    for file in files:
        lines = [line.strip().split("\t") for line in open(file).readlines()]

        # check number of hits < 5 (default)
        if len(lines) <= args.n_seqs:
            file_id = os.path.basename(file).replace(".tblout-filtered", "")
            get_genbank_fasta(lines, file_id, args.out_dir)



if __name__ == "__main__":
    main()

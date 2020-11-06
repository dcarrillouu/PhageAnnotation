'''
This scripts takes the genome annotation table and adds OG and taxa information.
Previous stage to plotting with gggenes. The table format is similat to that I
generated to compare variant callers.

The input is the directory with all the .faa, one per genome. And the file
"table_OGs_protein_names.txt" from Broccoli step3
'''

import argparse
import multiprocessing
from functools import partial
import pandas as pd
import numpy as np
from Bio import SeqIO
import glob
import os

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_faa_directory',
                               dest='in_dir',
                               required=True,
                               help='File "table_OGs_protein_counts.txt" from step3'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-og', '--broccoli_ogs_files',
                               dest='ogs_file',
                               required=True,
                               help='"table_OGs_protein_names.txt" from Broccoli step3'
                               )

    return parser.parse_args()


def genome_to_table(ogs_dict, out_dir, faa_file):

    # read faa file
    records = SeqIO.parse(faa_file, "fasta")

    table = list()
    for record in records:
        #print(record.id)
        split = record.id.split("|")
        genome = split[0]
        start  = split[2].split("..")[0]
        end    = split[2].split("..")[1]
        strand = split[4]

        tow = [genome, start, end, strand, record.id]

        if record.id in ogs_dict[genome]:
            tow.append(ogs_dict[genome][record.id])
        else:
            tow.append("")

        table.append(tow)


    outfile = out_dir + "/" + os.path.basename(faa_file).replace(".faa", ".genes-OG")
    with open(outfile, "w") as fout:
        header = ["genome", "start", "end", "strand", "complete_id", "OG"]
        fout.write("\t".join(header) + "\n")

        for line in table:
            fout.write("\t".join(line) + "\n")



def main():

    args = parse_args()

    # Read and transpose OG file and store to dict.
    # Notice the structure of the dict:
    #   k=genome, v= dict{k=protein, v=OG_id}
    OG_annot = dict()
    ogs_matrix = pd.read_csv(args.ogs_file, header=0, sep="\t", index_col=0, low_memory=False, keep_default_na=False)
    t_ogs_matrix = np.transpose(ogs_matrix)

    genomes_id = [genome.split("_prodigal")[0] for genome in t_ogs_matrix.index.tolist()]
    proteins_list = t_ogs_matrix.values.tolist()
    ogs_id     = t_ogs_matrix.columns.tolist()

    # iterate genomes and proteins lists simultaneously
    for genome, list in zip(genomes_id, proteins_list):
        OG_annot[genome] = dict()

        # now iterate one proteins_list and OGs_id
        for proteins, OG in zip(list, ogs_id):
            if proteins != "": # if there is a protein(s) annotated for that OG...
                split = proteins.split(" ") # proteins separated by " " in the broccoli file
                for protein in split:
                    OG_annot[genome][protein] = OG


    # list *.faa files
    faa_files = glob.glob(args.in_dir + "/*prodigal*.faa")

    # create partial function
    part = partial(genome_to_table, OG_annot, args.out_dir)

    pool = multiprocessing.Pool(processes=60)
    shared_content = pool.map(part, faa_files)
    pool.close()
    pool.join()









if __name__ == "__main__":
    main()

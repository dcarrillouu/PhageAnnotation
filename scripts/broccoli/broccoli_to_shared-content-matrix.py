'''
Builds a shared-protein-content matrix with the "table_OGs_protein_counts.txt" matrix
from Broccoli, step3. This matrix has one row per OG, and one genome per column.
Counts mean the number of OG proteins found in a particular genome (column). Notice
that is not binary, value can be higher that 1 if several proteins have been found.

Among the matrix above, this script also takes as input the "n-proteins_best-coding.txt"
file to know the total proteins in a genome.
'''

import argparse
from pathlib import Path
from functools import partial
import multiprocessing
import os
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--broccoli_presabs_matrix',
                               dest='broccoli_presabs',
                               required=True,
                               help='File "table_OGs_protein_counts.txt" from step3'
                               )
    requiredArgs.add_argument('-np', '--n_proteins_genomes',
                               dest='n_prots_file',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--output_matrix',
                               dest='out_matrix',
                               required=True,
                               help=''
                               )

    return parser.parse_args()

def compute_genome_shared(genomes_list, ogs_dict, n_prots_list, genome_query):
    # initialize the list for the genome that will contain the shared values
    shared_content = list()

    # iterate the genomes
    for genome_target in genomes_list:
        shared = 0
        # iterate the OGs list
        for i in range(0, len(ogs_dict[genome_query])):
            # if both of the genomes have that ortholog...
            if ogs_dict[genome_query][i] != 0 and ogs_dict[genome_target][i] != 0:
                # ... sum the n_prots of the query genome to shared
                shared += int(ogs_dict[genome_query][i])

        ### !! For now I left the matrix asymmetric.
        # lowest n_protein of the two genomes in the denominator
        #lowest = min(n_prots_list[genome_query], n_prots_list[genome_target])
        #shared_perc = float(shared/lowest)
        shared_perc = float(shared/n_prots_list[genome_query])
        shared_content.append("{:.2f}".format(shared_perc))


    shared_content.insert(0, genome_query)
    #print(shared_content)
    return shared_content


def main():
    args = parse_args()

    # read presabs_matrix. Rrecall that OGs are in the row and genomes in the columns
    # it will be easier if we transponse the matrix
    brocco_matrix = pd.read_csv(args.broccoli_presabs, header=0, sep="\t", index_col=0)
    t_brocco_matrix = np.transpose(brocco_matrix)

    # create a dict to store k=genome values=[protein counts]
    # first extract genomes and values to lists
    # it assumes that is PRODIGAL annotation
    genomes_id = [genome.split("_prodigal")[0] for genome in t_brocco_matrix.index.tolist()]
    counts  = t_brocco_matrix.values.tolist()
    # create the dictionary
    genomes_ogs = {genome:counts_list for genome, counts_list in zip(genomes_id, counts)}

    #print(genomes_ogs["OGQH01000028"])

    # store n_prots for each genome in a dict
    n_prots = {line.split("\t")[0]:int(line.strip().split("\t")[2]) for line in open(args.n_prots_file).readlines()}

    # prepare the partial function
    part = partial(compute_genome_shared, genomes_id, genomes_ogs, n_prots)

    pool = multiprocessing.Pool(processes=60)
    shared_content = pool.map(part, genomes_id)
    pool.close()
    pool.join()

    # prepare header and out table
    genomes_id.insert(0, "")
    with open(args.out_matrix, "w") as fout:
        fout.write("\t".join(genomes_id) + "\n")
        for line in shared_content:
            fout.write("\t".join(line) + "\n")

if __name__ == '__main__':
    main()

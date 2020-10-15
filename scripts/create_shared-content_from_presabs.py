'''
Creates the shared-protein-content matrix from the non-binary presabs matrix
generated with "create_presabs_from_mcl.py"
'''

import argparse
from pathlib import Path
from functools import partial
import multiprocessing
import os


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--presabs_matrix',
                               dest='presabs_file',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-np', '--n_proteins_genomes',
                               dest='n_prots_file',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--output_matrix',
                               dest='out_matrix',
                               type=lambda p: Path(p).resolve(),
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
            if ogs_dict[genome_query][i] != "0" and ogs_dict[genome_target][i] != "0":
                # ... sum the n_prots of the query genome to shared
                shared += int(ogs_dict[genome_query][i])

        # lowest n_protein of the two genomes in the denominator
        lowest = min(n_prots_list[genome_query], n_prots_list[genome_target])
        shared_perc = float(shared/lowest)
        shared_content.append("{:.2f}".format(shared_perc))


    shared_content.insert(0, genome_query)
    #print(shared_content)
    return shared_content




def main():
    args = parse_args()

    # read presabs_matrix. Discard first line (OGs header)
    lines = [line.strip().split("\t") for line in open(args.presabs_file).readlines()[1:]]

    # store OGs counts for each genome in a dict.
    genomes_ogs = {line[0]:line[1:] for line in lines}

    # create a genomes list
    genomes_id = [line[0] for line in lines]
    #print(genomes_id)

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

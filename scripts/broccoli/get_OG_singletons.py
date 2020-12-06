'''
Takes the file where each gene is annotated with an OG, looks for genes
that were not included in any OG, and search them in the hmmscan results to see
if there are significant hits to OGs.
'''

import argparse
import multiprocessing
from functools import partial
from Bio import SearchIO
import glob
import os



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-g', '--genomes-OG_tables_dir',
                               dest='tables_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-hmm', '--hmmscan_dir',
                               dest='hmmscan_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def inspect_genes(hmmscan_dir, genome_OG_table):

    # read genes_OG
    genes = [line.strip().split("\t") for line in open(genome_OG_table, "r").readlines()[1:]]

    # identify lines without OG at the end (list shorter than 6)
    check = [gene[-1] for gene in genes if len(gene) == 5]

    to_return = list()
    # if genes without OG are identified...
    if check:
        genome_id = os.path.basename(genome_OG_table).split(".genes-OG")[0]

        # read hmmscan results and store to dict
        hmmscan_file = f"{hmmscan_dir}/{genome_id}.tblout"
        records = SearchIO.parse(hmmscan_file, "hmmer3-tab")

        # store only evalue and bitscore (for now)
        results = dict()
        for record in records:
            results[record.id] = list()
            hits = record.hits
            for hit in hits:
                results[record.id] = [hit.id, str(hit.evalue), str(hit.bitscore)]


        # iterate through unassigned genes and provide with hmmscan results if there was a hit

        for gene in check:
            if gene in results:
                to_return.append([gene, results[gene][0], results[gene][1], results[gene][2]])
            else:
                to_return.append([gene, "", "", ""])

    return to_return





def main():

    args = parse_args()


    # get genome_OG tables.
    tables = glob.glob(f"{args.tables_dir}/*genes-OG")

    partial_func = partial(inspect_genes, args.hmmscan_dir)

    # parallel
    pool = multiprocessing.Pool()
    singletons = pool.map(partial_func, tables)
    pool.close()
    pool.join()

    tow = list()
    for genome in singletons:
        if genome: # if there are singletons...
            for singleton in genome:
                tow.append(singleton)

    tow.sort(key=lambda singleton: singleton[0])

    with open(args.out_file, "w") as fout:
        for gene in tow:
            fout.write("\t".join(gene) + "\n")







if __name__ == "__main__":
    main()

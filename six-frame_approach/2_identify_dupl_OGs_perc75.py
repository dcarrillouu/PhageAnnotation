'''
This script reads the filtered OG_xx.faa files, identify OGs with more than one protein
per genome, and substracts the sequences with length in percentile 75-100.
'''

import argparse
import glob
import os
import numpy as np
from functools import partial
import multiprocessing
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-of', '--ogs_faa_dir',
                               dest='ogs_faa_dir',
                               required=True,
                               help='directory with the filtered OG_xx.faa files'
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def compute_75_percentile(out_dir, faa_file):
    '''
    Computes the median of the proteins' length and stores those within a certain range.
    It also writes a table with statistics about the number of excluded sequences. Format:
     [0]        [1]          [2]         [3]       [4]       [5]
    OG_id  total_proteins   75th perc   range   included   excluded
    '''

    OG_id = os.path.basename(faa_file).split(".faa")[0]
    records = SeqIO.parse(faa_file, "fasta")

    # get all the genome_ids present in the file, do the uniq and check if any is more than one time
    genomes = [record.id.split("|")[0] for record in records]
    uniq = list(set(genomes))
    check = False
    for genome in uniq:
        if genomes.count(genome) > 1:
            check = True

    # if a genome has been seen more than one time, proceed to percentile estimation
    if check:
        records = SeqIO.parse(faa_file, "fasta")
        seq_lengths = [len(record.seq) for record in records]
        #print(seq_lengths)


        percentiles = np.percentile(seq_lengths, [0, 25, 50, 75, 100])
        min, max = percentiles[3], percentiles[4]

        in_range  = list()
        out_range = list()
        in_range_ids  = list()
        out_range_ids = list()

        for record in SeqIO.parse(faa_file, "fasta"):
            if min <= len(record.seq) <= max:
                in_range.append(record)
                in_range_ids.append(record.id)
            else:
                out_range.append(record)
                out_range_ids.append(record.id)

        total = len(in_range) + len(out_range)
        #print(OG_id, total, "perc75")
        fraction_incl = "{:.2f}".format(len(in_range)/total)
        fraction_excl = "{:.2f}".format(len(out_range)/total)

        # check the number of sequences in_range. If it's only one, move the min value to be the median (percentiles[2])
        if len(in_range) > 1:
            with open(f"{out_dir}/{OG_id}_perc75.faa", "w") as fout:
                SeqIO.write(in_range, fout, "fasta")


            stats_tow = [OG_id, str(total), str(percentiles[3]), f"{min}-{max}", str(len(in_range)), str(len(out_range)), fraction_incl, fraction_excl, ",".join(in_range_ids),",".join(out_range_ids)]
            return stats_tow

        # if there was only one sequence within range...
        else:
            min, max = percentiles[2], percentiles[4]

            in_range  = list()
            out_range = list()
            in_range_ids  = list()
            out_range_ids = list()


            for record in SeqIO.parse(faa_file, "fasta"):
                if min <= len(record.seq) <= max:
                    in_range.append(record)
                    in_range_ids.append(record.id)
                else:
                    out_range.append(record)
                    out_range_ids.append(record.id)
                    
            total = len(in_range) + len(out_range)
            #print(OG_id, total, "perc50")
            fraction_incl = "{:.2f}".format(len(in_range)/total)
            fraction_excl = "{:.2f}".format(len(out_range)/total)

            # check the number of sequences in_range. If it's only one, move the min value to be the median (percentiles[2])

            with open(f"{out_dir}/{OG_id}_perc50.faa", "w") as fout:
                SeqIO.write(in_range, fout, "fasta")

            stats_tow = [OG_id, str(total), str(percentiles[3]), f"{min}-{max}", str(len(in_range)), str(len(out_range)), fraction_incl, fraction_excl, ",".join(in_range_ids),",".join(out_range_ids)]
            return stats_tow
    else:
        return ""




def main():

    args = parse_args()

    # get faa files
    faa_files = glob.glob(f"{args.ogs_faa_dir}/*faa")

    # partial function
    partial_func = partial(compute_75_percentile, args.out_dir)

    pool = multiprocessing.Pool(processes=None)
    dupl_OGs_stats = pool.map(partial_func, faa_files)
    pool.close()
    pool.join()

    while "" in dupl_OGs_stats:
        dupl_OGs_stats.remove("")

    with open(f"{args.out_dir}/perc75_stats.txt", "w") as fout:
        for line in dupl_OGs_stats:
            fout.write("\t".join(line) + "\n")







if __name__ == "__main__":
    main()

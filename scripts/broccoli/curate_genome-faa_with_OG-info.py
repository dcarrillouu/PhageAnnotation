'''
Script to curate faa files for each genome, mostly by merging split genes.
It uses the OG information to know genes ORFs need to be merged.

    1)  Checks in the original prodigal.gff if the gene has both start-stop codons.
        If not, discard the gene.
'''

import argparse
import glob
import os
import multiprocessing
from functools import partial
from Bio import SeqIO



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--genome_tables',
                               dest='tables_dir',
                               required=True,
                               help='dir with genome_tables where each gene is '
                               'annotated with its OG'
                               )
    requiredArgs.add_argument('-p', '--prodigal_dir',
                               dest='prodigal_dir',
                               required=True,
                               help='directory with the prodigal annotation. I use '
                               '"2_ORF_calling/final_annotation/0_fna_faa_gff"'
                               )

    # requiredArgs.add_argument('-mcl', '--mcl_file',
    #                            dest='mcl_file',
    #                            required=True,
    #                            help='tab file from MCL, one cluster per row'
    #                            )
    # requiredArgs.add_argument('-stats', '--ogs_stats_file',
    #                            dest='ogs_stats_file',
    #                            required=True,
    #                            help='file "statistics_per_OG.txt"'
    #                            )


    return parser.parse_args()


def check_completeness(genome_table, prodigal_dir):
    '''
    Goes to the gff and check for the completeness of the gene. This step is mostly
    to check if I have to discard the first/last gene.
    '''

    # get genome_id, which already conteins caller info
    genome = os.path.basename(genome_table).split(".genes-OG-profiles")[0]

    # read gff file.
    gff_file = f"{prodigal_dir}/{genome}.gff"
    lines = [line for line in open(gff_file).readlines() if not line.startswith("#")]
    # iterate the lines looking for incomplete genes
    incomplete = list()
    for line in lines:
        # if I find an incomplete gene, split the line and store the gene
        if ";partial=10;" in line or ";partial=01;" in line:
            gene = line.split("\tID=")[1].split(";")[0]
            incomplete.append(gene)


    # read genome table and store to dict, k=short_gene_id (ie. 1_40), v=whole_line
    table = dict()

    # raw table is like this:
    # genome    start   end    strand    complete_id     yutin_profiles   OG
    # Insert a new column:
    # genome    start   end    strand  COMPLETENESS   complete_id     yutin_profiles   OG

    lines = [line.split("\t") for line in open(genome_table).readlines()]
    for line in lines[1:]:
        line[-1] = line[-1].replace("\n", "")
        short_gene_id = line[4].split("|")[-1]
        # check if the gene is incomplete
        if short_gene_id in incomplete:
            line.insert(4, "partial")
        else:
            line.insert(4, "")

        table[short_gene_id] = line

    return table





def check_duplicated_OGs(genome_dict):
    '''
    Looks for duplicated OG within the genome.
    '''


    for i in range(1, len(genome_dict)-3):
        if genome_dict[f"1_{i}"][-1] == genome_dict[f"1_{i+1}"][-1] and genome_dict[f"1_{i}"][-1] not in ["OG_75", "OG_310", ""]:
            print(genome_dict[f"1_{i}"][5:])
            print(genome_dict[f"1_{i+1}"][5:])
            print(genome_dict[f"1_{i+2}"][5:])















def curate(genome_table, prodigal_dir):

    genome_dict = check_completeness(genome_table, prodigal_dir)

    check_duplicated_OGs(genome_dict)








def main():

    args = parse_args()

    # get genome tables. I work with that tables and not with the faa files directly
    tables = glob.glob(f"{args.tables_dir}/*.genes-OG-profiles")

    for table in tables:
        curate(table, args.prodigal_dir)










if __name__ == '__main__':
    main()

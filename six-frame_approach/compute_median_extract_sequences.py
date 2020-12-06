'''

'''

import argparse
import glob
import os
import multiprocessing
import numpy as np
from functools import partial
from Bio import SeqIO
import statistics



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-b', '--broccoli_table',
                               dest='broccoli_table',
                               required=True,
                               help='file "table_OGs_protein_names.txt" from broccoli. '
                               'The script takes the protein names per OG from this file.'
                               )
    requiredArgs.add_argument('-f', '--faa_fle',
                               dest='faa_file',
                               required=True,
                               help='faa file with all the proteins, "all_final_proteins.faa"'
                               )

    requiredArgs.add_argument('-c', '--chimeric_file',
                               dest='chimeric_file',
                               required=True,
                               help='file with chimeric proteins from broccoli, "chimeric_proteins.txt"'
                               )

    requiredArgs.add_argument('-g', '--gff_dir',
                               dest='gff_dir',
                               required=True,
                               help='directory with the GFF prodigal annotation. I use '
                               '"2_ORF_calling/final_annotation/0_fna_faa_gff"'
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )


    return parser.parse_args()

def read_chimeric_file(chimeric_file):
    '''
    Reads the "chimeric_proteins.txt" file and returns a list with the chimeric
    protein names
    '''
    # names on column [1]. Discard header
    chim_names = [line.strip().split("\t")[1] for line in open(chimeric_file).readlines()[1:]]
    return chim_names

def faa_to_dict(faa_file):
    '''
    Parses the faa file with all the proteins and store in a dict, k=record.id, v=record
    '''
    # parse the faa_file
    records = SeqIO.parse(faa_file, "fasta")

    # store to dict
    protein_records = dict()
    for record in records:
        protein_records[record.id] = record

    return protein_records

def get_partial_proteins(protein_records, gff_dir):
    '''
    Takes the proteins records as input. From that, extract the genome id with
    the caller information. Then go to the gff and store its information.
    It stores the GENE SHORT ID (ie. 1_149) if it is partial.
    '''
    genome_caller = dict()
    for protein in protein_records:
        split = protein.split("|")
        genome_caller[split[0]] = split[1]

    partial = dict()
    for genome, caller in genome_caller.items():
        gff_file = f"{gff_dir}/{genome}_{caller}.gff"

        lines = [line for line in open(gff_file).readlines() if not line.startswith("#")]
        # iterate the lines looking for incomplete genes
        partial[genome] = list()
        for line in lines:
            # if I find an incomplete gene, split the line and store the gene short id
            if ";partial=10;" in line or ";partial=01;" in line:
                gene = line.split("\tID=")[1].split(";")[0]
                #print("\t", gene, line)
                partial[genome].append(gene)
        #print(partial[genome])

    return partial

def broccoli_to_dict(broccoli_table, chimeric_proteins, partial_proteins):
    '''
    Converts the broccoli "table_OGs_protein_names.txt" file to dict, k=OG_id, v=[protein_names]
    Don't consider chimeric and partial sequences.
    It returns a dict ONLY with the OGs that appear more than one time in the genome (after removing chimeric and partial seqs of course)
    '''
    # discard header line with the genome information
    lines = [line.strip().split("\t") for line in open(broccoli_table).readlines()[1:]]

    #
    OGs_proteins = dict()
    OGs_proteins_dupl = dict()
    for line in lines:
        # remove empty columns: genomes for which the OG was not detected
        while "" in line:
            line.remove("")
        # when 2 proteins are found in a genome, they are in the same column separated by a white space.
        # in those cases I have to split by " " and append to another list (flat_list)
        flat_list = list()
        for genome_column in line[1:]: # discard first column, is the OG_id
            proteins = genome_column.split(" ")

            for protein in proteins:
                flat_list.append(protein)

        # check if chimeric and remove if it is
        for protein in flat_list:
            if protein in chimeric_proteins:
                flat_list.remove(protein)


        # FOR WHATEVER THE REASON, IT WAS NOT REMOVING CORRECTLY THE PARTIAL PROTEINS
        # BUT IT DOES WHEN I ADD OTHER LOOP WITH THE SAME CODE. REALLY REALLY WEIRD
        # remove partial proteins. Recall that id of partial genes is stored as the SHORT ID, so first I have to obtain the short id
        for protein in flat_list:
            splited = protein.split("|")
            if splited[-1] in partial_proteins[splited[0]]:
                #print(split[0], partial_proteins[split[0]], split[-1])
                flat_list.remove(protein)


        # remove partial proteins. Recall that id of partial genes is stored as the SHORT ID, so first I have to obtain the short id
        for protein in flat_list:
            splited = protein.split("|")
            if splited[-1] in partial_proteins[splited[0]]:
                #print(protein)
                #print(split[0], partial_proteins[split[0]], split[-1])
                flat_list.remove(protein)


        OGs_proteins[line[0]] = flat_list


        # check if a genomes is more than one time in the OG. ie, if there are
        # more than one protein in the genome annotated with that OG
        genomes = [protein.split("|")[0] for protein in flat_list]
        uniq = list(set(genomes))
        for genome in uniq:
            if genomes.count(genome) > 1:
                # associate proteins with OG_id
                OGs_proteins_dupl[line[0]] = flat_list





    return OGs_proteins, OGs_proteins_dupl

def get_fasta_sequences(dictionary, fasta_file, outdir):
    '''
    Uses seqtk to extract given protein ids per og and write to file
    '''
    for OG, proteins in dictionary.items():
        #print(OG)
        tmpfile = f"{outdir}/{OG}.ids"
        outfile = f"{outdir}/{OG}.faa"

        with open(tmpfile, "w") as tmp:
            tmp.write("\n".join(proteins))

        os.system(f"seqtk subseq {fasta_file} {tmpfile} > {outfile}")

def calculate_median_statistics(OGs_proteins, protein_records, out_dir):
    '''
    Computes the median of the proteins' length and stores those within a certain range.
    It also writes a table with statistics about the number of excluded sequences. Format:
     [0]        [1]          [2]      [3]       [4]       [5]
    OG_id  total_proteins   median   range   included   excluded
    '''

    # initiate list to store statistics and lately write to file
    stats_tow = list()

    in_range  = dict()
    out_range = dict()

    for OG, proteins in OGs_proteins.items():
        seq_lengths = [len(protein_records[protein].seq) for protein in proteins]

        # calculate the median length
        median = statistics.median(seq_lengths)

        # calculate the inclusion range from the median.
        # 10% of deviation
        min, max = median-(median*0.1), median+(median*0.1)

        # check wether proteins are within median range
        in_range[OG]  = list()
        out_range[OG] = list()
        for protein in proteins:
            if min <= len(protein_records[protein].seq) <= max:
                in_range[OG].append(protein)
            else:
                out_range[OG].append(protein)


        if len(in_range[OG]) > 0:
            fraction_incl = "{:.2f}".format(len(in_range[OG])/len(proteins))
            fraction_excl = "{:.2f}".format(len(out_range[OG])/len(proteins))

            stats_tow.append([OG, str(len(proteins)), str(median), f"{min}-{max}", str(len(in_range[OG])), str(len(out_range[OG])), fraction_incl, fraction_excl, "median"])
        else:
            # calculate the percentiles
            percentiles = np.percentile(seq_lengths, [0, 25, 50, 75, 100])

            # calculate the inclusion range from the 75th percentile.
            # 10% of deviation
            min, max = percentiles[3]-(percentiles[3]*0.1), percentiles[3]+(percentiles[3]*0.1)

            # check wether proteins are within median range
            in_range[OG]  = list()
            out_range[OG] = list()
            for protein in proteins:
                if min <= len(protein_records[protein].seq) <= max:
                    in_range[OG].append(protein)
                else:
                    out_range[OG].append(protein)

            fraction_incl = "{:.2f}".format(len(in_range[OG])/len(proteins))
            fraction_excl = "{:.2f}".format(len(out_range[OG])/len(proteins))

            stats_tow.append([OG, str(len(proteins)), str(percentiles[3]), f"{min}-{max}", str(len(in_range[OG])), str(len(out_range[OG])), fraction_incl, fraction_excl, "percentile"])


    # write stats to file
    header = ["OG_id", "total_proteins", "median/perc.", "range", "included", "excluded", "\n"]
    with open(f"{out_dir}/median-percentile_statistics.txt", "w") as fout:
        fout.write("\t".join(header))

        for line in stats_tow:
            fout.write("\t".join(line) + "\n")

    return in_range

def compute_75_percentile(OGs_proteins, protein_records, out_dir):
    '''
    Computes the median of the proteins' length and stores those within a certain range.
    It also writes a table with statistics about the number of excluded sequences. Format:
     [0]        [1]          [2]      [3]       [4]       [5]
    OG_id  total_proteins   median   range   included   excluded
    '''

    # initiate list to store statistics and lately write to file
    stats_tow = list()

    in_range  = dict()
    out_range = dict()

    for OG, proteins in OGs_proteins.items():
        seq_lengths = [len(protein_records[protein].seq) for protein in proteins]

        # calculate the percentiles
        percentiles = np.percentile(seq_lengths, [0, 25, 50, 75, 100])

        # calculate the inclusion range from the 75th percentile.
        # 10% of deviation
        #min, max = percentiles[3]-(percentiles[3]*0.1), percentiles[3]+(percentiles[3]*0.1)
        min, max = percentiles[3], percentiles[4]

        # check wether proteins are within median range
        in_range[OG]  = list()
        out_range[OG] = list()
        for protein in proteins:
            if min <= len(protein_records[protein].seq) <= max:
                in_range[OG].append(protein)
            else:
                out_range[OG].append(protein)

        fraction_incl = "{:.2f}".format(len(in_range[OG])/len(proteins))
        fraction_excl = "{:.2f}".format(len(out_range[OG])/len(proteins))


        if len(in_range[OG]) > 1:
            stats_tow.append([OG, str(len(proteins)), str(percentiles[3]), f"{min}-{max}", str(len(in_range[OG])), str(len(out_range[OG])), fraction_incl, fraction_excl, ",".join(in_range[OG]),",".join(out_range[OG]) ])
        else:
            min, max = percentiles[2], percentiles[4]

            # check wether proteins are within median range
            in_range[OG]  = list()
            out_range[OG] = list()
            for protein in proteins:
                if min <= len(protein_records[protein].seq) <= max:
                    in_range[OG].append(protein)
                else:
                    out_range[OG].append(protein)

            fraction_incl = "{:.2f}".format(len(in_range[OG])/len(proteins))
            fraction_excl = "{:.2f}".format(len(out_range[OG])/len(proteins))

            stats_tow.append([OG, str(len(proteins)), str(percentiles[3]), f"{min}-{max}", str(len(in_range[OG])), str(len(out_range[OG])), fraction_incl, fraction_excl, ",".join(in_range[OG]),",".join(out_range[OG]) ])


    # write stats to file
    header = ["OG_id", "total_proteins", "75th perc.", "range", "included", "excluded", "\n"]
    with open(f"{out_dir}/perc75-100_statistics.txt", "w") as fout:
        fout.write("\t".join(header))

        for line in stats_tow:
            fout.write("\t".join(line) + "\n")

    return in_range




def main():

    args = parse_args()

    #
    chimeric_proteins = read_chimeric_file(args.chimeric_file)

    #
    protein_records = faa_to_dict(args.faa_file)

    #
    partial_proteins = get_partial_proteins(protein_records, args.gff_dir)

    #
    OGs_proteins, OGs_proteins_dupl = broccoli_to_dict(args.broccoli_table, chimeric_proteins, partial_proteins)

    #
    #outdir_all_OGs_filtered = f"{args.out_dir}/1_filtered_proteins/0_all_OG_fasta/"
    #get_fasta_sequences(OGs_proteins, args.faa_file, outdir_all_OGs_filtered)

    #
    #calculate_median_statistics(OGs_proteins, protein_records, args.out_dir)

    #
    in_range_seqs = compute_75_percentile(OGs_proteins_dupl, protein_records, args.out_dir)

    #
    outdir_dupl_OGs_filtered = f"{args.out_dir}/1_filtered_proteins/0_fasta/1_dupl_OGs_fasta/"
    get_fasta_sequences(in_range_seqs, args.faa_file, outdir_dupl_OGs_filtered)


if __name__ == '__main__':
    main()

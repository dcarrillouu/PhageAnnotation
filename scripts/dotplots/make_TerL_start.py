'''
Reorders the genome to make it starts from the TerL gene.
It takes the genome tables annotated with Broccoli families, looks for the
family specified by the user, and makes the genomes start from there.

Notice that if the protein is in the negative strand, the results will be
reverse complemented
'''


import argparse
import os
import glob
from functools import partial
import multiprocessing
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--tables_dir',
                               dest='tables_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-g', '--genomes_dir',
                               dest='genomes_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )



    return parser.parse_args()



def main():

    args = parse_args()

    # read and store genome tables with OG information. k=genome, v=[lines]
    tables = glob.glob(args.tables_dir + "/*.genes-OG")

    genomes_OG = dict()
    for table in tables:
        genome = os.path.basename(table).split("_prodigal")[0]
        lines = [line.strip().split("\t") for line in open(table).readlines()[1:]]
        genomes_OG[genome] = lines

    # detect the OG_60 and its start.
    # FOR NOW I GET ONLY THOSE WITH ONE TERMINASE
    terminase_position = dict()
    for genome, genes in genomes_OG.items():
        # count number of times OG_60 is in the genome
        ogs = [gene[-1] for gene in genes]

        # if it is more than 1 time...
        if ogs.count("OG_60") > 1:
            # the thing here is to see where the gene starts
            terl_genes = [gene for gene in genes if gene[-1] == "OG_60"]
            #print(terl_genes)
            strand = list(set([gene[3] for gene in terl_genes])) # collapse the strand info for all the genes
            # check that all the terL orfs are in the same strand. 2 genomes with more than 1. don't do nothing with them for now

            if len(strand) == 1 and strand[0] == "+": # forward strand
                for gene in terl_genes:
                    if gene[-2].split("|")[-1].strip() not in ["1_1", "1_2", "1_3"]: # discard the gene if it's at the beginning of the genome, it's a splited gene
                        terminase_position[genome] = [int(gene[1]), int(gene[2]), gene[3]] # retain the first gene that is not the 1_1, 1_2 or 1_3...
                        break # ...and break the loop to keep it and don't overwrite

                # print(terl_genes)
                # print(genome, terminase_position[genome])
                # print()

            if len(strand) == 1 and strand[0] == "-": # reverse strand
                #print(terl_genes)
                for gene in reversed(terl_genes): # start from the last terl gene in the genome
                    if gene != genes[-1] and gene != genes[-2] and gene != genes[-3]: # discard the gene if it's at the end of the genome, it's a splited gene
                        terminase_position[genome] = [int(gene[1]), int(gene[2]), gene[3]] # retain the first gene that is not the 1_1, 1_2 or 1_3...
                        break # ...and break the loop to keep it and don't overwrite
                print(terl_genes)
                print(genome, terminase_position[genome])
                print()

        # if it is 1 time...
        elif ogs.count("OG_60") == 1:
            for gene in genes:
                if gene[-1] == "OG_60":
                    # [start, end, strand]
                    terminase_position[genome] = [int(gene[1]), int(gene[2]), gene[3]]


    # process genomes with ther TerL information
    fastas = glob.glob(args.genomes_dir + "/*.fasta")
    fastas_id = {os.path.basename(fasta).replace(".fasta", ""):fasta  for fasta in fastas}
    #print(fastas_id)

    for genome_id, terl in terminase_position.items():
        ofile = f"{args.out_dir}/{genome_id}.fasta"
        record = SeqIO.read(fastas_id[genome_id], "fasta")

        if terl[2] == "+": # if is the forward strand
            seq1 = record.seq[terl[0]:]
            seq2 = record.seq[:terl[0]]
            seq = seq1 + seq2
            tow = SeqRecord(seq, id=record.id, description='')

            SeqIO.write(tow,ofile, "fasta")

        else: # if is the reverse strand
            seq1 = record.seq[:terl[1] + 1].reverse_complement()
            seq2 = record.seq[terl[1] + 1:].reverse_complement()
            seq = seq1 + seq2
            tow = SeqRecord(seq, id=record.id, description='')

            SeqIO.write(tow,ofile, "fasta")


if __name__ == "__main__":
    main()

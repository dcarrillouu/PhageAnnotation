'''
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
                               help='directory with final annotation .faa file for each genome'
                               )
    requiredArgs.add_argument('-gff', '--input_gff_directory',
                               dest='gff_dir',
                               required=True,
                               help='directory with the GFF prodigal annotation. I use '
                               '"2_ORF_calling/final_annotation/0_fna_faa_gff"'
                               )

    requiredArgs.add_argument('-og', '--broccoli_ogs_files',
                               dest='ogs_file',
                               required=True,
                               help='"table_OGs_protein_names.txt" from Broccoli step3'
                               )

    requiredArgs.add_argument('-ch', '--chimeric_proteins_file',
                               dest='chimeric_file',
                               required=True,
                               help='"chimeric_proteins.txt" from Broccoli step3'
                               )

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help='directory to put the genome tables, one file per genome'
                               )


    return parser.parse_args()


def OGs_to_dict(ogs_file):
    '''
    Reads the "table_OGs_protein_names.txt" from Broccoli step3 and store the protein names
    in a dict with key=genome
    '''

    # Read and transpose OG file and store to dict.
    # Notice the structure of the dict:
    #   k=genome, v= dict{k=protein, v=OG_id}
    OG_annot = dict()
    ogs_matrix = pd.read_csv(ogs_file, header=0, sep="\t", index_col=0, low_memory=False, keep_default_na=False)
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
                # some proteins are chimeric. Account for that
                for protein in split:
                    if protein in OG_annot[genome]:
                        OG_annot[genome][protein] += f",{OG}"
                    else:
                        OG_annot[genome][protein] = OG

    return OG_annot

def get_partial_proteins(gff_dir):
    '''
    It reads the gff files in the provided directory and returs a dictionary with the
    partial genes, if they is any. k=genome, v=[short_id of partial genes]
    '''

    # get the files
    gff_files = glob.glob(f"{gff_dir}/*.gff")

    partial_proteins = dict()
    for gff_file in gff_files:
        genome = os.path.basename(gff_file).split("_prodigal")[0]
        partial_proteins[genome] = list()

        # read the gff file of the genome
        lines = [line for line in open(gff_file).readlines() if not line.startswith("#")]

        for line in lines:
            # if I find an incomplete gene, split the line and store the gene short id
            if ";partial=10;" in line or ";partial=01;" in line:
                gene = line.split("\tID=")[1].split(";")[0]
                partial_proteins[genome].append(gene)

    return partial_proteins

def parse_chimeric_file(chimeric_file):
    '''
    Reads the "chimeric_proteins.txt" file and returns a list with the chimeric
    protein names
    '''
    # names on column [1]. Discard header
    chimeric_proteins = [line.strip().split("\t")[1] for line in open(chimeric_file).readlines()[1:]]

    return chimeric_proteins

def create_table(ogs_dict, chimeric_p, partial_p, yutin_annot, pfam_annot, nr_annot, out_dir, faa_file):
    '''
    Reads the faa file with ALL the proteins of the genome. Then it writes all
    associated information:

    Genome  Protein_long_id  Protein_start  Protein_end  Strand  Partial?  Chimeric?
      0          1                 2            3          4        5          6

    Protein_short_id  Yutin  PFAM  Blastp_NR  Assigned_OG
          7             8     9       10         11
    '''

    # read faa file
    records = SeqIO.parse(faa_file, "fasta")

    # iterate sequences and add to table
    table = list()
    for record in records:
        # get available information from the protein_id (protein header)
        split = record.id.split("|")
        genome = split[0]
        start  = split[2].split("..")[0]
        end    = split[2].split("..")[1]
        strand = split[4]
        short  = split[-1]

        # check for partialness
        if short in partial_p[genome]:
            partialness = "True"
        else:
            partialness = "False"

        # check for chimeric
        if record.id in chimeric_p:
            chimeric = "True"
        else:
            chimeric = "False"

        # check for yutin
        if record.id in yutin_annot:
            yutin = yutin_annot[record.id]
        else:
            yutin = ""

        # check for pfam
        if record.id in pfam_annot:
            pfam = pfam_annot[record.id]
        else:
            pfam = ""

        # check for nr
        if record.id in nr_annot:
            nr= nr_annot[record.id]
        else:
            nr = ""

        # check for OG_id
        if record.id in ogs_dict[genome]:
            og_id = ogs_dict[genome][record.id]
        else:
            og_id = ""


        # merge all the information and add to table
        tow = [genome, record.id, start, end, strand, partialness, chimeric, short, yutin, pfam, nr, og_id]
        table.append(tow)


    outfile = f"{out_dir}/" + os.path.basename(faa_file).replace(".faa", ".table")
    with open(outfile, "w") as fout:
        header = ["genome", "protein_id", "start", "end", "strand", "partial", "chimeric", "short_id", "yutin", "pfam", "nr", "OG"]
        fout.write("\t".join(header) + "\n")

        for line in table:
            fout.write("\t".join(line) + "\n")



def main():

    args = parse_args()

    # read broccoli table relating each protein to an OG
    OG_annot = OGs_to_dict(args.ogs_file)

    # get partial proteins info from gff files
    partial_proteins = get_partial_proteins(args.gff_dir)

    # get chimeric proteins
    chimeric_proteins = parse_chimeric_file(args.chimeric_file)

    # Now get faa files. I use these because they contain all the proteins, including
    # those that were not assigned to an OG.
    faa_files = glob.glob(args.in_dir + "/*prodigal*.faa")

    # I DONT HAVE YUTIN, PFAM AND NR ANNOT, SO FOR NOW I CREATE MOCK DICTIONARIES
    yutin_annot = dict()
    pfam_annot = dict()
    nr_annot = dict()

    # create partial function
    part = partial(create_table, OG_annot, chimeric_proteins, partial_proteins, yutin_annot, pfam_annot, nr_annot, args.out_dir)

    pool = multiprocessing.Pool(processes=10)
    shared_content = pool.map(part, faa_files)
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()

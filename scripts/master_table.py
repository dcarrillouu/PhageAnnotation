'''
Script to merge in table taxonomic, OG and functional information. It intends to
show tendencies within a OG. Also, it is useful to have all this information
together.

It takes as input the "OG proteins" faa files which contains all the proteins
assigned to an OG. From that I obtain the relation protein-OG.
It also takes the "contig_clustering.csv" file from Bas to obtain the taxonomy.
I take the annotation that Yutin obtained from her prot2profile and prot2domains
files.


For now:
    - [0]: Protein ID
    - [1]: Genome
    - [2]: Family
    - [3]: Subfamily
    - [4]: Protein length
    - [5]: OG
    - [6]: number of OG proteins in the genome
    #- []: Yutin's paper profile annotation
    #- []: Yutin's paper domain annotation
'''
import argparse
import os
import glob
from Bio import SeqIO
import pandas as pd



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-faa', '--OGs_faa',
                               dest='OGs_faa',
                               required=True,
                               help='directory with the OGs .faa files'
                               )
    requiredArgs.add_argument('-taxa', '--taxa_file',
                               dest='taxa_file',
                               required=True,
                               help='Bas taxa file "contigs_clustering.txt"'
                               )
    requiredArgs.add_argument('-n', '--n_og_file',
                               dest='n_og_file',
                               required=True,
                               help='Broccoli file with the number of proteins per OG in a genome'

                               )
    requiredArgs.add_argument('-o', '--output_file',
                               dest='out_file',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def parse_table_n_OG(broccoli_table):
    '''
    Reads the table_OGs_protein_counts.txt file from Broccoli and returns
    a table that can be accessed as table[<GENOME>][<OG_ID>] to get the number
    '''
    n_og_table = pd.read_csv(broccoli_table, header=0, sep="\t", index_col=0)
    newcols = [column.split("_prodigal")[0] for column in n_og_table.columns]
    n_og_table.columns = newcols

    return n_og_table



def parse_OGs_faas(OGs_faa_dir, n_og_table):
    '''
    Returns a dictionary with k=protein_id and v=[genome, {NA family}, {NA subfamily}, {NA genus}, length, n_og, OG_id]
    It also takes the input of parse_table_n_OG so the number can be incorporated already in this step.
    '''

    # get faa files
    faas = glob.glob(f"{OGs_faa_dir}/*faa")

    # initiate dictionary
    dictionary = dict()

    # iterate faa files and store to list
    for faa in faas:
        OG_id = os.path.basename(faa).replace(".faa", "")
        records = SeqIO.parse(faa, "fasta")

        # get proteins ids. It is the only information I need for now
        ids = [record.id for record in records]

        for id in ids:
            splits = id.split("|")
            genome = splits[0]
            length = splits[-3]
            #genomic = splits[-1].split("_")[-1] # genomic order
            n_og   = str(n_og_table[genome][OG_id])

            # chimeric proteins appear more than 1 time. Indicate their multiple OGs
            if id in dictionary:
                dictionary[id][3] += f",{OG_id}"
                #print(id, dictionary[id])

            else:
                dictionary[id] = [
                                    genome,
                                    "",
                                    "",
                                    "",
                                    length,
                                    n_og,
                                    OG_id
                                    #genomic
                                 ]

    return dictionary



def parse_add_taxonomy(taxa_file, dictionary):
    '''
    Takes a taxonomy file and adds Family and Subfamily information to the
    dictionary returned by read_OGs_faas
    '''
    # read taxa_file and store in a dictionary, key=genome_id, value=taxa_ranks
    lines = [line.strip().split("\t") for line in open(taxa_file, "r").readlines()[1:]] # discard header line
    # notice the column slices are very specific to this file. It will vary with other taxa file, like Yutin
    taxa = {line[1]:line[17:20] for line in lines}

    for protein, list in dictionary.items():
        if list[0] in taxa:
            #print(list[0])
            dictionary[protein][1] = taxa[list[0]][0]
            dictionary[protein][2] = taxa[list[0]][1]
            dictionary[protein][3] = taxa[list[0]][2]


def yutin_paper_profiles_domains(profiles, domains, dictionary):
    '''
    Ads the profile and/or domain information from the yutin's paper. It uses
    the genomic order of the protein (1_185, 1_47...). This should be OK, all
    the genes I have inspected match with her order and start-end.
    It takes into account if the genome was not in the study
    '''
    pass



def main():

    args = parse_args()

    # Read the n_og information
    n_og_table = parse_table_n_OG(args.n_og_file)

    # Parse the info from the faa files. Also add the n_og information from above
    proteins = parse_OGs_faas(args.OGs_faa, n_og_table)

    # Parse and add the taxonomy information
    parse_add_taxonomy(args.taxa_file, proteins)

    #print(proteins["11_lib200753_1|prodigal-11|80308..82356|683|-|1_95"])

    # transfer dictionary to list, and sort the list before writing
    tow = [[protein] + list for protein, list in proteins.items()]
    tow.sort(key=lambda line: line[6]) # sort per OG

    with open(args.out_file, "w") as fout:
        for line in tow:
            fout.write("\t".join(line) + "\n")




if __name__ == "__main__":
    main()

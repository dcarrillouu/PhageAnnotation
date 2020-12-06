'''
This script reads the new genome tables, the all_final_proteins.faa and the broccoli
names table to create OG_xx.faa files filtered of chimeras and partial proteins.
'''

import argparse
import glob
import os

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

    requiredArgs.add_argument('-gt', '--genome_tables_dir',
                               dest='genome_tables_dir',
                               required=True,
                               help='directory with the new genome tables. It is '
                               '"4_six-frame_approach/2_genome_tables/0_raw"'
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()

def parse_genome_tables(tables_dir, out_dir):
    '''
    Reads genome tables and store in a list chimeric and partial proteins.
    It also returns a file with all the proteins identified as partial or chimeric
    '''

    tables = glob.glob(f"{tables_dir}/*table")

    to_remove = list()
    to_write  = list()
    for table in tables:
        # read the table
        lines = [line.split("\t") for line in open(table).readlines()[1:]]
        # iterate the proteins
        for line in lines:
            if line[5] == "True":
                to_remove.append(line[1])
                to_write.append([line[1], "partial", line[-1].replace("\n", "")])
            if line[6] == "True":
                to_remove.append(line[1])
                to_write.append([line[1], "chimeric", line[-1].replace("\n", "")])

    with open(f"{out_dir}/removed_proteins.txt", "w") as fout:
        for protein in to_write:
            fout.write("\t".join(protein) + "\n")

    return to_remove

def broccoli_to_dict(broccoli_table, to_remove):
    '''
    Converts the broccoli "table_OGs_protein_names.txt" file to dict, k=OG_id, v=[protein_names]
    '''
    # discard header line with the genome information
    lines = [line.strip().split("\t") for line in open(broccoli_table).readlines()[1:]]

    #
    OGs_proteins = dict()
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

        list_to_return = list()
        # check if chimeric/partial and remove if it is
        for protein in flat_list:
            if protein not in to_remove:
                list_to_return.append(protein)


        OGs_proteins[line[0]] = list_to_return

    return OGs_proteins

def dict_to_faa(OG_proteins, faa_file, out_dir):
    '''
    Takes the dictionary of k=OG_id, v=[filtered proteins] and extract these sequences
    from the all_final_proteins.faa file using seqtk
    '''

    for OG, proteins in OG_proteins.items():
        # create an IDs file for seqtk
        with open(f"{out_dir}/{OG}.ids", "w") as fout:
            fout.write("\n".join(proteins))

        # run seqtk
        os.system(f"seqtk subseq {faa_file} {out_dir}/{OG}.ids > {out_dir}/{OG}.faa")

def main():

    args = parse_args()

    # read genome_tables and store which proteins have to be removed
    to_remove = parse_genome_tables(args.genome_tables_dir, args.out_dir)

    # read broccoli protein_names table and store to dict (k=OG_id) while removing chimeric and partial proteins
    #OG_proteins = broccoli_to_dict(args.broccoli_table, to_remove)

    # write sequences of the dictionary to faa file
    #dict_to_faa(OG_proteins, args.faa_file, args.out_dir)












if __name__ == "__main__":
    main()

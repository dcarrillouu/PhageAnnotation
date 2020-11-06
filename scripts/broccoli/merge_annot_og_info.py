'''
Creates a table, one row per protein, with the information about the obtained
annotation and OG formatted as follows:

protein  genome  OG  best_nickname  all_domains
'''

import argparse
import glob
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-pd', '--proteins_dir',
                               dest='proteins_dir',
                               required=True,
                               help='directory with all the proteins that were analyzed.'
                               'Only ".faa" files are read.'
                               )


    requiredArgs.add_argument('-og', '--og_file',
                               dest='og_file',
                               required=True,
                               help='"orthologous_groups.txt" from Broccoli step 3'
                               )
    requiredArgs.add_argument('-a', '--annot_file',
                               dest='annot_file',
                               required=True,
                               help='tabular annotation file'
                               )

    requiredArgs.add_argument('-ad', '--annot-detail_file',
                               dest='annotd_file',
                               required=True,
                               help='detailed annotation file, with all the domains'
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    # Dict to store all the proteins
    proteins = dict()

    # first step is to get all the proteins_ids that were analyzed
    # store the genome as first item of the list
    protein_files = glob.glob(args.proteins_dir + "/*.faa")
    for file in protein_files:
        handle = SeqIO.parse(file, "fasta")
        for record in handle:
            genome = record.id.split("|")[0]
            proteins[record.id] = [genome]


    ### read orthologous_groups.txt and store in dict. k=protein, value=OG_ids
    prot2og = dict()
    lines = [line.strip().split("\t") for line in open(args.og_file).readlines()[1:]]

    # proteins of the OG are in line[1] separated by " "
    for line in lines:
        og_id = line[0]
        proteins_og = line[1].split(" ")
        for protein in proteins_og:
            if protein in prot2og:
                prot2og[protein].append(og_id)
            else:
                prot2og[protein] = [og_id]

    # add OG info to main dict (proteins).
    # iterate main dict. If protein has OG, add it. Else, add empty "".
    for protein, list in proteins.items():
        if protein in prot2og:
            proteins[protein].append(";".join(prot2og[protein]))
        else:
            proteins[protein].append("")


    ### Read best nickname from annotation.txt
    best_nickn = dict()
    lines = [line.strip().split("\t") for line in open(args.annot_file).readlines()] # don't discard header, there isn't

    for line in lines:
        best_nickn[line[0]] = line[1]

    # add best nickname info to main dict (proteins).
    # iterate main dict. If protein has nickname, add it. Else, add empty "".
    for protein, list in proteins.items():
        if protein in best_nickn:
            proteins[protein].append(best_nickn[protein])
        else:
            proteins[protein].append("")


    ### Read all_hits domain from detailed_annotation.txt
    hits = dict()
    lines = [line.strip().split("\t") for line in open(args.annotd_file).readlines()] # don't discard header, there isn't

    for line in lines:
        to_add = f"{line[1]}({line[2]})"
        if line[0] in hits:
            hits[line[0]].append(to_add)
        else:
            hits[line[0]] = [to_add]

    # add best nickname info to main dict (proteins).
    # iterate main dict. If protein has nickname, add it. Else, add empty "".
    for protein, list in proteins.items():
        if protein in hits:
            proteins[protein].append(";".join(hits[protein]))
        else:
            proteins[protein].append("")



    # write results to file
    with open(args.out_dir + "/proteins_annot_og.txt", "w") as fout:
        for protein, list in proteins.items():
            fout.write("{}\t{}\n".format(protein, "\t".join(list)))




















if __name__ == '__main__':
    main()

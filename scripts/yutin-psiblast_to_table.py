'''
Parses the annotation obtained from Yutin's custom profiles (VP_xxxx) and psiblast.

The files are expected to be tabular and formatted as follows:

   0       1      2      3      4       5        6        7     8      9     10   11     12   13    14    15    16     17
qaccver saccver length pident nident positive mismatch gapopen gaps qcovhsp qlen qstart qend slen sstart send evalue bitscore

It also requires the file 'yutin_nicknames.txt' to collapse some annotations into
'HNH', 'DNA_lig' or whatever.

It generates 2 files:
    - annotation.txt, with only one hit per query and no more information
    - annotation_detail.txt: with all the hits information
'''

import argparse
import os
import glob


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_directory',
                               dest='in_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-ext', '--annot_extension',
                               dest='ext',
                               required=True,
                               help='extension of the psiblast files'
                               )

    requiredArgs.add_argument('-n', '--nicknames',
                               dest='in_nickn',
                               required=True,
                               help='tabular "DOM_ID" "NICKNAME"'
                               )

    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def main():

    args = parse_args()

    # read nicknames and store in a dict (notice the upper() )
    nickn = {line.split("\t")[0].upper():line.strip().split("\t")[1] for line in open(args.in_nickn, "r").readlines()}

    # get psiblast files
    files = glob.glob(os.path.abspath(args.in_dir) + "/*" + args.ext)

    # create a dict to store the final hits from all the files
    annot_dict = dict()

    # iterate files. Keep only the last hit for each query in the file
    # Recall that query = yutin_profile, hit = my_proteins
    # one dict per file, so the last file overwrite previous hits
    for file in files:
        file_dict = dict()

        lines = [line.strip().split("\t") for line in open(file, "r").readlines()]

        # remove empty lines or convergence message
        while [""] in lines:
            lines.remove([""])
        while ["Search has CONVERGED!"] in lines:
            lines.remove(["Search has CONVERGED!"])

        # now, iterate lines after cleaning
        for line in lines:
            # recall that queries are like VP00019_consensus, cd00085_0. I need the first part.
            domain = line[0].split("_")[0].upper()

            # key= protein(hit), value=[domain nickname, domain, bitscore, hit coverage]
            # compute hit coverage. Why? Why not.
            hcov = truncate((int(line[2]) - int(line[8])) / int(line[13]), 2)
            file_dict[line[1]] = [nickn[domain], domain, float(line[17]), hcov, line[14], line[15]]

        for protein, annot in file_dict.items():
            if protein in annot_dict:
                annot_dict[protein] += [annot]
                #print(annot_dict[protein])
            else:
                annot_dict[protein] = [annot]


    # Write the detailed file with all the domains
    with open(args.out_dir + "/annotation_detail.txt", "w") as fout:
        for protein, annots in annot_dict.items():
            for annot in annots:
                annot[2], annot[3] = str(annot[2]), str(annot[3])
                fout.write("{}\t{}\n".format(protein, "\t".join(annot)))


    # write the annotation, best domains per protein
    with open(args.out_dir + "/annotation.txt", "w") as fout:
        for protein, annots in annot_dict.items():
                # write the one with the best bitscore
                sorted_hits = sorted(annots, key=lambda x: x[2], reverse=True) # sort by bitscore
                # write the first line only (highest bitscore)
                sorted_hits[0][2], sorted_hits[0][3] = str(sorted_hits[0][2]), str(sorted_hits[0][3])
                fout.write("{}\t{}\n".format(protein, "\t".join(sorted_hits[0])))


if __name__ == "__main__":
    main()

'''
Adds to the shared-content matrix the next items:
  - Taxonomic annotation from Bas
  - Best-coding caller (so I can choose to show it in the heatmap/dendrogram)
  - number of proteins in the genome
'''


import argparse
from pathlib import Path
import os


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--shared-content_table',
                               dest='shared_table_file',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-taxa', '--taxa_file',
                               dest='taxa_file',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help='Bas taxa file "contigs_clustering.txt"'
                               )

    requiredArgs.add_argument('-caller', '--caller_info_file',
                               dest='caller_file',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help='File with the caller used for each genome'
                               )

    requiredArgs.add_argument('-o', '--out_annot-matrix',
                               dest='out_matrix',
                               type=lambda p: Path(p).resolve(),
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def main():
    args = parse_args()

    # read taxa_file and store in a dictionary, key=genome_id, value=taxa_ranks
    lines = [line.strip().split("\t") for line in open(args.taxa_file, "r").readlines()[1:]] # discard header line
    # notice the column slices are very specific to this file. It will vary with other taxa file, like Yutin
    taxa = {line[1]:line[17:20] for line in lines}

    # caller info. I obtain it from the 'n-proteins_best-coding.txt' generated with get_best_coding.py
    caller = {line.split("\t")[0]:line.split("\t")[1] for line in open(args.caller_file).readlines()}
    n_prot = {line.split("\t")[0]:line.strip().split("\t")[2] for line in open(args.caller_file).readlines()}

    with open(args.out_matrix, "w") as fout:
        lines = [line.strip().split("\t") for line in open(args.shared_table_file).readlines()]

        lines[0].insert(0, "") # the strip() removed the first empty item of the header, so I insert it again
        lines[0] += ["Family", "Subfamily", "Genus", "Caller", "n_prots"]
        fout.write("\t".join(lines[0]) + "\n")

        # iterate rows and assign taxonomy if available
        for line in lines[1:]:
            if line[0] in taxa:
                line += taxa[line[0]]
            else:
                line += ["", "", ""]

            line.append(caller[line[0]])
            line.append(n_prot[line[0]])

            fout.write("\t".join(line) + "\n")



if __name__ == '__main__':
    main()

'''
The goal is to assess how supp tRNAs are distributed along the genomes accounting
for their STOP codon reassignments.


'''

import argparse
import glob
import os

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_dir',
                               dest='in_dir',
                               required=True,
                               help='dir with the tRNA prediction files'
                               )

    requiredArgs.add_argument('-ext', '--extension',
                               dest='ext',
                               required=True,
                               help='extension of the tRNAscan files in the input dir '
                               )

    requiredArgs.add_argument('-c', '--coding',
                               dest='coding_file',
                               required=True,
                               help='file describing which was the selected coding for the genome'
                               )

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def parse_tRNAscan_file(in_file, geno_coding, trna_coding):
    genome_id = os.path.basename(in_file).split(".trna")[0]

    lines = [line.split("\t") for line in open(in_file).readlines()[3:]] # remove the header, 3 lines

    add = ""
    if lines: # if the list is not empty, ie if tRNAs were found in that genome
        for line in lines:
            if float(line[-2]) >= 35:
                if line[4] == "Sup":
                    add = line[5]
                    if geno_coding[genome_id] == "prodigal-11":
                        print(genome_id, line[5])

    trna_coding[geno_coding[genome_id]].append(add)

def main():

    args = parse_args()

    # read the coding info and associate it to the genome
    lines = [line.strip().split("\t") for line in open(args.coding_file).readlines()]
    geno_coding = {line[0]:line[1] for line in lines}

    # list input files
    files = glob.glob(args.in_dir + "/*" + args.ext)

    trna_coding = {"prodigal-11":list(),
                   "prodigal-TAG":list(),
                   "prodigal-TGA":list()}

    for file in files:
        parse_tRNAscan_file(file, geno_coding, trna_coding)


    fout = open(args.out_dir + "/supp-tRNA_per_coding.txt", "w")
    for coding, genomes in trna_coding.items():
        total_genomes = len(genomes)
        tow = [coding, str(total_genomes)]

        uniq = list(set(genomes))

        for codon in uniq:
            n = genomes.count(codon)
            tow.append(f"{codon}|{n}")

        fout.write("\t".join(tow) + "\n")

    fout.close()



















if __name__ == "__main__":
    main()

import argparse
from pathlib import Path
import os


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-mcl', '--mcl_outfile',
                               dest='mcl_out',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-gl', '--genomes_lengths',
                               dest='geno_lengths',
                               type=lambda p: Path(p).resolve(strict=True),
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--out_presabs-matrix',
                               dest='out_matrix',
                               type=lambda p: Path(p).resolve(),
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-ogs', '--out_orthologous-groups',
                               dest='out_ogs',
                               type=lambda p: Path(p).resolve(),
                               required=True,
                               help=''
                               )



    return parser.parse_args()


def main():
    args = parse_args()

    # WARNING: if used diamond, be careful with the first line. If you don't remove
    # the "*" lines from the diamond output, they will also appear in the output
    # of mcl, in the first line. So remove it from diamond.
    lines = [line.strip().split("\t") for line in open(args.mcl_out).readlines()]

    # give each cluster an ID and store in a dict
    ogs = dict()
    cont = 1
    for line in lines:
        ogs["OG_" + str(cont)] = line
        cont += 1

    # write to file the OGs. The file is the same as the mcl one, but with the
    # first column as OG_identifier
    with open(args.out_ogs, "w") as fout:
        for og_id, og_proteins in ogs.items():
            fout.write("{}\t{}\n".format(og_id, "\t".join(og_proteins)))



    genomes_id = [line.split("\t")[0] for line in open(args.geno_lengths).readlines()]

    # initialice matrix
    handle_matrix = dict()
    print(cont)
    # one column per each OG: 8,432
    # notice that cont is 8,433 so we will put the genomeid in the first column
    # also notice that this matrix won't be binary: it will contain the number
    # of proteins included in that OG for that genome
    for genome in genomes_id:
        handle_matrix[genome] = [0] * cont

    # go through the orthologs groups
    for og_id, og_proteins in ogs.items():
        n = int(og_id.split("_")[1])

        for protein in og_proteins:
            genome = protein.split("|")[0]
            #print(protein, n, og_id)
            handle_matrix[genome][n] += 1


    # Create header and write to the file
    header = ["OG_{}".format(count) for count in range(1,cont)]
    header.insert(0, "")
    print(len(header))
    with open(args.out_matrix, "w") as fout:
        fout.write("\t".join(header) + "\n")
        for genome, presabs in handle_matrix.items():
            # insert gneomeid in the first (empty) column
            presabs[0] = genome
            tow = [str(item) for item in presabs]

            fout.write("\t".join(tow) + "\n")








if __name__ == '__main__':
    main()

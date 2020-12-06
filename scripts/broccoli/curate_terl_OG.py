from Bio import SeqIO
import argparse

'''

'''


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-faa', '--OG_faa',
                               dest='OG_faa',
                               required=True,
                               help='Terl OG faa file'
                               )
    requiredArgs.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=True,
                               help='deduped out TerL faa'
                               )
    return parser.parse_args()


def main():

    args = parse_args()

    done = list()
    dup_genomes = list()

    records = SeqIO.parse(args.OG_faa, "fasta")

    # obtain the genomes where more than one protein of the OG is found
    for record in records:
        genome = record.id.split("|")[0]
        if genome in done:
            dup_genomes.append(genome)
        else:
            done.append(genome)

    print(dup_genomes)

    # write only those proteins from genomes with only one copy
    records = SeqIO.parse(args.OG_faa, "fasta")
    tow = list()
    for record in records:
        genome = record.id.split("|")[0]
        if genome not in dup_genomes:
            tow.append(record)

    with open(args.out_file, "w") as fout:
        SeqIO.write(tow, fout, "fasta")






if __name__ == "__main__":
    main()

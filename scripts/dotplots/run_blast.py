'''
Runs Blast analysis to check for colinearity.

v1.0:   
        - uses only the 869 genomes from Bas study. Not all of them have
          taxonomic assignment
'''


import argparse
import os
import glob
from functools import partial
import multiprocessing


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_dir',
                               dest='in_dir',
                               required=True,
                               help='Directory with the genomes fasta files'
                               )
    requiredArgs.add_argument('-out', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help='Directory to store the results'
                               )

    requiredArgs.add_argument('-t', '--taxa_file',
                               dest='taxa_file',
                               required=True,
                               help='Tabular file with taxonomic assignments for the genomes'
                               )



    return parser.parse_args()


def makeblastdb(genome, in_dir, out_dir):
    fasta = f"{in_dir}/{genome}.fasta"
    out = f"{out_dir}/1_blastdbs/{genome}"

    os.system(f"/home/danielc/software/ncbi-blast-2.10.1+/bin/makeblastdb -in {fasta} -dbtype nucl -out {out}")



def run_blastn(in_dir, out_dir, genomes):
    # query: genomes[0], reference: genomes[1]

    query = f"{in_dir}/{genomes[0]}.fasta"
    ref   = f"{out_dir}/1_blastdbs/{genomes[1]}"
    out   = f"{out_dir}/2_blast-results/{genomes[0]}'___vs___'{genomes[1]}.blastn"
    #print(out)
    if not os.path.exists(out):
        os.system(f"/home/danielc/software/ncbi-blast-2.10.1+/bin/blastn -query {query} -out {out} -task blastn -db {ref} -outfmt '6 qaccver saccver length pident qlen qstart qend slen sstart send evalue bitscore'")
    else:
        print(f"{out} already exists")


def main():

    args = parse_args()

    # read taxa file and store to dict: k=genome , v=[family, subfamily, genus]
    lines = [line.strip().split("\t") for line in open(args.taxa_file, "r").readlines()[1:]] # discard header line
    # notice the column slices are very specific to this file. It will vary with other taxa file, like Yutin
    taxa = {line[1]:line[17:20] for line in lines if line[17] != "NA"} # discard those with NA assignment
    print(len(taxa))
    # get also the length of the genomes:
    lengths = {line[1]:int(line[2]) for line in lines if line[17] != "NA"}

    # get input fasta files
    files = glob.glob(args.in_dir + "/*.fasta")
    # retain only those with taxa assignment
    for file in files:
        name = os.path.basename(file).split(".fasta")[0]
        if name not in taxa:
            files.remove(file)


    # get subfamilies and store in dict. It will contain largest genome per subfam as value
    subfams = {subfam:"" for subfam in list(set(ranks[1] for genome,ranks in taxa.items()))}


    # get the largest genome per subfamily
    for subfam, largest in subfams.items():
        max = ["genome", 0]
        for genome, ranks in taxa.items():
            if ranks[1] == subfam:
                if lengths[genome] > max[1]:
                    max = [genome, lengths[genome]]
                    subfams[subfam] = [genome, lengths[genome]]

    # create blastdb for the largest genome. Folder blastdbs
    for subfam, largest in subfams.items():
        outdb_path = f"{args.out_dir}/1_blastdbs/{largest[0]}.nhr"
        if not os.path.exists(outdb_path):
            makeblastdb(largest[0], args.in_dir, args.out_dir)
        else:
            print(f"{outdb_path} already exists")


    # prepare blast analysis for all the genomes
    to_run = list()
    for subfam, largest in subfams.items():
        for genome, ranks in taxa.items():
            if ranks[1] == subfam and genome != largest[0]: # avoid self aligning genomes
                to_run.append([genome, largest[0]]) # genome vs largest_subfamily_genome


    blast = partial(run_blastn, args.in_dir, args.out_dir )
    pool = multiprocessing.Pool(processes=8)
    pool.map(blast, to_run)
    pool.close()
    pool.join()








if __name__ == "__main__":
    main()

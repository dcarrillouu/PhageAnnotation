'''
This script takes the tables generated with "yutin_to_genome_tables.py" and the
taxonomic information to generate a unique file containing all the profiles and
taxonomies found within every OG. In this way I can see if a OG is specific for
a function/taxonomy.

I also grab OG information from "statistics_per_OG.txt"

'''

import argparse
import glob
import os

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--tables_directory',
                               dest='tables_dir',
                               required=True,
                               help='directory with genome tables annotated with OG and yutin.'
                               'Generated with "yutin_to_genome-tables.py", extension ".genes-OG-profiles"'
                               )
    requiredArgs.add_argument('-o', '--output_file',
                               dest='out_file',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-taxa', '--taxa_file',
                               dest='taxa_file',
                               required=True,
                               help='file with the taxonomic information'
                               )
    requiredArgs.add_argument('-stats', '--ogs_stats_file',
                               dest='ogs_stats_file',
                               required=True,
                               help='file "statistics_per_OG.txt"'
                               )


    return parser.parse_args()



def main():

    args = parse_args()

    # read taxa_file and store in a dictionary, key=genome_id, value=taxa_rank
    lines = [line.strip().split("\t") for line in open(args.taxa_file, "r").readlines()[1:]] # discard header line
    # notice the column slices are very specific to this file. It will vary with other taxa file, like Yutin
    families = {line[1]:line[17] for line in lines}
    subfamilies = {line[1]:line[18] for line in lines}
    genera = {line[1]:line[19] for line in lines}


    # read OG statistics
    OGs_stats = {line.split("\t")[0]:line.strip().split("\t")[1:] for line in open(args.ogs_stats_file).readlines()[1:]}
    OGs_sorted = [line[0] for line in open(args.ogs_stats_file).readlines()[1:]]


    # create two dictionaries to store OG information: one for profiles, another for taxa
    # start them with the OGs (keys) of OGs_stats
    OGs_profiles  = {og:list() for og, stats in OGs_stats.items()}
    OGs_family    = {og:list() for og, stats in OGs_stats.items()}
    OGs_subfamily = {og:list() for og, stats in OGs_stats.items()}
    OGs_genus     = {og:list() for og, stats in OGs_stats.items()}

    for table in glob.glob(f"{args.tables_dir}/*.genes-OG-profiles"):
        lines = [line.split("\t") for line in open(table).readlines()[1:]]
        # itereate through the proteins
        for line in lines:
            # check if the protein has OG
            if "OG_" in line[-1]:
                og_id = line[-1].strip() # reomove the "\n"

                # add profile information
                if line[-2] != "":
                    OGs_profiles[og_id] += line[-2].split(",")

                # add taxa information. Check if genome is annotated with taxa
                if line[0] in families:
                    OGs_family[og_id].append(families[line[0]])
                    OGs_subfamily[og_id].append(subfamilies[line[0]])
                    OGs_genus[og_id].append(genera[line[0]])




    for OG, stats in OGs_stats.items():
        # format profiles info
        tow = list()
        if OGs_profiles[OG]: # if list is not empty
            nicknames = list(set(OGs_profiles[OG]))
            for nickname in nicknames:
                n = OGs_profiles[OG].count(nickname)
                tow.append(f"{nickname} ({n})")
        if tow:
            tow = ", ".join(tow)
            OGs_stats[OG].append(tow)
        else:
            OGs_stats[OG].append("")

        # format family info
        tow = list()
        if OGs_family[OG]:
            taxas = list(set(OGs_family[OG]))
            for taxa in taxas:
                n = OGs_family[OG].count(taxa)
                tow.append(f"{taxa} ({n})")
        if tow:
            tow = ", ".join(tow)
            OGs_stats[OG].append(tow)
        else:
            OGs_stats[OG].append("")

        # format subfamily info
        tow = list()
        if OGs_subfamily[OG]:
            taxas = list(set(OGs_subfamily[OG]))
            for taxa in taxas:
                n = OGs_subfamily[OG].count(taxa)
                tow.append(f"{taxa} ({n})")
        if tow:
            tow = ", ".join(tow)
            OGs_stats[OG].append(tow)
        else:
            OGs_stats[OG].append("")

        # format subfamily info
        tow = list()
        if OGs_genus[OG]:
            taxas = list(set(OGs_genus[OG]))
            for taxa in taxas:
                n = OGs_genus[OG].count(taxa)
                tow.append(f"{taxa} ({n})")
        if tow:
            tow = ", ".join(tow)
            OGs_stats[OG].append(tow)
        else:
            OGs_stats[OG].append("")


    with open(args.out_file, "w") as fout:

        for OG, stats in OGs_stats.items():
            fout.write(f"{OG}\t" + "\t".join(stats) + "\n")





if __name__ == '__main__':
    main()

'''
Create genome tables with yutin profiles and OG information.
It takes the genome_tables obtained with "broccoli_to_gggenes-tables.py" and
the "annotation_detail.txt" from "yutin-psiblast_to_table.py".

It returns all the hits for the protein indicating start and stop positions in
the protein. If the same domain has several hits of the same profile and they
overlap, I keep one with the best bitscore. For instance, gene49 here, I would
take the first hit:

OGOH01000022|prodigal-TAG|79851..81692|614|-|1_78       gene49  VP02751 467.0   0.98    1       607
OGOH01000022|prodigal-TAG|79851..81692|614|-|1_78       gene49  VP00394 210.0   0.93    2       613


UPDATE: I end up not collapsing by start-stop position. Is not always straightfoward with
profiles like the RNApBprime. Given that profiles are redundant, I obtain hits from
all of them and can not distinguish wich is good and don't overlap with other realiable
profiles like RNApB. Look at OLDR01000012|prodigal-TAG|16535..29857|4441|-|1_17.
So I annotate the tables with only the non-redundant nicknames.
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
                               help='directory with genome tables annotated with OG.'
                               'Generated with "broccoli_to_gggenes-tables.py", extension ".genes-OG"'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-annot', '--annot_detail_file',
                               dest='annot_detail_file',
                               required=True,
                               help='file with the yutins domain/profile psiblast annot'
                               )


    return parser.parse_args()



def main():

    args = parse_args()

    # read annotation file and store to dict. k=protein, v=[[annot1], [annot2]]
    # each annot contains profile_nickname, profile_ID, bitscore, hit_cov, start, stop
    annot_file = [line.strip().split("\t") for line in open(args.annot_detail_file).readlines()]

    annot = dict()
    for line in annot_file:
        if line[0] in annot:
            annot[line[0]].append(line[1])
        else:
            annot[line[0]] = [line[1]]

    for protein, nicknames in annot.items():
        annot[protein] = ",".join(list(set(nicknames)))


    #print(annot["OLDR01000012|prodigal-TAG|16535..29857|4441|-|1_17"])

    # read genomes tables and annotate them
    tables = glob.glob(f"{args.tables_dir}/*.genes-OG")
    for table in tables:
        outfile = f"{args.out_dir}/{os.path.basename(table)}-profiles"
        lines = [line.split("\t") for line in open(table).readlines()]

        # insert the profile annotation in the column previous to OG info
        lines[0].insert(-1, "yutin_profiles")

        for line in lines[1:]:
            # protein_id in line[4]
            if line[4] in annot: # if there is any profile annotated for this protein
                line.insert(-1, annot[line[4]])
            else: # else write empty column
                line.insert(-1, "")

        # write to fileº
        with open(outfile, "w") as fout:
            for line in lines:
                fout.write("\t".join(line))




















if __name__ == '__main__':
    main()

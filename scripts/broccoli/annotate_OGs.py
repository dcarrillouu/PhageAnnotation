'''
It takes the 'orthologous_groups.txt' from Broccoli, step3, and assigns the functional
annotation of all its proteins.
It returns two files:
    - broccoli_ogs_nicknames.txt, with the nicknames associated with the og proteins
    - broccoli_ogs_domains.txt, with the domains associated with the og proteins

Annotation file needs to be a tabular file with 'protein' \t 'annot'.
'''

import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

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

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    # Read annot file and store in a dict. Store nickname and profile_id
    annot = {line.split("\t")[0]:line.split("\t")[1:3] for line in open(args.annot_file).readlines()}

    # Read OGs file, remove header and save to list
    # [0] = OG_id, [1] = proteins in the OG
    ogs = [line.strip().split("\t") for line in open(args.og_file).readlines()[1:]]

    # init dicts for og_nicknames and og_domains
    og_nicknm = dict()
    og_domain = dict()

    # init matrix to store lines that will be writen to file
    matrix_nickn  = list()
    matrix_domain = list()

    # iterate ogs and count the number of annots
    for og in ogs:
        og_id = og[0]
        og_nicknm[og_id] = list()
        og_domain[og_id] = list()

        # put proteins on a list to iterate them
        proteins = og[1].split(" ") # proteins separated by " "
        og_length = len(proteins)

        for protein in proteins:
            if protein in annot:
                og_nicknm[og_id].append(annot[protein][0])
                og_domain[og_id].append(f"{annot[protein][1]}|{annot[protein][0]}") # notice I also store the nickname along the domain

        nicknames = list(set(og_nicknm[og_id]))
        domains   = list(set(og_domain[og_id]))

        # process nickname and store in a sorted way
        counts = dict()
        for nickname in nicknames:
            counts[nickname] = og_nicknm[og_id].count(nickname)

        tow = [og_id, og_length]
        for nickn in sorted(counts, key=counts.get, reverse=True):
            #print(og_id, nickn, counts[nickn])
            tow.append(f"{nickn}|{counts[nickn]}")
        matrix_nickn.append(tow)


        # process domain and store in a sorted way
        counts = dict()
        for domain in domains:
            counts[domain] = og_domain[og_id].count(domain)

        tow = [og_id, og_length]
        for domain in sorted(counts, key=counts.get, reverse=True):
            #print(og_id, nickn, counts[nickn])
            tow.append(f"{domain}|{counts[domain]}")
        matrix_domain.append(tow)

        #print(matrix_nickn)
        # sort matrices by most populated og and write to files
        sort_matrix_nickn  = sorted(matrix_nickn, key=lambda x: int(x[1]), reverse=True)
        sort_matrix_domain = sorted(matrix_domain, key=lambda x: int(x[1]), reverse=True)

        with open(args.out_dir + "/broccoli_ogs_nicknames.txt", "w") as fout:
            for line in sort_matrix_nickn:
                line[1] = str(line[1]) # number of proteins in og, string to write to file
                fout.write("\t".join(line) + "\n")


        with open(args.out_dir + "/broccoli_ogs_domains.txt", "w") as fout:
            for line in sort_matrix_domain:
                line[1] = str(line[1]) # number of proteins in og, string to write to file
                fout.write("\t".join(line) + "\n")












if __name__ == "__main__":
    main()

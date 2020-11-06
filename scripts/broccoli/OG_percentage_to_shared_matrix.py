'''
This script takes a shared-content matrix and adds the annotation column 'perc_OG'
meaning the percentage of the total proteins that were assigned to an OG. For that,
it requires the "statistics_per_species.txt" from Broccoli step3
'''

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_shared-content-matrix',
                               dest='in_matrix',
                               required=True,
                               help=''
                               )


    requiredArgs.add_argument('-og', '--broccoli_perc-ogs_file',
                               dest='ogs_file',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    # read broccoli file and assign to dict
    lines = [line.strip().split("\t") for line in open(args.ogs_file).readlines()[1:]]
    percs_ogs = {line[0].split("_prodigal")[0]:float(line[1]) for line in lines}

    # some of the percs are > 100 (don't know why yet). Change them to 100.
    for genome, perc in percs_ogs.items():
        if perc > 100:
            percs_ogs[genome] = "100"
        else:
            percs_ogs[genome] = str(perc)


    # read shared-content matrix
    lines = [line.strip().split("\t") for line in open(args.in_matrix).readlines()]

    # Create header
    lines[0].insert(0, "")
    lines[0].append("perc_OG")


    # iterate lines and add the perc associtated to the genome
    with open(args.in_matrix.replace(".txt", "_OGperc.txt"), "w") as fout:
        fout.write("\t".join(lines[0]) + "\n")
        for line in lines[1:]:
            line.append(percs_ogs[line[0]])
            fout.write("\t".join(line) + "\n")
















if __name__ == "__main__":
    main()

'''

'''


import argparse
import os
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from randomcolor import RandomColor



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-og', '--og_id',
                               dest='og_id',
                               required=True,
                               help='"orthologous_groups.txt" from Broccoli step 3'
                               )
    requiredArgs.add_argument('-t', '--genome-tables_dir',
                               dest='gtables_dir',
                               required=True,
                               help='tabular annotation file'
                               )

    requiredArgs.add_argument('-n', '--n_context',
                               dest='n_context',
                               required=True,
                               help='number of context genes to show'
                               )


    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def grep_to_dictionary(grep_file):
    '''
    1) Removes empty "--" lines produced by "end of file" in grep.
    2) Formats the start-stop positions of the genes so the region starts at 0
    3) Stores each genome as key in a dictionary
    '''
    # read the file and remove duplicates. Can not split yet to do this.
    lines_raw   = [line for line in open(grep_file).readlines()[1:]] # discard header
    lines_dedup = list(set(lines_raw))

    # split the lines
    lines = [line.split("\t") for line in lines_dedup]

    # remove end_of_file lines
    while ["--\n"] in lines:
        lines.remove(["--\n"])

    # remove "\n" from the OG column. convert start and stop to int. Change strand values to +1 or -1
    for line in lines:
        line[-1] = line[-1].replace("\n", "")
        line[1]  = int(line[1])
        line[2]  = int(line[2])
        if line[3] == "+":
            line[3] = 1
        else:
            line[3] = -1

    # Store genomes in a dict
    genomes = {line[0]:list() for line in lines}
    for line in lines:
        genomes[line[0]].append(line)

    # sort genes in a genome and make the first start at 20
    for genome, genes in genomes.items():
        sorted_genes = sorted(genes, key=lambda gene: gene[1])
        # the start position of the first sets up the displacement for the rest of genes
        displacement = sorted_genes[0][1] - 20

        for gene in sorted_genes:
            gene[1] -= displacement
            gene[2] -= displacement

        genomes[genome] = sorted_genes

    return genomes


def choose_angle(strand):
    angle = 45
    if strand == -1:
        angle = 160
    return angle

def create_colors(genomes_dict, og_id):
    '''

    '''
    OGs = list()
    for genome, genes in genomes_dict.items():
        for gene in genes:
            OGs.append(gene[-1])
    print(OGs)

    OGs = list(set(OGs))

    colors = {OG:RandomColor().generate()[0] for OG in OGs}

    # set colors for non_OG genes and target OG_ID
    colors[og_id] = "red"
    if "" in colors:
        colors[""] = "white"

    return colors

def draw_diagram(genomes_dict, og_id, n_context, colors):
    '''

    '''

    gdd = GenomeDiagram.Diagram('Test Diagram')

    i = 0
    for genome, genes in genomes_dict.items():
        track = gdd.new_track(i, greytrack=True, name=genome)
        features = track.new_set()

        for gene in genes:
            feature = SeqFeature(FeatureLocation(gene[1], gene[2]), strand=gene[3])
            features.add_feature(feature, name=f"{gene[-1]}   {gene[-2]}", label=True, color=colors[gene[-1]], label_angle=choose_angle(gene[3]), sigil="ARROW",
                                    arrowshaft_height=1, label_size=10, label_position="middle",)
        i += 1

    gdd.draw(format='linear', fragments=1)
    gdd.write(f"{og_id}_{n_context}.png", "png")





def main():

    args = parse_args()

    # put the OG and its context into a file
    os.system(f"grep -h -C {args.n_context} '{args.og_id}$' {args.gtables_dir}/*genes-OG-profiles | grep -v 'complete_id' > {args.out_dir}/{args.og_id}_{args.n_context}.grep")

    # read grep file and store genomes in a gene' sorted dict. Displace start-stop positions
    genomes = grep_to_dictionary(f"{args.out_dir}/{args.og_id}_{args.n_context}.grep")

    # create colors
    colors = create_colors(genomes, args.og_id)

    #
    draw_diagram(genomes, args.og_id, args.n_context, colors)


if __name__ == "__main__":
    main()

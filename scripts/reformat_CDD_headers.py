'''
Reformat the headers of MSAs downloaded from CDD because some of them make
psiblast fail. Recall that I don't introduce consensus sequence because it
decreases the performance of psiblast in this case.
'''

import argparse
import os
import glob


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_dir',
                               dest='in_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-ext', '--files_extension',
                               dest='ext',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()





def main():

    args = parse_args()

    # get input MSA files
    print(args.in_dir)
    files = glob.glob(os.path.abspath(args.in_dir) + "/*" + args.ext)

    for file in files:
        out = open(os.path.abspath(args.out_dir) + "/" + os.path.basename(file).replace("." + args.ext, "_rehead.fasta"), "w")

        # split the file by its sequences
        seqs = [seq.split("\n") for seq in open(file, "r").read().split(">")]
        seqs.remove([""])
        id = os.path.basename(file).replace("." + args.ext, "")

        cont = 0
        for seq in seqs:
            #print(seq)
            seq[0] = ">{}_{}".format(id, cont) # reformat header
            out.write("\n".join(seq))
            cont += 1
        out.close()















if __name__ == '__main__':
    main()

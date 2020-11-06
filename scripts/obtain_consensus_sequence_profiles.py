'''
Computes the consensus sequence from a multiple sequence alignment, and puts it
in the first line of a FASTA.

I use it to obtain the consensus of the MSA downloaded from CDD

Dependencies:
    - emboss (Installed in "emboss" env in my case)
'''

import argparse
import os
import glob
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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

    # obtain consensus
    for file in files:
        cons_file = file.replace("." + args.ext, ".consensus")
        id = os.path.basename(file).replace("." + args.ext, "")
        tmp_file = os.path.abspath(args.out_dir) + "/" + id + ".afac"
        out_file = os.path.abspath(args.out_dir) + "/" + id + ".afac"

        # obtain consensus with emboss cons, read the file and write to tmp
        os.system("cons -sequence {} -outseq {} -name {}_consensus".format(file, cons_file, id))
        cons = open(cons_file).read()
        tmp = open(tmp_file, "w")
        tmp.write(cons)

        # change headers of CDD MSA, some give error when blast
        seqs = [seq.split("\n") for seq in open(file, "r").read().split(">")]
        seqs.remove([""])

        cont = 0
        for seq in seqs:
            #print(seq)
            seq[0] = ">Seq_{}".format(cont) # reformat header
            tmp.write("\n".join(seq))
            cont += 1

        tmp.close()

        # parse tmp file with AlignIO and write final .afac file
        handle = AlignIO.read(tmp_file, "fasta")
        AlignIO.write(handle, out_file, "fasta")












if __name__ == '__main__':
    main()

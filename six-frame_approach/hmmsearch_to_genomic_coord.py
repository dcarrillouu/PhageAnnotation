'''

'''

import argparse
import glob
import os
import multiprocessing
from functools import partial
from Bio import SeqIO, SearchIO
import statistics



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hmmsearch_file',
                               dest='hmmsearch_file',
                               required=True,
                               help='plain text file with the hmmsearch results'
                               )
    requiredArgs.add_argument('-l', '--lengths_file',
                               dest='lengths_file',
                               required=True,
                               help='file with the genome lengths "coding_statistics.txt"'
                               )

    requiredArgs.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=True,
                               help='file to write the genomic coordinates of the hits'
                               )

    return parser.parse_args()


def translate_to_genomic_coords(start, end, frame, genome_size):
    """
    Translate the coordinates of the protein from transeq to genomic
    coordinates.

    Strand is used here as orientation and not really [non-]conding.
    If the frame is 1,2 or 3 (-->) I call this (+) strand.
    Else, it is the (-) strand

    param: int: start The starting coordinate on the protein
    param: int: end The ending coordinate on the protein
    param: int: frame The frame on which it is found (1-6)
    param: int: genome_size The size of the genomes

    return: tuple: (genomic start, genomic end, strand)
    """
    nucleic_start = start * 3
    nucleic_end = end * 3
    if frame == 1:
        genomic_start = nucleic_start - 2
        genomic_end = nucleic_end - 2
    if frame == 2:
        genomic_start = nucleic_start - 1
        genomic_end = nucleic_end - 1
    if frame == 3:
        genomic_start = nucleic_start
        genomic_end = nucleic_end
    if frame == 4:
        genomic_start = genome_size - (nucleic_start - 2)
        genomic_end = genome_size -  (nucleic_end - 2)
    if frame == 5:
        genomic_start = genome_size - (nucleic_start - 1)
        genomic_end = genome_size - (nucleic_end -1)
    if frame == 6:
        genomic_start = genome_size - nucleic_start
        genomic_end = genome_size - nucleic_end

    if frame in [1,2,3]:
        strand = '+'
    elif frame in [4,5,6]:
        strand = '-'
    else:
        raise ValueError("frame should be one of 1,2,3,4,5,6")

    return genomic_start, genomic_end, strand


def parse_hmmsearch_record(hmmsearch_file, lenghts):
    '''
    Each record contains the hits of the genomes to one OG.
    '''

    to_return = dict()

    records = SearchIO.parse(hmmsearch_file, "hmmer3-text")

    for record in records:
        # get the OG for which genome hits were detected
        OG_id = record.id
        to_return[OG_id] = list()
        # iterate the genome hits
        for hit in record.hits:
            if hit.is_included:
                # get genome id and its frame
                frame = int(hit.id.split("_")[-1])
                genome_id = "_".join(hit.id.split("_")[:-1])
                # iterate the different hsps of the hit
                for hsp in hit.hsps:
                    if hsp.is_included:
                        genomic_coords = translate_to_genomic_coords(hsp.env_start,
                                                                     hsp.env_end,
                                                                     frame,
                                                                     lenghts[genome_id])

                        to_return[OG_id].append([genome_id,
                                                 str(frame),
                                                 str(hsp.evalue),
                                                 str(hsp.bitscore),
                                                 str(hsp.query_start),
                                                 str(hsp.query_end),
                                                 str(min(genomic_coords[0], genomic_coords[1])),
                                                 str(max(genomic_coords[0], genomic_coords[1])),
                                                 str(genomic_coords[2])
                                                 ])


    return to_return


def main():

    args = parse_args()

    # read lengths file and store to dict
    lengths = {line.split("\t")[0]:int(line.split("\t")[1]) for line in open(args.lengths_file).readlines()}


    # parshe hmmsearch file and store records to list
    OGs_hits = parse_hmmsearch_record(args.hmmsearch_file, lengths)

    with open(args.out_file, "w") as fout:
        for OG, hits in OGs_hits.items():
            # sort the hits by genome
            genome_sorted = sorted(hits, key=lambda x: x[0])
            for hit in genome_sorted:
                fout.write(f"{OG}\t" + "\t".join(hit) + "\n")



if __name__ == '__main__':
    main()

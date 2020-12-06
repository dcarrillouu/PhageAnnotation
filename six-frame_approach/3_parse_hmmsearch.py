'''
Script to parse results from hmmsearch. If -t is provided, results of a transeq
scan are supposed as the input.

'''

import argparse
import glob
import os
from Bio import SearchIO




def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--hmmsearch_file',
                               dest='hmmsearch_file',
                               required=True,
                               help='plain text file with the hmmsearch results. '
                               'It can be the results of scanning proteins or transeq genomes'
                               )

    parser.add_argument('-l', '--lengths_file',
                               dest='lengths_file',
                               required=False,
                               help='file with the genome lengths "coding_statistics.txt"'
                               )

    parser.add_argument('-transeq', '--transeq',
                               dest='transeq_check',
                               required=False,
                               default=False,
                               const=True,
                               nargs='?',
                               help='marks the input as he result of a transeq search'
                               )

    parser.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=False,
                               help=''
                               )

    return parser.parse_args()

def translate_to_genomic_coords(start, end, frame, genome_size):
    """
    Translate the coordinates of the protein from transeq to genomic
    coordinates.

    Strand is used here as orientation and not really [non-]conding.
    If the frame is 1,2 or 3 (-->) I call this (+) strand.
    Else, it is the (-) strand

    param: int: start The starting coordinate on the transeq protein
    param: int: end The ending coordinate on the transeq protein
    param: int: frame The frame on which it is found (1-6)
    param: int: genome_size The size of the genomes
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

    return genomic_start, genomic_end

def parse_transeq_hmmsearch(hmmsearch_file, lenghts):
    '''
    Each record contains the hits of the genomes to one OG.
    '''

    to_return = dict()

    records = SearchIO.parse(hmmsearch_file, "hmmer3-text")

    for record in records:
        # get the OG for which genome hits were detected
        OG_id = record.id
        OG_len = record.seq_len
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
                                                 str(OG_len),
                                                 str(hsp.query_start),
                                                 str(hsp.query_end),
                                                 str(min(genomic_coords[0], genomic_coords[1])),
                                                 str(max(genomic_coords[0], genomic_coords[1]))
                                                 ])


    return to_return

def main():

    args = parse_args()

    # checker wether the results are from transeq or from proteins
    if args.transeq_check:

            # read lengths file and store to dict
            lengths = {line.split("\t")[0]:int(line.split("\t")[1]) for line in open(args.lengths_file).readlines()}

            # parse hmmsearch file and store records to list
            OGs_hits = parse_transeq_hmmsearch(args.hmmsearch_file, lengths)

            with open(args.out_file, "w") as fout:
                for OG, hits in OGs_hits.items():
                    # sort the hits by genome
                    genome_sorted = sorted(hits, key=lambda x: x[0])
                    for hit in genome_sorted:
                        fout.write(f"{OG}\t" + "\t".join(hit) + "\n")

    else:
        pass





















if __name__ == "__main__":
    main()

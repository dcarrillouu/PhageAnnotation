'''
Separates the OGs .hmm in different folders to run hmmsearch in different mutants.
It first separates the OGs according to the number of sequences, for which it looks
at the .summary file generated with hmmbuild
'''

import glob
import os

not_align = []

folders = ["hmmsearch1", "hmmsearch2", "hmmsearch3"] * 5000

summaries = glob.glob("*.summary")
n_seqs = list()
for file in summaries:
    # reading the file and discarding the "#" starting sequences returns this kind of line:
    # 1     OG_1001                 39   338   275     0.80  0.589
    lines = open(file).readlines()
    for line in lines:
        if not line.startswith("#") and not line.startswith("\n"):
                split = line.split(" ")
                while "" in split:
                    split.remove("")

                # split is now like this:
                # ['1', 'OG_485', '30', '108', '90', '0.99', '0.635', '\n']
                n_seqs.append([split[1], int(split[2])])

order = sorted(n_seqs, key=lambda OG: int(OG[1]), reverse=True)
print(order)
# #print(order)
#
i = 0
for og in order:
    os.system(f"cp {og[0]}.hmm {folders[i]}")
    i += 1

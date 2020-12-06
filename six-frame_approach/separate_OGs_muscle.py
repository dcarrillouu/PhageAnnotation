'''
Separates the OGs.faa in different folders to run MUSCLE in different mutant
simultaneously.

There are some OGs that I don't want to align yet. They are in the list "not_align"
'''

import glob
import os

not_align = []

folders = ["muscle1", "muscle2", "muscle3"] * 5000

os.system("grep -c '>' *faa > n-seqs.txt")
lines = [line.strip().split(".faa:") for line in open("n-seqs.txt").readlines()]

order = sorted(lines, key=lambda OG: int(OG[1]), reverse=True)

#print(order)

i = 0
for og in order:
    if og[0].split("OG_")[1] not in not_align:
        os.system(f"cp {og[0]}.faa {folders[i]}")
        i += 1
    else:
        print(og)

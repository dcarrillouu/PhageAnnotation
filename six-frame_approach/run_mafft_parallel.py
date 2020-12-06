import glob
import os
import multiprocessing
import time

# count the number of sequences in each OG_xx.faa and sort them from max to min
os.system("grep -c '>' *faa > n-seqs.txt")
lines = [line.strip().split(".faa:") for line in open("n-seqs.txt").readlines()]
order = sorted(lines, key=lambda OG: int(OG[1]), reverse=True)

# define out dir for the alignments
out_dir = "~/projects/Bas_phages/3_Broccoli/run1/4_six-frame_approach/1_filtered_proteins/1_mafft/2_easy_OGs/"

files = glob.glob("*faa")

for OG in order:
    start = time.time()
    print(OG)
    os.system(f"mafft-einsi --thread 40 {OG[0]}.faa > {out_dir}/{OG[0]}.mafft-einsi 2> {out_dir}/{OG[0]}.mafft-einsi.log")
    end = time.time()
    print("\t",end-start)

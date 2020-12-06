import glob
import os
import multiprocessing
import time

# count the number of sequences in each OG_xx.faa and sort them from max to min
os.system("grep -c '>' *.mafft-einsi > n-seqs.txt")
lines = [line.strip().split(".mafft-einsi:") for line in open("n-seqs.txt").readlines()]
order = sorted(lines, key=lambda OG: int(OG[1]), reverse=True)

# define out dir for the alignments
out_dir = "/home/danielc/projects/Bas_phages/3_Broccoli/run1/4_six-frame_approach/1_filtered_proteins/5_nr_scan/0_easy_OGs/psiblast/"

files = glob.glob("*mafft-einsi")

for OG in order:
    print(OG)
    os.system(f"time /home/danielc/software/ncbi-blast-2.10.1+/bin/psiblast -in_msa {OG[0]}.mafft-einsi -db /net/phage/linuxhome/tina/BLAST/NR/nr -out {out_dir}/{OG[0]}.nr-psiblast -evalue 1e-4 -outfmt '6 qaccver saccver pident length mismatch gapopen gaps qcovhsp qlen qstart qend slen sstart send evalue bitscore' -num_iterations 5 -num_threads 70")
    print("\n")

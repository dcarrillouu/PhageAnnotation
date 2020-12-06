'''
Separates the OGs .hmm in different folders to run hmmsearch in different mutants.
It first separates the OGs according to the number of sequences, for which it looks
at the .summary file generated with hmmbuild
'''

import glob
import os
import multiprocessing


summaries = glob.glob("../*.summary")
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
                n_seqs.append([f"{split[1]}.hmm", int(split[2])])

order = sorted(n_seqs, key=lambda OG: int(OG[1]), reverse=True)
print(order[:10])

##
hmms = glob.glob("*.hmm")

hmms_to_run = list()

for og in order:
    if og[0] in hmms and og[1] > 20:
        hmms_to_run.append(og[0])

print(hmms_to_run[:10])

def run_hmm_hmmsearch(hmm_file):
    og_id = os.path.basename(hmm_file).split(".hmm")[0]
    print(og_id)
    os.system(f"time hmmsearch --cpu 5 -o /home/danielc/projects/Bas_phages/3_Broccoli/run1/4_six-frame_approach/1_filtered_proteins/5_nr_scan/0_easy_OGs/hmmsearch/{og_id}.hmmsearch --domtblout /home/danielc/projects/Bas_phages/3_Broccoli/run1/4_six-frame_approach/1_filtered_proteins/5_nr_scan/0_easy_OGs/hmmsearch/{og_id}.domtblout {hmm_file} /net/phage/linuxhome/tina/BLAST/NR/fasta/nr")


print(len(hmms_to_run))

pool = multiprocessing.Pool(processes=12)
shared_content = pool.map(run_hmm_hmmsearch, hmms_to_run)
pool.close()
pool.join()

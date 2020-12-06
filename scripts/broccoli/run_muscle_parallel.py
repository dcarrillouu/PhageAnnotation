import glob
import os
import multiprocessing


files = glob.glob("*faa")


def run_muscle(infile):
    outfile = "/home/danielc/projects/Bas_phages/3_Broccoli/run1/4_six-frame_approach/1_filtered_proteins/1_muscle/0_dupl_OGs/" + os.path.basename(infile).replace(".faa", ".msa")
    log = "/home/danielc/projects/Bas_phages/3_Broccoli/run1/4_six-frame_approach/1_filtered_proteins/1_muscle/0_dupl_OGs/" + os.path.basename(infile).replace(".faa", ".log")
    os.system(f"muscle -in {infile} -out {outfile} -log {log}")



pool = multiprocessing.Pool(processes=4)
pool.map(run_muscle, files)
pool.close()
pool.join()

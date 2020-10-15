'''
Isolates the proteomes of the 869 clustered genomes in a folder, so I can
run OrthoFinder and Broccoli only on those genomes
'''

genomes = [line.split("\t")[1] for line in open("/home/danielc/projects/Bas_phages/0_raw/Bas/contig_clustering.csv").readlines()[1:]]

import glob
import os
faas = glob.glob("/home/danielc/projects/Bas_phages/ORF_calling/final_annotation/*faa")

tow = list()
cont = 0
for genome in genomes:
    for faa in faas:
        if os.path.basename(faa).startswith(genome+ "_prodigal"):
                os.system("cp {} /home/danielc/projects/Bas_phages/ORF_calling/final_annotation/proteomes_869_Bas/".format(faa))
                cont += 1

print(cont)

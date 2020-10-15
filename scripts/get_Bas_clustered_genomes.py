genomes = [line.split("\t")[1] for line in open("/home/danielc/projects/Bas_phages/0_raw/Bas/contig_clustering.csv").readlines()[1:]]

import glob
import os
faas = glob.glob("/home/danielc/projects/Bas_phages/ORF_calling/final_annotation/*faa")

tow = list()
cont = 0
for genome in genomes:
    for faa in faas:
        if os.path.basename(faa).startswith(genome+ "_prodigal"):
                print(faa)
                tow.append(faa)
                cont += 1

print(cont)
os.system("cat {} > /home/danielc/projects/Bas_phages/ORF_calling/final_annotation/diamond_db/869_bas_clustered_genomes.faa".format(" ".join(tow)))

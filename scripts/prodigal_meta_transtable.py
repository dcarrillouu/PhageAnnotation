import glob
import os

gffs = glob.glob("/home/danielc/projects/Bas_phages/ORF_calling/test_prodigal-meta/*..gff")

out = open("/home/danielc/projects/Bas_phages/ORF_calling/trcode_test-prodigal-meta.txt", "w")
for gff in gffs:
    code = open(gff).readlines()[2].split("transl_table=")[1].split(";")[0]
    genome = os.path.basename(gff).replace("..gff", "")
    out.write(f'{genome}\t{code}\n')

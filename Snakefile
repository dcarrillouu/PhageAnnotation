import pathlib
import glob
import os
from Bio import SeqIO
import sys

'''
DEPENDENCIES:
    - trnascan2
    - infernal
    - checkv

TODO:
    - trnascan2.0 dependency (which version is in conda?)
'''


input_genomes = glob.glob("../genomes_fasta/*fasta")
genomes_id    = [os.path.basename(genome).replace(".fasta", "") for genome in input_genomes]
input_path    = os.path.dirname(input_genomes[0])
scripts_path  = os.path.abspath("/home/dani/programs/PhageAnnotation_devel/PhageAnnotation/scripts")
sys.path.append(scripts_path)
import process_phanotate
import process_prodigal
import tables
import pVOG

prodigal_path  = "/home/dani/programs/software.prodigal/tmp.prodigal/code.prodigal/prodigal"
phanotate_path = "/home/dani/programs/PHANOTATE/phanotate.py"
prodigal_extensions  = [".gff", ".faa", ".fna"]
phanotate_extensions = [".tab", ".faa", ".fna" ]
callers = ["phanotate", "prodigal-11", "prodigal-TAG", "prodigal-TGA"]

pVOG_database_path   = "/home/dani/programs/databases/AllvogHMMprofiles/all_pVOGs.hmm"
pVOG_functional_info = "/home/dani/programs/PhageAnnotation_devel/PhageAnnotation/files/pvogs_annotations.tsv"



rule all:
    input:
        "summary_table_0.txt",
        expand("ORF_calling/prodigal-11/{genome}_prodigal-11{ext}", genome=genomes_id, ext=prodigal_extensions),
        expand("ORF_calling/prodigal-TAG/{genome}_prodigal-TAG{ext}", genome=genomes_id, ext=prodigal_extensions),
        expand("ORF_calling/prodigal-TGA/{genome}_prodigal-TGA{ext}", genome=genomes_id, ext=prodigal_extensions),
        expand("ORF_calling/phanotate/{genome}_phanotate{ext}", genome=genomes_id, ext=phanotate_extensions),
        expand("ORF_calling/tRNAscan/{genome}.trna", genome=genomes_id),
        expand("ORF_calling/{caller}/genes_tables/{genome}_{caller}.genestbl", caller=callers, genome=genomes_id),
        expand("Functional_annotation/VOG/{genome}_{caller}.tablvog", caller=callers, genome=genomes_id),
        expand("Functional_annotation/VOG/genes_tables/{genome}_{caller}.genestbl", caller=callers, genome=genomes_id),
        expand("gggenes_plots/{genome}.png", genome=genomes_id)





rule get_genome_lengths:
    input: input_genomes
    output: "summary_table_0.txt"
    threads: 1
    run:
        with open(output[0], "w") as fout:
            for genome in input:
                record = SeqIO.read(genome, "fasta")
                length = len(record.seq)

                fout.write(f'{record.id}\t{length}\n')

# rule check_circular:
#     input: input: input_path + "/{genome}.fasta"

rule run_prodigal_11:
    input: input_path + "/{genome}.fasta"
    output:
        gff   = "ORF_calling/prodigal-11/{genome}_prodigal-11.gff",
        fasta = multiext("ORF_calling/tmp/{genome}_prodigal-11",
                 ".faa",
                 ".fna"
                 )
    threads: 1
    params: prodigal_path
    shell:
        "{params} -f gff -p meta -i {input} -o {output.gff} -a {output.fasta[0]} -d {output.fasta[1]}"

rule run_prodigal_TAG:
    input: input_path + "/{genome}.fasta"
    output:
        gff   = "ORF_calling/prodigal-TAG/{genome}_prodigal-TAG.gff",
        fasta = multiext("ORF_calling/tmp/{genome}_prodigal-TAG",
                 ".faa",
                 ".fna"
                 )
    threads: 1
    params: prodigal_path
    shell:
        "{params} -f gff -p meta -i {input} -o {output.gff} -a {output.fasta[0]} -d {output.fasta[1]} -TAG Q"

rule run_prodigal_TGA:
    input: input_path + "/{genome}.fasta"
    output:
        gff   = "ORF_calling/prodigal-TGA/{genome}_prodigal-TGA.gff",
        fasta = multiext("ORF_calling/tmp/{genome}_prodigal-TGA",
                 ".faa",
                 ".fna"
                 )
    threads: 1
    params: prodigal_path
    shell:
        "{params} -f gff -p meta -i {input} -o {output.gff} -a {output.fasta[0]} -d {output.fasta[1]} -TGA W"

rule run_phanotate:
    input: input_path + "/{genome}.fasta"
    output: "ORF_calling/phanotate/{genome}_phanotate.tab"
    threads: 1
    params: phanotate_path
    shell:
        "{params} -o {output} {input}"

rule phanotate_ORF_sequences:
    input:
        table  = "ORF_calling/phanotate/{genome}_phanotate.tab",
        genome = input_path + "/{genome}.fasta"
    output:
        multiext("ORF_calling/phanotate/{genome}_phanotate",
                 ".faa",
                 ".fna"
                 )
    threads: 1
    run:
        process_phanotate.extract_ORF_sequences(
                                                input.genome,
                                                input.table
                                               )

rule prodigal_headers:
    input:  "ORF_calling/tmp/{genome}_prodigal-{mode}.{extension}"
    output: "ORF_calling/prodigal-{mode}/{genome}_prodigal-{mode}.{extension}"
    params: "{mode}"
    threads: 1
    run:
        process_prodigal.fix_headers(
                                     input[0],
                                     output[0],
                                     params[0]
                                    )

rule run_tRNAscan2:
    input: input_path + "/{genome}.fasta"
    output: "ORF_calling/tRNAscan/{genome}.trna"
    threads: 1
    shell:
        "tRNAscan-SE -B -o {output} {input}"

rule genes_tables:
    input:
        "ORF_calling/{caller}/{genome}_{caller}.faa",
        "ORF_calling/tRNAscan/{genome}.trna"
    output: "ORF_calling/{caller}/genes_tables/{genome}_{caller}.genestbl"
    threads: 1
    params: "{genome}", "{caller}"
    run:
        tables.genes2table(
                           input[0],
                           input[1],
                           output[0],
                           params[0],
                           params[1]
                          )

rule hmmsearch_pVOG:
    input: "ORF_calling/{caller}/{genome}_{caller}.faa", pVOG_database_path
    output: "Functional_annotation/VOG/{genome}_{caller}.tablvog",
            "Functional_annotation/VOG/{genome}_{caller}.vog"
    threads: 2
    shell:
        "hmmsearch --cpu 2 --tblout {output[0]} {input[1]} {input[0]} > {output[1]}"

rule parse_pVOG:
    input:
        "Functional_annotation/VOG/{genome}_{caller}.tablvog",
        pVOG_functional_info,
        "ORF_calling/{caller}/genes_tables/{genome}_{caller}.genestbl"
    output:
        "Functional_annotation/VOG/genes_tables/{genome}_{caller}.genestbl"
    threads: 1
    run:
        pVOG.pvogs2genetable(input[0],
                             input[1],
                             input[2],
                             output[0]
                            )

rule gggenes_plots:
    input: expand("Functional_annotation/VOG/genes_tables/{{genome}}_{caller}.genestbl", caller=callers)
    output: "gggenes_plots/{genome}.png"
    threads: 1
    script:
        "/home/dani/programs/PhageAnnotation_devel/PhageAnnotation/scripts/plot_gggenes.R"





# rule genes_statistics:
#     input: expand("ORF_calling/{caller}/genes_tables/{genome}_{caller}.genestbl")
#     output: "summary_table_1"
#     threads:1
#     run:
#         tables.

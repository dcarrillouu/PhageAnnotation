#!/usr/bin/env python
# Encoding: utf8

'''

'''

from Bio import SearchIO


def pvogs2genetable(pvog_table, pvog_annot_ref, genes_table, out_table):
    '''

    '''
    # save genes table in a dictionary with KEY=gene_header
    lines_genestable = [line.strip().split("\t") for line in open(genes_table).readlines()]
    genestable = {line[7]:line for line in lines_genestable[1:]}

    # parse pVOGs table
    targets = dict()
    with open(pvog_table) as handle:
        records = SearchIO.parse(handle, "hmmer3-tab")

        for record in records:
            for hit in record.hits:
                toadd = [record.id, hit.evalue, hit.bitscore]
                if hit.id not in targets:
                    targets[hit.id] = [toadd]
                else:
                    targets[hit.id].append(toadd)

    # retain only the best pVOG hit
    # Filter evalue < 1e-5 and bitscore > 20
    targets_filtered = dict()
    for protein, vogs in targets.items():
        #print(hits)
        # move the VOG with the lowest evalue to the top of the list
        sorted_vogs = sorted(vogs, key=lambda vog: float(vog[1]))
        # for hit in sorted_hits:
        #     print("\t", hit)
        if sorted_vogs[0][1] <= 1e-5 and sorted_vogs[0][2] >= 20:
            targets_filtered[protein] = [sorted_vogs[0][0]]


    # Obtain and assign functional information about the VOG term
    annot = dict()
    lines_annot = [line.strip().split("\t") for line in open(pvog_annot_ref).readlines()[1:]]
    for line in lines_annot:
        annot[line[0]] = [line[1]]

    # add pvog information to genes_table
    for gene_header, line in genestable.items():
        if gene_header in targets_filtered:
            line += targets_filtered[gene_header] + annot[targets_filtered[gene_header][0]]
        else:
            line += ["", ""]



    # create list to put output table and sort it
    tow = list()
    for gene_header, line in genestable.items():
        tow.append(line)

    sorted_tow = sorted(tow, key=lambda gene: int(gene[1]))


    # create out_table and write
    with open(out_table, "w") as fout:
        header = lines_genestable[0] + ["VOG", "VOG_annot"]
        fout.write("\t".join(header) + "\n")

        for line in sorted_tow:
            fout.write("\t".join(line) + "\n")

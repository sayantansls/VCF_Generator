#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 21:59:56 2020

@author: sayantan
"""
import time as tm
import utils
import twobitreader
import re
import copy
import sys

HEADERS = ['#CHROM',
           'POS',
           'ID',
           'REF',
           'ALT',
           'QUAL',
           'FILTER',
           'INFO',
           'FORMAT',
           'STRAN-404_S4_DMiSeq02-Run0059']

ENTRY_T = {'#CHROM': '',
           'POS': '',
           'ID': '.',
           'REF': '',
           'ALT': '',
           'QUAL': '.',
           'FILTER': 'PASS',
           'INFO': '.',
           'FORMAT': 'GT:GQ:DP:SR:VR:VA:SB:ABQ:AMQ',
           'STRAN-404_S4_DMiSeq02-Run0059': '1/1:1000.00:608:54.28:54.28:0:0.41:36.11:254.01'}

sep = '\t'

def load_human_genome_sequence(genome_file):
    global genome
    genome = twobitreader.TwoBitFile(genome_file)
    
def load_metadata(vcf_template):
    metadata = list()
    f = open(vcf_template, 'r')
    for line in f.readlines():
        if line.startswith('##'):
            metadata.append(line)
    return metadata
                
#The headers in the genes.tsv file are as follows:
#0 -- Strand_gene_id
#1 -- ChrName
#2 -- Strand
#3 -- Gene_start
#4 -- Gene_end
#5 -- Symbol
#6 -- Entrez_id
                           
def create_gene_chromosome_map(genesfile):
    global gene_chrom_dict
    gene_chrom_dict = dict()
    
    for gene in utils.records_iterator(genesfile):
        if gene['ChrName'] not in gene_chrom_dict:
            gene_chrom_dict[gene['ChrName']] = list()
        gene_chrom_dict[gene['ChrName']].append(gene['Symbol'])
        
#The headers in the input file are as follows:
#0 -- Gene name
#1 -- genomicHGVS

def process_input_file(input_file):
    global variants_list
    variants_list = list()
    for variant in utils.records_iterator(input_file):
        lst = [variant['Gene name'], variant['genomicHGVS']]
        variants_list.append(lst)
    return variants_list

def segregate_variants():
    global subs, others
    subs, others = [list(), list()]
    for variant in variants_list:
        if '>' in variant[1]:
            subs.append(variant)
        else:
            others.append(variant)
    return subs,others

def get_chromosome(gene):
    for k,v in gene_chrom_dict.items():
        if gene in v:
            return k

def del_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    if len(positions) == 1:
        pos = int(positions[0]) - 1
        ref = genome[chrom][pos-1:int(positions[0])].upper()
        alt = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start) - 1
        ref = genome[chrom][pos-1:int(end)].upper()
        alt = genome[chrom][pos-1].upper()
    return pos, ref, alt

def other_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    ins_bases = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    pos = positions[0]
    ref = genome[chrom][int(pos)-1].upper()
    alt = ref + ins_bases
    return pos, ref, alt

def create_substitution_entries(output):
    for variant in subs:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene, gHGVS = [variant[0], variant[1]]
        ref_alt = re.findall(r'[A-Z]', gHGVS)
        position = ''.join(re.findall(r'[0-9]+', gHGVS))
        ref, alt = [ref_alt[0], ref_alt[1]]

        ENTRY['#CHROM'] = get_chromosome(gene)
        ENTRY['POS'] = position
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [ENTRY[i] for i in HEADERS]
        output.write(sep.join(field_values))
        output.write('\n')

def create_non_substitution_entries(output):
    for variant in others:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene, gHGVS = [variant[0], variant[1]]
        chrom = get_chromosome(gene)

        if "del" in gHGVS:
            (pos, ref, alt) = del_handling(gHGVS, chrom)
        else:
            (pos, ref, alt) = other_handling(gHGVS, chrom)

        ENTRY['#CHROM'] = chrom
        ENTRY['POS'] = pos
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [str(ENTRY[i]) for i in HEADERS]
        output.write(sep.join(field_values))
        output.write('\n')

def main(input_file, output_file, genome_file, vcf_template, genesfile):
    print("Start of code:", tm.ctime(tm.time()))

    load_human_genome_sequence(genome_file)
    metadata = load_metadata(vcf_template)

    output = open(output_file, 'w')
    for line in metadata:
        output.write(line)

    output.write(sep.join(HEADERS))
    output.write('\n')

    create_gene_chromosome_map(genesfile)
    variants_list = process_input_file(input_file)

    print("Total variants in the file:", len(variants_list))

    (subs, others) = segregate_variants()

    print("#Substitution variants:", len(subs))
    print("#INDEL variants:", len(others))

    create_substitution_entries(output)
    create_non_substitution_entries(output)

    print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
"""
author : sayantan (sayantan.ghosh@strandls.com)
"""
import re, vcf_utils

def set_human_genome(genomic_build):
    global genome
    genome = vcf_utils.load_human_genome_sequence(genomic_build)

# Gets the pos, ref, alt for a Deletion Variant (g.1124566delG)
def process_deletion_variant(genomicHgvs, chrom):
    positions = re.findall(r'[0-9]+', genomicHgvs)
    if len(positions) == 1:
        pos = int(positions[0]) - 1
        ref = genome[chrom][pos-1:int(positions[0])].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start) - 1
        ref = genome[chrom][pos-1:int(end)].upper()
        
    alt = genome[chrom][pos-1].upper()
    return pos, ref, alt

# Gets the pos, ref, alt for an Insertion Variant (g.14587156insAA)
def process_insertion_variant(genomicHgvs, chrom):
    positions = re.findall(r'[0-9]+', genomicHgvs)
    ins_bases = ''.join(re.findall(r'[A-Z]+', genomicHgvs))
    pos = int(positions[0])
    ref = genome[chrom][pos-1].upper()
    alt = ref + ins_bases
    return pos, ref, alt

# Gets the pos, ref, alt for a Duplication Variant (g.98745415dupA)
def process_duplication_variant(genomicHgvs, chrom):
    positions = re.findall(r'[0-9]+', genomicHgvs)
    refFromHgvs = ''.join(re.findall(r'[A-Z]+', genomicHgvs))
    if len(positions) == 1:
        pos = int(positions[0])
        refFromGenome = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start)
        refFromGenome = genome[chrom][pos-1:int(end)].upper()
        
    ref = vcf_utils.validate_reference({"refFromHgvs": refFromHgvs, "refFromGenome": refFromGenome, "genomicHgvs": genomicHgvs})
    alt = ref * 2
    return pos, ref, alt

# Gets the pos, ref, alt for a Delins Variant (g.54669745_54669748delinsCTGG)
def process_delins_variant(genomicHgvs, chrom):
    positions = re.findall(r'[0-9]+', genomicHgvs)
    if len(positions) == 1:
        pos = int(positions[0])
        ref = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start)
        ref = genome[chrom][pos-1:int(end)].upper()

    alt = ''.join(re.findall(r'[A-Z]+', genomicHgvs))
    return pos, ref, alt
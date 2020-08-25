"""
author : sayantan (sayantan.ghosh@strandls.com)
"""
import re, vcf_utils

# Gets the pos, ref, alt for a Deletion Variant (g.1124566delG)
def deletion_handling(genomicHGVS, chrom, genome):
    positions = re.findall(r'[0-9]+', genomicHGVS)
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
def insertion_handling(genomicHGVS, chrom, genome):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    ins_bases = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    pos = positions[0]
    ref = genome[chrom][int(pos)-1].upper()
    alt = ref + ins_bases
    return pos, ref, alt

# Gets the pos, ref, alt for a Duplication Variant (g.98745415dupA)
def duplication_handling(genomicHGVS, chrom, genome):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    refFromHgvs = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    if len(positions) == 1:
        pos = int(positions[0])
        refFromGenome = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start)
        refFromGenome = genome[chrom][pos-1:int(end)].upper()
        
    ref = vcf_utils.check_ref_hgvs_genome(refFromHgvs, refFromGenome, genomicHGVS)
    alt = ref * 2
    return pos, ref, alt

# Gets the pos, ref, alt for a Delins Variant (g.54669745_54669748delinsCTGG)
def delins_handling(genomicHGVS, chrom, genome):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    if len(positions) == 1:
        pos = int(positions[0])
        ref = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start)
        ref = genome[chrom][pos-1:int(end)].upper()

    alt = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    return pos, ref, alt
# author -- sayantan
# The purpose of this code is as follows:
# 1 -- Variant calling format (VCF) is a type of file in which variations in the human gene are recorded
# 2 -- StrandOmics, the variant prioritization tool, consumes files only in this format
# 3 -- This is a simple utility script which consumes a file of gene names and gHGVS nomeclature and converts it into the corresponfing VCF file format
# 4 -- The input file should have tab separated values of [Gene name] and [Genomic HGVS]
# 5 -- The output will be a complete VCF file including all variant information to be consumed by StrandOmics tool

from __future__ import print_function
import time as tm
import re
import getpass
import os
import zipfile
import twobitreader
import shutil

seconds_start = tm.time()
print("Start of code:",tm.ctime(seconds_start))

username = getpass.getuser()
path = "/home/"+username+"/Downloads/temp_dir"
try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s " % path)

print("Unzipping the files into the directory created")

with zipfile.ZipFile('/home/'+username+'/Downloads/VCF_creation.zip','r') as zip_ref:
    zip_ref.extractall(path)
    
#---------------------------------------------------GENE CHROMSOME DICTIONARY-------------------------------------------------
#The headers in the genes.tsv file are as follows:
#0 -- Strand_gene_id
#1 -- ChrName
#2 -- Strand
#3 -- Gene_start
#4 -- Gene_end
#5 -- Symbol
#6 -- Entrez_id

fhand = open(path+'/genes.tsv','r')
gene_symbol = 5
chromosome = 1
total_list = list()
gene_chrom_dict = dict()
for line in fhand.readlines():
    lst = line.split('\t')[chromosome:gene_symbol+1:4]
    total_list.append(lst)

#Creating a chromosome gene dictionary
result = {}
for k,v in total_list:
    result.setdefault(k,[]).append(v)
gene_chrom_dict = {k: "|".join(v) for k,v in result.items()}

#-----------------------------------------------------PROCESSING INPUT FILE-------------------------------------------------
#The headers in the input file are as follows:
#0 -- Gene name
#1 -- genomicHGVS

gene_name = 0
gHGVS = 1
variants_list = list()
filename = input("Enter the input file name:")
#Taking the entries from the input file
f_input = open(filename,'r')
for line in f_input.readlines():
    lst = line.split('\t')
    lst[gHGVS] = lst[gHGVS].strip()
    variants_list.append(lst)
    
print("Total number of variants in the file:",len(variants_list))

#----------------------------------------------------STORING GENERIC DATA--------------------------------------------------
rs_id = 'rs143035002'
format_id = 'GT:GQ:DP:SR:VR:VA:SB:ABQ:AMQ'
sample = '0/1:1000.00:694:75.22:75.22:0:0.20:37.07:254.00'
with open(path+'/metadata.txt','r') as readfile:
	meta_data = readfile.read()
#-------------------------------------------SEGREGATING THE INPUT FILE VARIANTS------------------------------------------
#Segregating the variants into substitutions and not-substitutions
subs = list()
others = list()

for entry in variants_list:
    if ">" in entry[gHGVS]:
        subs.append(entry)
    else:
        others.append(entry)

print("Number of substitution variants:",len(subs))
print("Number of INDEL variants:", len(others))

def get_chromosome(string): #Function to map the chromosome to the gene name
    for k,v in gene_chrom_dict.items():
        if gene in v:
            return k

genome = twobitreader.TwoBitFile(path+'/hg19.2bit') #Accessing the human genome sequence

#-----------------------------------------HANDLING SUBSTITUTION VARIANTS------------------------------------------------
all_sub_variants = ''
for entry in subs:
    gene = entry[gene_name]
    genomicHGVS = entry[1]
    ref_alt = re.findall(r'[A-Z]',genomicHGVS)
    position = ''.join(re.findall(r'[0-9]+',genomicHGVS))
    ref = ref_alt[0]
    alt = ref_alt[1]
    chrom = get_chromosome(gene)
    final_string = chrom + '\t' + position + '\t' + rs_id + '\t' + ref + '\t' + alt + '\t' + '.' + '\t' + 'PASS' + '\t' + '.' + '\t' + format_id + '\t' + sample
    all_sub_variants += final_string + '\n'


#----------------------------------------HANDLING NON-SUBSTITUTION VARIANTS-------------------------------------------
#Generating the VCF entries for the non-substitution variants
def del_handling(genomicHGVS,chrom):
    positions = re.findall(r'[0-9]+',genomicHGVS)
    if len(positions) == 1:
        pos = int(positions[0]) - 1 
        ref = genome[chrom][pos-1:int(positions[0])].upper()
        alt = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0] ,positions[1]]
        pos = int(start) - 1
        ref = genome[chrom][pos-1:int(end)].upper()
        alt = genome[chrom][pos-1].upper()
    return pos, ref, alt

def other_handling(genomicHGVS,chrom):
    positions = re.findall(r'[0-9]+',genomicHGVS)
    ins_bases = ''.join(re.findall(r'[A-Z]+',genomicHGVS))
    pos = positions[0]
    ref = genome[chrom][int(pos)-1].upper()
    alt = ref + ins_bases
    return pos, ref, alt

all_other_variants = ''
for entry in others:
    gene = entry[gene_name]
    genomicHGVS = entry[gHGVS]
    chrom = get_chromosome(gene)
    if "del" in genomicHGVS:
        (pos, ref, alt) = del_handling(genomicHGVS,chrom)
    else:
        (pos, ref, alt) = other_handling(genomicHGVS,chrom)

    final_string = chrom + '\t' + str(pos) + '\t' + rs_id + '\t' + ref + '\t' + alt + '\t' + '.' + '\t' + 'PASS' + '\t' + '.' + '\t' + format_id + '\t' + sample
    all_other_variants += final_string + '\n'

#---------------------------------------------WRITING OUTPUT FILE----------------------------------------------------
with open('outfile.vcf','w') as writefile:
    writefile.write(meta_data)
    writefile.write(all_sub_variants)
    writefile.write(all_other_variants)

shutil.rmtree(path) #Deleting the extracted files and directory

seconds_end = tm.time()
print("End of code:",tm.ctime(seconds_end))
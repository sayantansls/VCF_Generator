import gzip
import csv

def open_file_for_read(infile):
    if infile.endswith('.gz'):
        return gzip.open(infile, 'rb')
    else:
        return open(infile, 'r')

def open_file_for_write(outfile):
    if outfile.endswith('.gz'):
        return gzip.open(outfile, 'wb')
    else:
        return open(outfile, 'w')

def records_iterator(infile, col_names=None, sep='\t'):
    keys = col_names
    with open_file_for_read(infile) as f:
        for line in f:
            #line = line.strip()
            if not keys:
                if line[0] == '#':
                    line = line[1:]
                keys = [x.strip() for x in line.split(sep)]
            else:
                values = [x.strip() for x in line.split(sep)]
                rec = dict(zip(keys,values))
                rec['__line__'] = line.strip()
                yield rec

def csv_records_iterator(infile):
    reader = csv.DictReader(open(infile, 'r'))
    for rec in reader:
        yield rec

def load_table(infile, key_column_name = None, col_names=None, sep='\t'):
    if key_column_name:
        table = {}
        for rec in records_iterator(infile, col_names, sep):
            key = rec[key_column_name]
            if key in table:
                raise Exception('Duplicate key ' + key)
            table[key] = rec
        return table
    else:
        table = []
        return [rec for rec in records_iterator(infile, col_names)]

# vcf field ids
CHROM   = 0
POS     = 1
ID      = 2
REF     = 3
ALT     = 4
QUAL    = 5
FILTER  = 6
INFO    = 7
FORMAT  = 8

def parse_vcf_info_field(info_str):
    info_map = {}
    fields = info_str.split(';')
    for f in fields:
        if f.find('=') < 0:
            info_map[f] = True
            continue
        k,v = f.split('=')
        info_map[k] = v
    return info_map
    
def vcf_iterator(infile):
    with open_file_for_read(infile) as f:
        in_comments_section = True
        for line in f:
            if in_comments_section:
                if line.startswith('#CHROM'):
                    in_comments_section = False
                continue
            line = line.strip()
            rec = line.split('\t')
            rec[INFO] = parse_vcf_info_field(rec[INFO])
            rec[POS] = int(rec[POS])
            yield rec

def float_equal(f1, f2, allowed_error=0.00001):
    return abs(f1-f2) <= allowed_error

#def print_counts(counts):
 #   for k,v in counts.iteritems():
  #      print '{}\t{:,}\t{}'.format(k, v, v)

def incr(counts, key, n=1):
    if key not in counts:
        counts[key] = n
    else:
        counts[key] = counts[key] + n

def equals_ic(a, b):
    return a.lower() == b.lower()

def inlist_ic(v, inlist):
    lower_v = v.lower()
    for x in inlist:
        if lower_v == x.lower():
            return True
    return False

def same_piped_lists(list1, list2):
    return sorted(list1.split('|')) == sorted(list2.split('|'))

def in_list_as_substr(lst, s):
    for e in lst:
        if s in e:
            return True
    return False

# [ xml
def text(node):
    if node is not None and node.text:
        return node.text.encode('utf-8').strip()
    else:
        return ''
# ] xml

# fields in cancer_terms_mapping_file
# 1          strand_cancer_term: Colorectal Cancer
# 2            bskn_cancer_term: Colorectal Cancer
# 3 bskn_cancer_external_source: ILLUMINA-CUSTOM
# 4     bskn_cancer_external_id: 191711
def load_cancer_terms_mapping(cancer_terms_mapping_file):
    mapping = {}
    for r in records_iterator(cancer_terms_mapping_file):
        key = r['strand_cancer_term']
        value = [r['bskn_cancer_term'], r['bskn_cancer_external_source'], r['bskn_cancer_external_id']]
        mapping[key] = value
    return mapping

# example fields in drug_terms_mapping_file
# fields in drug_terms_mapping_file
# 1            ssdb_drug_term: Everolimus
# 2            bskn_drug_name: everolimus
# 3      bskn_drug_group_name: ILLUMINA-CUSTOM
# 4        bskn_drug_group_id: 484035
def load_drug_terms_map(drug_terms_mapping_file):
    mapping = {}
    for r in records_iterator(drug_terms_mapping_file):
        key = r['ssdb_drug_term']
        value = [r['bskn_drug_term'], r['bskn_drug_group_name'], r['bskn_drug_group_id']]
        mapping[key] = value
    return mapping 

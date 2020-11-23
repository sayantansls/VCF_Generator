"""
author: sayantan (sayantan.ghosh@strandls.com)
Unit Testing for the functions in the script - variant_functions.py
"""

import unittest
import variant_functions
import vcf_utils

test_genomic_build = '37'
(gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(test_genomic_build)
variant_functions.set_human_genome(test_genomic_build)

class TestVariantFunctions(unittest.TestCase):

	def test_process_deletion_variant(self):
		variant = {"gene_name": "BRCA1", "genomicHGVS": "g.41209107delG"}
		expected = {"pos": 41209106, "ref": "TG", "alt": "T"}
		(pos, ref, alt) = variant_functions.process_deletion_variant(variant["genomicHGVS"], 
			vcf_utils.get_chromosome(variant["gene_name"], gene_chrom_dict))
		self.assertEqual(pos, expected["pos"])
		self.assertEqual(ref, expected["ref"])
		self.assertEqual(alt, expected["alt"])


	def test_process_insertion_variant(self):
		variant = {"gene_name": "BRCA1", "genomicHGVS": "g.41244808_41244809insA"}
		expected = {"pos": 41244808, "ref": "C", "alt": "CA"}
		(pos, ref, alt) = variant_functions.process_insertion_variant(variant["genomicHGVS"], 
			vcf_utils.get_chromosome(variant["gene_name"], gene_chrom_dict))
		self.assertEqual(pos, expected["pos"])
		self.assertEqual(ref, expected["ref"])
		self.assertEqual(alt, expected["alt"])

	def test_process_duplication_variant(self):
		variant = {"gene_name": "CFHR5", "genomicHGVS": "g.196963265dupA"}
		expected = {"pos": 196963265, "ref": "A", "alt": "AA"}
		(pos, ref, alt) = variant_functions.process_duplication_variant(variant["genomicHGVS"], 
			vcf_utils.get_chromosome(variant["gene_name"], gene_chrom_dict))
		self.assertEqual(pos, expected["pos"])
		self.assertEqual(ref, expected["ref"])
		self.assertEqual(alt, expected["alt"])

	def test_process_delins_variant(self):
		variant = {"gene_name": "APC", "genomicHGVS": "g.112174576_112174577delinsAT"}
		expected = {"pos": 112174576, "ref": "GC", "alt": "AT"}
		(pos, ref, alt) = variant_functions.process_delins_variant(variant["genomicHGVS"], 
			vcf_utils.get_chromosome(variant["gene_name"], gene_chrom_dict))
		self.assertEqual(pos, expected["pos"])
		self.assertEqual(ref, expected["ref"])
		self.assertEqual(alt, expected["alt"])

if __name__ == '__main__':
	unittest.main()



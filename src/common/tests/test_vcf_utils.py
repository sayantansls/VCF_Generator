"""
author: sayantan (sayantan.ghosh@strandls.com)
Unit Testing for the functions in the script - vcf_utils.py
"""

import unittest
import vcf_utils

test_genomic_build = '37'
genome = vcf_utils.load_human_genome_sequence(test_genomic_build)
(gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(test_genomic_build)

class TestVcfUtils(unittest.TestCase):
	def test_load_genome(self):
		test_cases = [{"position": 41209106, "gene": "BRCA1", "base": "G"}, {"position": 196963264, "gene": "CFHR5", "base": "A"}]

		for test_case in test_cases:
			self.assertEqual(genome[vcf_utils.get_chromosome(test_case["gene"], gene_chrom_dict)][test_case["position"]], test_case["base"])

	def test_get_chromosome(self):
		test_variants = [
			{"gene": "ALK", "chrom": "chr2"}, 
			{"gene": "BRCA1", "chrom": "chr17"}, 
			{"gene": "BRCA2", "chrom": "chr13"}]

		for test_variant in test_variants:
			self.assertEqual(vcf_utils.get_chromosome(test_variant["gene"], gene_chrom_dict)[:3], "chr")
			self.assertTrue(vcf_utils.get_chromosome(test_variant["gene"], gene_chrom_dict)[3:].isdigit())
			self.assertEqual(vcf_utils.get_chromosome(test_variant["gene"], gene_chrom_dict), test_variant["chrom"])

	def test_segregate_variants(self):
		test_variants = [
			{"gene": "BRCA1", "genomicHgvs": "g.41215907C>T"}, 
			{"gene": "BRCA2", "genomicHgvs": "g.32929233C>G"}, 
			{"gene": "BRCA2", "genomicHgvs": "g.32968884_32968885delTT"}, 
			{"gene": "BRCA1", "genomicHgvs": "g.41244808_41244809insA"}]
		(subs, others) = vcf_utils.segregate_variants(test_variants)
		self.assertEqual(len(subs), 2)
		self.assertEqual(len(others), 2)

	def test_validate_gene(self):
		test_genes = {"correct": ["ALK", "BRCA1"], "incorrect": ["POLOP", "GHOSH"]}

		for gene in test_genes["correct"]:
			self.assertFalse(vcf_utils.validate_gene(gene, genes_set))

		for gene in test_genes["incorrect"]:
			self.assertTrue(vcf_utils.validate_gene(gene, genes_set))

	def test_validate_position(self):
		test_cases = {
			"correct" :[{"chrom": "chr17", "pos": 41215907, "ref": "C"}, {"chrom": "chr13", "pos": 32929233, "ref": "C"}],
			"incorrect": [{"chrom": "chr17", "pos": 41215907, "ref": "A"}, {"chrom": "chr13", "pos": 32929233, "ref": "T"}]}

		for test_case in test_cases["correct"]:
			(isPositionInvalid, ref_from_genome) = vcf_utils.validate_position(test_case, genome)
			self.assertFalse(isPositionInvalid)
			self.assertEqual(test_case["ref"], ref_from_genome)

		for test_case in test_cases["incorrect"]:
			(isPositionInvalid, ref_from_genome) = vcf_utils.validate_position(test_case, genome)
			self.assertTrue(isPositionInvalid)
			self.assertNotEqual(test_case["ref"], ref_from_genome)

if __name__ == '__main__':
	unittest.main()
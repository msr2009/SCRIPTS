#\usr\bin\python

"""
translate_fasta.py

translates input FASTA file.
"""

# lookup table for codon translation
def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
             'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
             'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
             'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
             'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
             'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
             'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
             'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
             'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
             'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
             'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
             'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
             'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
             'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
             'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
             'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' }
	return lookup[codon.lower()]

# translate DNA -> amino acid
def translate_sequence(seq):
	translated_seq = ''
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq

def main():

	f_out = open(option.fasta+'.translated', 'w')
		
	handle = open(option.fasta, 'r')
	for record in SeqIO.parse(handle, 'fasta'):
		if record.seq.startswith('ATG') == True:
			print >> f_out, '>' + record.id
			print >> f_out, record.seq.translate()
			
	handle.close()	
	f_out.close()


	
if __name__ == '__main__':
	from optparse import OptionParser
	from math import log
	import sys, os
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'fasta', help = 'path to fasta file')
	parser.add_option('-o', '--output', action = 'store', type = 'string', dest = 'f_output_folder', help = 'path to output folder')
	parser.add_option('--require_atg', action = 'store', type = 'string', dest = 'atg', help = 'require a start codon for each sequence?')	
	(option, args) = parser.parse_args()

	main()
	
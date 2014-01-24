"""
get_synonymous_variants.py

From enrich ratios file, finds all synonymous variants, outputs new ratios file containing only syn vars.

Matt Rich, 8/2013
"""


def main(ratiosfile, wtseq):
	
	#translate wtseq
	wt = translate_sequence(wtseq)
	
	#open outfile
	f_out = open(ratiosfile+".synonymous", "w")
	
	#open ratiosfile
	for line in open(ratiosfile, "r"):
		#print out the header
		if line.startswith("seqID"):
			print >> f_out, line.strip()
		
		#for every other line, check if the translated sequence is the same as the WT translation
		l = line.strip().split('\t')
		if translate_sequence(l[1]) == wt:
			print >> f_out, line.strip()

	f_out.close()
	
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
	try:
		return lookup[codon.lower()]
	except KeyError:
		return 'X'

# translate DNA -> amino acid
def translate_sequence(seq):
	return ''.join([ lookup_codon(seq[b:b+3]) for b in xrange(0, len(seq)-3, 3) ])
	
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'ratios', help = 'path to ratios file')
	parser.add_option('--wtseq', action = 'store', type = 'string', dest = 'wtseq', help = 'wildtype DNA sequence')
	(option, args) = parser.parse_args()

	main(option.ratios, option.wtseq)

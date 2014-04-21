"""
optimize_codons.py

Takes input of a coding DNA sequence or protein sequence, outputs codon-optimized 
version based on codon table input

Matt Rich 4/2014
"""

def main(in_seq, usage_table, sortby):
	# read usage table and make lookup table with ranked codons
	# e.g., 
	# tab = {"A": [['gct','gcc','gca','gcg'],[.65, .25, .07, .03]], etc...}
	
	codons = {}
	r = re.compile("\s?([ACTUGactug]{3})\s+([ACDEFGHIKLMNPQRSTVWY\*])\s+(\d\.\d\d)\s+[\s\d\.]{4}\(.{1,8}\)\s+")
	
	utab = open(usage_table, 'r').read()
	for c in re.findall(r, utab):
		if c[1] in codons: 
			codons[c[1]][0].append(c[0])
			codons[c[1]][1].append(float(c[2]))
			codons[c[1]][2].append( (c[0].upper().count("C") + c[0].upper().count("G"))/3.0 ) 
		else:
			codons[c[1]] = [[c[0]], [float(c[2])], [ (c[0].upper().count("C") + c[0].upper().count("G"))/3.0 ]]

	#sort codons
	s_codons = sort_codons(codons, sortby)

	#choose optimal codon for each position
	return "".join([ codons[x][0][0] for x in in_seq ])

def sort_codons(ctab, sortby):
        for c in ctab:
		if sortby == "usage":
                	ctab[c][0] = [ x for (y,x) in sorted(zip(ctab[c][1], ctab[c][0])) ][::-1]
                	ctab[c][1] = sorted(ctab[c][1])[::-1]
		elif sortby == "gc":
                	ctab[c][0] = [ x for (y,x) in sorted(zip(ctab[c][2], ctab[c][0])) ][::-1]
                	ctab[c][2] = sorted(ctab[c][2])[::-1]
					
# FASTA parser
def readFASTA(f):
	f_out = {}
	k = ''
	val = ''
	
	#find first line that isn't a comment
	comment = True
	while comment == True:
		l = f.next()
		if l.startswith("#") != True:
			#initialize the first key
			k = l.strip()[1:]
			comment = False
	
	while True:
		try:
			l = f.next()
						
			if l.startswith('>'):
				f_out[k] = val
				k = l.strip()[1:]
				val = ''			
			else:
				val += l.strip()

		except StopIteration:
			val += l.strip()
			f_out[k] = val
			return f_out

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
	return ''.join([ lookup_codon(seq[b:b+3]) for b in xrange(0, len(seq), 3) ])

# is the sequence protein or DNA?
def check_sequence(seq, alphabet=set("ACTGU")):
	leftover = set(seq.upper()) - alphabet
	return not leftover
	

if __name__ == "__main__":
	
	from optparse import OptionParser
	import re
	
	parser = OptionParser()
	parser.add_option('--in', action = 'store', type = 'string', dest = 'in_seq', help = 'FASTA file of sequences to be codon-optimized')
	parser.add_option('--usage_table', '-u', action = 'store', type = 'string', dest = 'usage_tab', help = "tab-delimited table of codon usage")
	parser.add_option('--dna', action = 'store_false', dest = 'pro', help = "input is a dna sequence")
	parser.add_option('--protein', action = 'store_true', dest = 'pro', help = "input is a protein sequence", default=True)
	parser.add_option('--sortby', action = 'store', dest = 'sortby', help = "how to choose optimal codon? options: usage, gc", default = "usage")
	(option, args) = parser.parse_args()
	
	fastas = readFASTA(open(option.in_seq,'r'))
	for f in fastas:
		if option.pro == False:
			print ">"+f
			print main(translate_sequence(fastas[f]), option.usage_tab, option.sortby)
		else:
			print ">"+f	
			print main(fastas[f], option.usage_tab, option.sortby)


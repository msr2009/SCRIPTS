"""
genome_viability.py

check annotated variants for hits in essential genes, synthetic lethal pairs, etc...


"""

def main(infile, essential_file, synthetic_file):
	#first populate sets of essential genes and synthetic lethal pairs
	essentials = set()
	if type(essential_file) == set:
		essentials = essential_file
	else:
		for line in open(essential_file, 'r'):
			l = line.strip().split('\t')
			try:
				if l[1].startswith('Y') == True:
					essentials.add(l[1].upper())	
			except IndexError:
				pass
	
	synthetic = set()
	if type(synthetic_file) == set:
		synthetic = synthetic_file
	else:
		for line in open(synthetic_file, 'r'):
			l = line.strip().split('\t')
			if float(l[4]) < -0.5:
				synthetic.add(l[0]+','+l[2])
				synthetic.add(l[2]+','+l[0])
		
	#create a list of possibly-inactivated genes
	#negative substitution matrix score
	bad_genes = []
	for line in open(infile, 'r'):
		if line.startswith('chr') == True:
			continue
		l = line.strip().split('\t')
		if l[5] == 'coding-nonsynonymous':
			if matrix(l[7][0], l[7][-1]) < -1:
				bad_genes.append(l[6])
		elif l[5] == 'splice-site':
			bad_genes.append(l[6])
				
	#cross-reference bad genes list against essentials
	for gene in bad_genes:
		if gene in essentials:
#			print "False: essential - " + gene
			return [False, len(bad_genes)]
			
	
	#create all-by-all combinations of bad_genes and cross-ref against synthetic lethals
	for gene1 in bad_genes:
		for gene2 in bad_genes:
			if gene1+','+gene2 in synthetic:
#				print "False: synthetic - " + gene1 + ', ' + gene2
				return [False, len(bad_genes)]
	
	#if neither essential genes nor synthetic lethal interactions hit, return True (viable)
#	print "True"
	return [True, len(bad_genes)]


#takes single amino acid change and returns matrix score
#this is a BLOSUM80 matrix
def matrix(res1, res2):
	mat = { 'XXX': ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*'],
			'A': [5,-2,-2,-2,-1,-1,-1,0,-2,-2,-2,-1,-1,-3,-1,1,0,-3,-2,0,-6],
			'R': [-2,6,-1,-2,-4,1,-1,-3,0,-3,-3,2,-2,-4,-2,-1,-1,-4,-3,-3,-6],
			'N': [-2,-1,6,1,-3,0,-1,-1,0,-4,-4,0,-3,-4,-3,0,0,-4,-3,-4,-6],
			'D': [-2,-2,1,6,-4,-1,1,-2,-2,-4,-5,-1,-4,-4,-2,-1,-1,-6,-4,-4,-6],
			'C': [-1,-4,-3,-4,9,-4,-5,-4,-4,-2,-2,-4,-2,-3,-4,-2,-1,-3,-3,-1,-6],
			'Q': [-1,1,0,-1,-4,6,2,-2,1,-3,-3,1,0,-4,-2,0,-1,-3,-2,-3,-6],
			'E': [-1,-1,-1,1,-5,2,6,-3,0,-4,-4,1,-2,-4,-2,0,-1,-4,-3,-3,-6],
			'G': [0,-3,-1,-2,-4,-2,-3,6,-3,-5,-4,-2,-4,-4,-3,-1,-2,-4,-4,-4,-6],
			'H': [-2,0,0,-2,-4,1,0,-3,8,-4,-3,-1,-2,-2,-3,-1,-2,-3,2,-4,-6],
			'I': [-2,-3,-4,-4,-2,-3,-4,-5,-4,5,1,-3,1,-1,-4,-3,-1,-3,-2,3,-6],
			'L': [-2,-3,-4,-5,-2,-3,-4,-4,-3,1,4,-3,2,0,-3,-3,-2,-2,-2,1,-6],
			'K': [-1,2,0,-1,-4,1,1,-2,-1,-3,-3,5,-2,-4,-1,-1,-1,-4,-3,-3,-6],
			'M': [-1,-2,-3,-4,-2,0,-2,-4,-2,1,2,-2,6,0,-3,-2,-1,-2,-2,1,-6],
			'F': [-3,-4,-4,-4,-3,-4,-4,-4,-2,-1,0,-4,0,6,-4,-3,-2,0,3,-1,-6],
			'P': [-1,-2,-3,-2,-4,-2,-2,-3,-3,-4,-3,-1,-3,-4,8,-1,-2,-5,-4,-3,-6],
			'S': [1,-1,0,-1,-2,0,0,-1,-1,-3,-3,-1,-2,-3,-1,5,1,-4,-2,-2,-6],
			'T': [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-2,-1,-1,-2,-2,1,5,-4,-2,0,-6],
			'W': [-3,-4,-4,-6,-3,-3,-4,-4,-3,-3,-2,-4,-2,0,-5,-4,-4,11,2,-3,-6],
			'Y': [-2,-3,-3,-4,-3,-2,-3,-4,2,-2,-2,-3,-2,3,-4,-2,-2,2,7,-2,-6],
			'V': [0,-3,-4,-4,-1,-3,-3,-4,-4,3,1,-3,1,-1,-3,-2,0,-3,-2,4,-6],
			'*': [-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,1],
	 		}
	return mat[res1.upper()][mat['XXX'].index(res2.upper())]


if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	
	parser = OptionParser()
	parser.add_option('-f', '--input', action = 'store', type = 'string', dest = 'inputfile', help = 'file with annotated mutations')
	parser.add_option('-e', '--essential', action = 'store', type = 'string', dest = 'essentials', help = 'list of essential genes')
	parser.add_option('-s', '--synthetic', action = 'store', type = 'string', dest = 'synth', help = 'synthetic lethal dataset (constanzo, 2010; lenient cutoff)')
	parser.add_option('-m', '--matrix', action = 'store', type = 'string', dest = 'synth', help = 'substitution matrix for determining mutation effect (default=BLOSUM80')
	(option, args) = parser.parse_args()

	viable = main(option.inputfile, option.essentials, option.synth)


#/usr/bin/python

"""
compare_ems_mutations.py 

takes input of folder containing output files from ems_genome.py; compares mutations in each strain, calculates statistics (overlap in genes, coverage of genome, etc...)

python compare_ems_mutations.py -f <folder containing files>
"""

def main():
	
	# grab list of files from folder
	file_list = os.listdir(option.f_input)
	
	# open each file, populate dictionary
	# of mutations
	mutations = {}
	# and of annotations (counting hits)
	annotations = {}
	# loop through list of files
	for f in file_list:
		# open and read file
		for line in open(f, 'r'):
			l = line.strip().split('\t')
			# add mutations
			mutations[l[0][3:]+':'+l[1]] = l
			# add annotations
			for g in l[3:]:
				gene = g.split(':')
				if gene[0] == 'gene':
					if gene[1] in annotations:
						annotations[gene] += 1
					else: 
						annotations[gene] = 1 
				else:
					if g in annotations:
						annotations[g] += 1
					else:
						annotations[g] = 1

	



if __name__ == '__main__':
	
	from optparse import OptionParser
	import sys, math, os, re

	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'f_input', help = 'path to folder containing mutation files')
	(option, args) = parser.parse_args()
	
	main()
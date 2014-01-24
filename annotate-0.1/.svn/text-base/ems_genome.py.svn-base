#/usr/bin/python

"""
ems_genome.py

Script that simulates performing EMS mutagenesis on a yeast genome. Input genome and expected number of mutations per genome.
Randomly selects positions to mutate.

python ems_genome.py -f <fasta file of genome> -e <expected number of mutations per genome> -i <number of iterations> -o <output folder>
"""

# lookup table for changing roman numerals in annotation file to decimal 
def roman_lookup(rnumeral):
	roman = { 'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
			'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16 }
	return roman[rnumeral.upper()]

# reverse complement 	
def revcomp(sequence):
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	rc = ''
	for base in sequence[::-1]:
		rc += comp[base]
	return rc

# annotate mutations (coding, non-coding, functional element, etc...)	
def annotate_mutations(mutations, annotations):
	#exclude these features from annotation
	exclude_set = set(['chromosome', 'CDS', 'repeat_region', 'region', 'telomere', 'nucleotide_match', 'long_terminal_repeat'])
	#read down lines in annotation file, checking if any mutations fall in annotated regions
	for line in open(annotations, 'r'):
		#return annotations when reach bottom of annotations
		if line.strip() == '##FASTA': 
			#after comparing mutations to annotations, mark all unannotated mutations as 'intergenic'
			for c in mutations:
				for mut in mutations[c]:
					if len(mutations[c][mut]) == 1:
						tmp = mutations[c][mut]
						tmp.append('intergenic')
						mutations[c][mut] = tmp
			return mutations
						
		#skip the first few lines that don't have annotations			
		elif line.startswith('#') != True:
			#split lines
			tokens = line.strip().split('\t')
			#don't worry about some annotations
			if tokens[2] not in exclude_set:
				try:
					chrom = roman_lookup(tokens[0][3:])
				except KeyError:	
					if tokens[0] == 'chrMito':
						chrom = 17
		
				start = int(tokens[3])
				end = int(tokens[4])
				c = 'chr'+str(chrom)
							
				#check if this annotation corresponds to any of the mutations in dictionary
				for mut in mutations[c]:
					# found a mutation in an annotated region?
					if mut >= start and mut <= end:
						tmp = mutations[c][mut]
						if tokens[2] == 'gene':
							gene_name = tokens[8].split(';')[0].split('=')[1]
							tmp.append(tokens[2]+':'+gene_name)
						else:	
							tmp.append(tokens[2])
						mutations[c][mut] = tmp
			
# mutate genome, given reference			
def mutate_genome(reference_genome, expected):
	
	# read fasta genome into list
	genome = []
	for record in SeqIO.parse(reference_genome, 'fasta'):
		genome.append(record.seq)
		
	# calculate some statistics
	# total length of genome
	total_length = 0.0
	gc = 0.0
	mut = 0
	for i in genome:
		gc += i.count('G')
		gc += i.count('C')
	
	# lengths of each chromosome
	chrom_lengths = map(lambda x: len(x), genome)
	total_length = sum(chrom_lengths)	

	# frequency of mutations (mutation/length)
	mutate = expected/gc

	#threshold of non-G>A mutations
	threshold = .1

	# loop through genome, mutate
	# create new genome to store mutations
	mutated_genome = []
	
	for chromosome in genome:
		mutated_chromosome = ''
		for position in range(len(chromosome)):
			# randomly decide if base should be mutated
			trans = random()
			
			#G > X
			if chromosome[position] == 'G':
				#if below mutation threshold, G>A
				if trans > mutate*threshold and trans < mutate:
					mut += 1
					mutated_chromosome += 'A'
				#if in 1% below threshold, any other
				elif trans < mutate*threshold:
					mut += 1
					r1, r2 = random(), random()
					if max([r1,r2]) == r1:
						mutated_chromosome += 'T'
					else:
						mutated_chromosome += 'C'
				else:
					mutated_chromosome += 'G'
			
			#C > X			
			elif chromosome[position] == 'C':
				#if below mutation threshold, C>T
				if trans > mutate*threshold and trans < mutate:
					mut += 1
					mutated_chromosome += 'T'
				#if below threshold, any other
				elif trans < mutate*threshold:
					mut += 1
					r1, r2 = random(), random()
					if max([r1,r2]) == r1:
						mutated_chromosome += 'A'
					else:
						mutated_chromosome += 'G'
				else:
					mutated_chromosome += 'C'
		
			# A > X
			elif chromosome[position] == 'A':
				#if below threshold below mutation
				if trans < mutate*threshold:
					mut += 1
					r1, r2, r3 = random(), random(), random()
					if max([r1,r2,r3]) == r1:
						mutated_chromosome += 'C'
					elif max([r1,r2,r3]) == r2:
						mutated_chromosome += 'G'
					else:
						mutated_chromosome += 'T'
				else:
					mutated_chromosome += 'A'
			
			#T > X
			elif chromosome[position] == 'T':
				#if below threshold below mutation
				if trans < mutate*threshold:
					mut += 1
					r1, r2, r3 = random(), random(), random()
					if max([r1,r2,r3]) == r1:
						mutated_chromosome += 'C'
					elif max([r1,r2,r3]) == r2:
						mutated_chromosome += 'G'
					else:
						mutated_chromosome += 'A'
				else:
					mutated_chromosome += 'T'

		# add new chromosome to genome
		mutated_genome.append(mutated_chromosome)
#	print mut
	return mutated_genome									
				
def main(times, expected, outpath, prefix):
	filenames = set()
	while len(filenames) < times:
		filenames.add(outpath + prefix + str(randint(100000,999999))+'_ems_strain.txt')
	
	while len(filenames) > 0:
		genome = mutate_genome(option.f_input, option.expect)
		f_out = open(filenames.pop(), 'w')
		#print new genome
		for c in range(len(genome)):
			#print fasta header line
			print >> f_out, '>chr'+str(c+1)
			#print sequence lines (60bases/line)
			for j in [genome[c][i:i+60] for i in range(0, len(genome[c]), 60)]:
				print >> f_out, j
		

if __name__ == '__main__':
			
	from optparse import OptionParser
	from Bio import SeqIO
	import sys, math, os, re
	from random import random, randint	

	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'f_input', help = 'path to genome fasta file')
	parser.add_option('-e', '--expect', action = 'store', type = 'float', dest = 'expect', help = 'expected number of mutations per genome')
	parser.add_option('-i', '--iterations', action = 'store', type = 'int', dest = 'iterations', help = 'how many genomes to produce?')
	parser.add_option('-p', '--prefix', action = 'store', type = 'string', dest = 'prefix', help = 'prefix for filename (not path)')	
	parser.add_option('-o', '--output', action = 'store', type = 'string', dest = 'output', help = 'path to output folder')
	(option, args) = parser.parse_args()

	try:
		iterations = int(option.iterations)
	except:
		sys.exit('Iterations cannot be cast as integer')
	try:
		expect = int(option.expect)
	except:
		sys.exit('Expected number of mutations cannot be cast as integer')
	
	if option.prefix == None:
		prefix = ''
	else:
		try:
			prefix = str(option.prefix)
		except:
			prefix = ''
			
	output = option.output
	print iterations, expect, output, prefix

	main(iterations, expect, output, prefix)

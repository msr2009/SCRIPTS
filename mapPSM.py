#!/usr/bin/env python
#Warning: arbitrary, hard-coded limits on ratio below
#This script takes mapunlink files (it can create a ratio between two files as in roundX/roundY) and then generates a fake list of sequences for input into logoPlot or some other sequence logo generator.
#This script written by Doug Fowler (UW)

import sys, os, time, optparse, yapseq, general, pdb

def main():
	
	print time.asctime(time.localtime())
	
	parser = optparse.OptionParser()
	parser.add_option('--type', action = 'store', type = 'string', dest = 'type', help = 'DNA or protein')
	parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to files')
	parser.add_option('--numerator', action = 'store', type = 'string', dest = 'numerator', help = 'input file, is mapUnlink data')
	parser.add_option('--denominator', action = 'store', type = 'string', dest = 'denominator', help = 'input file, is mapUnlink data')
	parser.add_option('--mode', action = 'store', type = 'string', dest = 'mode', help = 'should I generate a fasta file?')
	parser.add_option('--ratio_mode', action = 'store', type = 'string', dest = 'ratio_mode', help = 'mode defining how data is handled.  ratio is the vanilla mode, which generates a simple ratio between two mapunlink files.  num_only just the numerator file, and num_only_absval uses just the numerator file and takes the absolute value of the fitness measurement therein')
	parser.add_option('--threshold', action = 'store', type = 'float', dest = 'threshold', help = 'threshold for accepting the fitness measurement.  for ratios, this will be 1 (e.g. data are included if the ratio is > 1).  for other fitness measurements (slope, for example) ratio could be 0.')
	parser.add_option('--wtseq', action = 'store', type = 'string', dest = 'wtseq', help = 'wt sequence')
	(option, args) = parser.parse_args()

	# Set input path and output file:
	f_numerator = open(option.path + option.numerator, 'U')
	
	if 'num_only' not in option.ratio_mode:
		f_denominator = open(option.path + option.denominator, 'U')
	
	f_outfile = open(option.path + '/weblogo/' + option.mode + '_' + option.numerator, 'w')

	# Set characters for molecule type:
	if option.type == 'protein':
		symbols = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']
		size = 25
	elif option.type == 'DNA':
		symbols = ['A','T','G','C']
		size = 75
		
	# Generate dictionaries of unlinked values:
	num_dict = yapseq.build_unlinked_dict(option.path + option.numerator, 'NA-fix')
	
	if 'num_only' not in option.ratio_mode:
		den_dict = yapseq.build_unlinked_dict(option.path + option.denominator, )
	
	# Obtain the wildtype sequence:
	wildtype = option.wtseq
	
	# Generate ratios dictionaries:
	if option.ratio_mode == 'ratio':
		ratio_dict = {}
		for pos in num_dict:
			ratio_dict[pos] = {}
			for aa in num_dict[pos]:
				if not den_dict[pos][aa] == 0 and not aa == "*":
					ratio = float(num_dict[pos][aa])/den_dict[pos][aa]
					if ratio > option.threshold:
						ratio_dict[pos][aa] = ratio
					else:
						ratio_dict[pos][aa] = 0
				else:
					ratio_dict[pos][aa] = 0
	#Use fitness measurements from the numerator file, without generating a ratio, and take the absolute value of the fitness measurment
	elif option.ratio_mode == 'num_only_absval':
		ratio_dict = {}
		for pos in num_dict:
			ratio_dict[pos] = {}
			for aa in num_dict[pos]:
				if not aa == "*":
					ratio = float(num_dict[pos][aa])
					
					if abs(ratio) > option.threshold:
						ratio_dict[pos][aa] = abs(ratio)
					
					else:
						ratio_dict[pos][aa] = 0
				else:
					ratio_dict[pos][aa] = 0
	
	#Use fitness measurements from the numerator file, without generating a ratio
	elif option.ratio_mode == 'num_only':
		ratio_dict = {}
		for pos in num_dict:
			ratio_dict[pos] = {}
			for aa in num_dict[pos]:
				if not aa == "*":
					ratio = float(num_dict[pos][aa])
					
					if ratio > option.threshold:
						ratio_dict[pos][aa] = ratio
					
					else:
						ratio_dict[pos][aa] = 0
				else:
					ratio_dict[pos][aa] = 0
	
	else:
		sys.exit('Error: choose a valid ratio_mode')
						
	emptyline = {}
	# Generate scaled residue list:
	chars_dict = {}
	for pos in ratio_dict:
		emptyline[pos] = 0
		chars_dict[pos] = {}
		vsum = sum(ratio_dict[pos].values())
		
		for aa in ratio_dict[pos]:
			if vsum != 0:
				chars_dict[pos][aa] = int(10000*float(ratio_dict[pos][aa])/vsum)
			else:
				chars_dict[pos][aa] = 0
				emptyline[pos] = 1
	
	# Fill in to the 1000 mark:
	source_dict = {}
	for pos in chars_dict:
		aas = general.valuesort(chars_dict[pos])
		aas.reverse()
		#print pos, aas[0], aas[1], aas[2], chars_dict[pos][aas[0]], chars_dict[pos][aas[1]], chars_dict[pos][aas[2]]
		diff = 10000-sum(chars_dict[pos].values())
		chars_dict[pos][aas[0]] += diff
		#print sum(chars_dict[pos].values())
		source_dict[pos] = ""
		
		if emptyline[pos] == 0:
			
			for aa in chars_dict[pos]:
				source_dict[pos] += aa*chars_dict[pos][aa] 
			
		elif emptyline[pos] == 1:
				source_dict[pos] = 'X'*10000
	
	# Generate sequences:
	for i in range(0,10000):
		seqID = "mapPSM." + str(i)
		sequence = ""
		for pos in sorted(source_dict.keys()):
			sequence += source_dict[pos][i]
		print >>f_outfile, ">" + seqID + "\n" + sequence + "\n"
	
	print time.asctime(time.localtime())
		
	f_outfile.close()
			
if __name__ == '__main__':
	main()
	
#python YapSeq/Python/mapPSM.py --path YapSeq/Python/output/ --numerator mapunlink_counts_mapcounts_yap6_qc1_PRO_qc2.txt --denominator mapunlink_counts_mapcounts_yap0_qc1_PRO_qc2.txt --mode fasta

#python YapSeq/Python/mapPSM.py --path YapSeq/Python/output/ --numerator mapunlink_wild_mapcounts_yap6_qc1_PRO_qc2.txt --denominator mapunlink_wild_mapcounts_yap0_qc1_PRO_qc2.txt --mode fasta

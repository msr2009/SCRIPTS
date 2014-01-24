#/usr/bin/python

"""
compare_read_depth.py

plots read depths for two Enrich datasets, like input vs. selected

python compare_read_depth.py --file1 FILE1 --file2 FILE2 --plot Y/N --output OUTPUTFOLDER
"""

def main():
	
	from optparse import OptionParser
	from math import log
	import sys, os
	
	parser = OptionParser()
	parser.add_option('--file1', action = 'store', type = 'string', dest = 'file1', help = 'path to 1st file')
	parser.add_option('--file2', action = 'store', type = 'string', dest = 'file2', help = 'path to 2nd file')
	parser.add_option('--plot', action = 'store', type = 'string', dest = 'plot', help = 'use R to plot? (y/n)')
	parser.add_option('--output', action = 'store', type = 'string', dest = 'f_output_folder', help = 'path to output folder')

	(option, args) = parser.parse_args()

	# hash to store data for both files
	# seqID: [data1, data2]
	reads = {}
	
	#open file1 and parse data
	for line in open(option.file1, 'r').readlines()[1:]:
		line = line.strip().split('\t')
		reads[line[0]] = [line[8]]
	#open file2 and parse data
	for line in open(option.file2, 'r').readlines()[1:]:
		line = line.strip().split('\t')
		if line[0] in reads:
			reads[line[0]].append(line[8])

	#print points
	outname = option.f_output_folder + 'read_depths_' + option.file1.split('counts_')[1] + '_' + option.file2.split('counts_')[1]
	f_out = open(outname, 'w')
	print >> f_out, '\t'.join(['seqID', option.file1, option.file2])
	
	for variant in reads:
		if len(reads[variant]) == 2:
			print >> f_out, '\t'.join([variant, reads[variant][0], reads[variant][1]])
				
if __name__ == '__main__':
	main()
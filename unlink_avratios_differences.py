#/usr/bin/python

"""
unlink_avratios_differences.py

takes input and selected unlink_counts file from Enrich, calculates ratios for each position
outputs ratios at each position for each mutation

python unlink_ratios.py --input INPUTLIB --s1 SELECTEDLIB1 --s2 SELECTEDLIB2

"""

def main():

	from optparse import OptionParser
	import sys, os
	from numpy import mean
		
	parser = OptionParser()
	parser.add_option('--input', action = 'store', type = 'string', dest = 'f_input', help = 'path to input unlinked counts file')
	parser.add_option('--s1', action = 'store', type = 'string', dest = 'f_s1', help = 'path to one selected unlinked counts file')
	parser.add_option('--s2', action = 'store', type = 'string', dest = 'f_s2', help = 'path to 2nd selected unlinked counts file')
	(option, args) = parser.parse_args()
	
	# grab header line
	header = []
	
	# store data in two arrays 
	d_input = []
	for line in open(option.f_input, 'r').readlines()[1:]:
		d_input.append(map(lambda x: float(x), line.strip().split()[1:]))
	
	d_sel1 = []
	for line in open(option.f_s1, 'r').readlines()[1:]:
		d_sel1.append(map(lambda x: float(x), line.strip().split()[1:]))
			
	d_sel2 = []
	for line in open(option.f_s2, 'r').readlines()[1:]:
		d_sel2.append(map(lambda x: float(x), line.strip().split()[1:]))
			
	# loop through each line in both arrays, dividing selected by input to calculate enrichment
	sel1_ratios = []
	sel2_ratios = []
	
	for position in range(len(d_input)):
		s2m = mean(d_sel2[position])
		s1m = mean(d_sel1[position])
		sim = mean(d_input[position])
		print '\t'.join([str(position), str((s2m-s1m)/sim)])
		
if __name__ == '__main__':
	main()
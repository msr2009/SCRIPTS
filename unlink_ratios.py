#/usr/bin/python

"""
unlink_ratios.py

takes input and selected unlink_counts file from Enrich, calculates ratios for each position
outputs ratios at each position for each mutation

python unlink_ratios.py --input INPUTLIB --selected SELECTEDLIB --output OUTPUTFOLDER

"""

def main():

	from optparse import OptionParser
	from math import log
	import sys, os
	
	parser = OptionParser()
	parser.add_option('--input', action = 'store', type = 'string', dest = 'f_input', help = 'path to input unlinked counts file')
	parser.add_option('--selected', action = 'store', type = 'string', dest = 'f_selected', help = 'path to selected unlinked counts file')
	parser.add_option('--output', action = 'store', type = 'string', dest = 'f_output_folder', help = 'path to output folder')
	parser.add_option('--log')
	(option, args) = parser.parse_args()
	
	# grab header line
	header = []
	
	# store data in two arrays 
	d_input = []
	for line in open(option.f_input, 'r').readlines():
		if line.startswith('position') == True:
			header = open(option.f_input, 'r').readlines()[0].strip().split()
		else:
			d_input.append(line.strip().split())
	
	d_selected = []
	for line in open(option.f_selected, 'r').readlines()[1:]:
		d_selected.append(line.strip().split())
		
		
	# open output file
	f_out_name = ''
	if option.f_output_folder != None:
		f_out_name = option.f_output_folder + 'unlink_ratios_' + option.f_input.split('_counts_')[1] + '_' + option.f_selected.split('_counts_')[1]
	else:
		f_out_name = 'unlink_ratios_' + option.f_input.split('_counts_')[1] + '_' + option.f_selected.split('_counts_')[1]	
	f_out = open(f_out_name, 'w')
	# print header
	print >> f_out, '\t'.join(header)

	# loop through each line in both arrays, dividing selected by input to calculate enrichment
	for position in range(len(d_input)):
		pos_ratios = [str(position)]
		for aa in range(1,22):
			if d_input[position][aa] == '0':
				pos_ratios.append('0')
			elif d_selected[position][aa] == '0':
				pos_ratios.append('NA')
			else:
				pos_ratios.append( str(log(float(d_selected[position][aa]) / float(d_input[position][aa]), 2)) )
		print >> f_out, '\t'.join(pos_ratios)
	f_out.close()
	
if __name__ == '__main__':
	main()
	
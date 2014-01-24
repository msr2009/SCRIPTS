#/usr/bin/python

"""
extract_position.py

extracts ratios from single position in Enrich ratios file

python extract_positon.py -f RATIOSFILE -p POSITION (in variable region) --output OUTPUTFILE

"""

def main():
	
	from optparse import OptionParser
	from math import log
	import sys, os
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'f_input', help = 'path to file')
	parser.add_option('-p', '--position', action = 'store', type = 'string', dest = 'pos', help = 'position to extract (from variable region)')
	parser.add_option('--output', action = 'store', type = 'string', dest = 'f_output_folder', help = 'path to output folder')

	(option, args) = parser.parse_args()
	
	f_output = open(option.f_input+'.pos'+option.pos, 'w')
	
	for line in open(option.f_input, 'r').readlines()[1:]:
		if option.pos in set( line.strip().split('\t')[0].split('-')[0].split(',') ):
			print >> f_output, line.strip()
#			print line.strip()
	f_output.close()		
	
if __name__ == '__main__':
	main()	
#/usr/bin/python

"""
mutation_histogram.py

outputs histogram of number of mutations per variant in Enrich sequenced library

python -f FILE

"""

def main():
	
	from optparse import OptionParser
	from math import log
	import sys, os
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'f_input', help = 'path to file')
	(option, args) = parser.parse_args()
	
	f_output = open(option.f_input+'_histogram', 'w')
	hist = {}
	
	total = 0
	for line in open(option.f_input, 'r').readlines()[1:]:
		mut_count = int(line.strip().split('\t')[3])
		n = int(line.strip().split()[8])
		if mut_count in hist:
			hist[mut_count][0] += n
			hist[mut_count][1] += 1
		else:
			hist[mut_count] = [n, 1]
		
	for m in sorted(hist.keys()):
		#print '\t'.join([str(m), str(hist[m]), str(hist[m]/total)])
		print >> f_output, '\t'.join([str(m), str(hist[m][0]), str(hist[m][1])]) 
		
	f_output.close()
			
if __name__ == '__main__':
	main()	
"""
read-vs-mutation_quality.py

based on count_histogram.py, reads quality file and randomly selects reads to output for use in a scatterplot

python read-vs-mutation_quality.py -f FILEPATH -m MODE --read_columns READCOLUMNS --mutation_columns MUTCOLUMNS --maxmut MAXIMUMMUTSPERREAD -n N
"""

def main():
	rcols = map(lambda x: int(x), option.rcol.split(','))
	mcols = map(lambda x: int(x), option.mcol.split(','))
	rcolnames = []
	mcolnames = []
	#dictionary to store reads
	reads = {}
	firstline = 1
	
#	file_length = getLength(option.infile)-1
	lines_read = 0
	lines_out = 0
	
	#open file, read lines
	for line in open(option.infile, 'r'):
		line = line.strip().split()
		#skip wildtype lines
		if '100' in line:
			continue
		#skip lines that have more mutations than maxmut
		if True in set([len(x.split(',')) > option.maxmut for x in line]):
			continue
		
		#on the first line, create appropriate number of hist dictionaries in [hists] if not inputted already
		if firstline == 1:
			for i in rcols:
				rcolnames.append(line[i])
			for i in mcols:
				mcolnames.append(line[i])
			firstline = 0
		
		#count data from each line
		elif firstline == 0:
			#a list to keep track of outputs
			output_list = []
			for rc in range(len(rcols)):
				#take appropriate value for each column
				tmp_dat = map(lambda x: int(float(x)), line[rcols[rc]].split(','))
				if option.mode == 'max':
					output_list.append(str(max(tmp_dat)))
				elif option.mode == 'min':
					output_list.append(str(min(tmp_dat)))
				elif option.mode == 'mean':
					output_list.append(str(mean(tmp_dat)))
				elif option.mode == 'median':
					output_list.append(str(median(tmp_dat)))
			for mc in range(len(mcols)):
				#take appropriate value for each column
				tmp_dat = map(lambda x: int(float(x)), line[mcols[mc]].split(','))
				if option.mode == 'max':
					output_list.append(str(max(tmp_dat)))
				elif option.mode == 'min':
					output_list.append(str(min(tmp_dat)))
				elif option.mode == 'mean':
					output_list.append(str(mean(tmp_dat)))
				elif option.mode == 'median':
					output_list.append(str(median(tmp_dat)))
			reads[lines_read] = output_list		
		lines_read += 1
	
	#open outfile and print header line
	f_out = open(option.infile+'-scatterplot-' + option.mode + '-' + str(option.number), 'w')
	print >> f_out, '\t'.join(rcolnames+mcolnames)
	#select reads to output
	for r in sample(range(len(reads)), option.number):
		print >> f_out, '\t'.join(reads[r])
	f_out.close()
	
# get length of file by calling wordcount		
def getLength(infile):
	import os
	#output wordcount command to shell
	os.system('wc -l ' + infile + ' > tmpwc.txt')
	l = 0
	for line in open('tmpwc.txt'):
		l = int(line.strip().split()[0])
	os.remove('tmpwc.txt')
	return float(l)
		
if __name__ == '__main__':

	from numpy import mean, median
	from random import sample
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'path to file')
	parser.add_option('-m', '--mode', action = 'store', type = 'string', dest = 'mode', help = 'output mode: mean,min,max,median')
	parser.add_option('--read_columns', action = 'store', type = 'string', dest = 'rcol', help = 'comma-delimited list of which columns contain read quality data')
	parser.add_option('--mutation_columns', action = 'store', type = 'string', dest = 'mcol', help = 'comma-delimited list of which columns contain mutation quality data')
	parser.add_option('--maxmut', action = 'store', type = 'int', dest = 'maxmut', help = 'max number of mutations in each read to use in counting')
	parser.add_option('-n', action = 'store', type = 'int', dest = 'number', help = 'number of reads to output data')
	(option, args) = parser.parse_args()

	main()
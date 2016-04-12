#/usr/bin/python

"""
find_reads.py

script takes input of file containing readIDs (one per line; can also contain other data, but this data won't be used unless specified),
searches through other files to find reads, prints lines.

python find_reads.py -f <READSLIST> -s <FILETOSEARCH>

"""

def main(infile, search):
	#open readID list, add IDs to set
	readIDs = set()
	for line in open(search, 'r'):
		readIDs.add(line.strip().split('\t')[0])
	
	#open file to search, read lines and check against set
	outfile = open(infile + ".subset", 'w')
	
	for read in read_fastq(infile):
		if read[0].split()[0].lstrip("@") in readIDs:
			print_fastq(read)
	
	outfile.close()
	

if __name__ == '__main__':
	from optparse import OptionParser
	import sys, time
	from fastq_tools import read_fastq, print_fastq

	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'file containing reads to be found')
	parser.add_option('-s', '--search', action = 'store', type = 'string', dest = 'search', help = 'file to search for reads (usually read_fuser output (*_B)')
	(option, args) = parser.parse_args()
	
	print time.asctime()
	main(option.infile, option.search)
	print time.asctime()

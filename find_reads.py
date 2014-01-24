#/usr/bin/python

"""
find_reads.py

script takes input of file containing readIDs (one per line; can also contain other data, but this data won't be used unless specified),
searches through other files to find reads, prints lines.

python find_reads.py -f <READSLIST> -s <FILETOSEARCH> -n <NUMBEROFLINESTOOUTPUT> -o <OVERWRITE?>

"""

def main():
	#open readID list, add IDs to set
	readIDs = set()
	for line in open(option.infile, 'r'):
		readIDs.add(line.strip().split('\t')[0])
	
	#open file to search, read lines and check against set
	outfile = open('/'.join(option.infile.strip().split('/')[0:-1]) +'/' + option.infile.strip().split('/')[-1]+'-'+option.search.strip().split('/')[-1] + '_data', 'w')
	
	firstline = True
	for line in open(option.search, 'r'):
		if firstline == True:
			print >> outfile, line.strip()
			firstline = False
			
		if line.strip().split('\t')[0] in readIDs:
			print >> outfile, line.strip()
			
	outfile.close()
	

if __name__ == '__main__':
	from optparse import OptionParser
	import sys, time

	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'file containing reads to be found')
	parser.add_option('-s', '--search', action = 'store', type = 'string', dest = 'search', help = 'file to search for reads (usually read_fuser output (*_B)')
	parser.add_option('-n', action = 'store', type = 'string', dest = 'lines', help = 'number of lines to output (for use with FASTQ files needing to output 4 lines each', default=0)
	(option, args) = parser.parse_args()
	
	print time.asctime()
	main()
	print time.asctime()
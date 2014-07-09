"""
comments

Matt Rich, DATE
"""

def main(f_in):

	barcodes = {}	
	
	for read in read_fastq(f_in):
		if read[1] in barcodes:
			barcodes[read[1]] += 1
		else:
			barcodes[read[1]] = 1

	for b in barcodes:
		print "\t".join([b, str(barcodes[b])])


if __name__ == "__main__":
	
	from optparse import OptionParser
	from fastq_tools import read_fastq

	parser = OptionParser()
	parser.add_option('--fq', action = 'store', type = 'string', dest = 'fastq', help = "HELP")
	(option, args) = parser.parse_args()
	
	main(option.fastq)	


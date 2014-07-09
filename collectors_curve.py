"""
collectors_curve.py

Given a FASTQ file containing barcode-ID reads, samples the reads and outputs subset FASTQs 
into a supplied folder

Matt Rich, 5/2014

"""

def main(fq, outfolder, diff, n):
	#store all read barcodes
	reads = []
	for record in read_fastq(fq):
		reads.append(record)	
	
	#sample barcodes
	for x in [ y*diff for y in range(1, n+1)]:
		outfile = open(outfolder + ".".join(fq.split(".")[:-1]) + "_" + str(x) + ".fq", "w")
		for i in range(x):
			print_fastq(choice(reads), outfile)
		outfile.close()

if __name__ == "__main__":
	
	from fastq_tools import read_fastq, print_fastq
	from random import choice
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('--fq', action = 'store', type = 'string', dest = 'fq', help = "FASTQ file containing barcode-id'd sequences")
	parser.add_option('-o', action = "store", type = "string", dest = "outfolder", help = "folder to store all subset FASTQs")
	parser.add_option("-n", action = "store", type = "int", dest = "samples", help = "number of sampled FASTQs to generate", default = 10)
	parser.add_option("-d", action = "store", type = "int", dest = "diff", help = "increase in reads per sampled FASTQ", default = 1000000)
	(option, args) = parser.parse_args()
	
	main(option.fq, option.outfolder, option.diff, option.samples)


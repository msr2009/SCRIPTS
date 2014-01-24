"""
splitFASTQ_bySeqID.py

Matt Rich 11/13
"""

def main(fastq, split, outfile):
	#store seqIDs
	seqID = set()
	#open outfile
	f_out = open(outfile, "w")
	
	#read each fastq record
	for record in read_fastq(fastq, None, 250000):
		seqID.add(":".join(record[0].split(":")[0:5]))
	
	print len(seqID)
		
	#read each fastq record and decide if it should be split or not
	for record in read_fastq(split, None, 250000):
		if ":".join(record[0].split(":")[0:5]) in seqID:
			print_fastq(record, f_out)
			
	
if __name__ == "__main__":
	from fastq_util import read_fastq, print_fastq
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option("-f", "--fastq", action = "store", type = "string", dest = "fastq", help = "fastq file with splitting sequences")
	parser.add_option("-s", "--split_fastq", action = "store", type = "string" , dest = "split", help = "fastq file to be split")
	parser.add_option("-o", "--out_fastq", action = "store", type = "string", dest = "outfile", help = "name of fastq file to output")
	
	(option, args) = parser.parse_args()

	main(option.fastq, option.split, option.outfile)
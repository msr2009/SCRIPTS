"""
faqual2fq.py

merges fasta and quality files (a la, phred output)
into a single fastq file (phred+33)

Matt Rich, DATE
"""

from Bio import SeqIO

def main(fasta, qual, phred64):

	#use phred64 if necessary
	p_conv = 33
	if phred64:
		p_conv = 64
	
	#dictionary to store read and quality data
	fq_dict = {}
	
	#read each fasta entry
	for record in SeqIO.parse(open(fasta, 'rU'), "fasta"):
		fq_dict[ record.id ] = [ record.seq ]
 	
	#same for qual file entries
	for record in SeqIO.parse(open(qual, 'rU'), "qual"):
		fq_dict[ record.id ].append( 
			record.letter_annotations["phred_quality"] )

	#print new fastq file
 	for read in fq_dict:
		print "@" + "#".join(read.split())
		print fq_dict[read][0]
		print "+"
		print qual2phred(fq_dict[read][1], p_conv)

def qual2phred(qual_list, p):
	"""returns phred scores converted to ASCII"""
	return "".join([ chr(p+x) for x in qual_list ])

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--fa', action = 'store', type = str,
		dest = 'fasta', help = "fasta file")
	parser.add_argument('--qual', action = 'store', type = str, 
		dest = 'qual', help = "qual file")
	parser.add_argument('--phred64', action = 'store_true', 
		dest = 'phred64', help = "use phred+64 instead of \
		phred+33 (default = False)", default=False)	
	args = parser.parse_args()
	
	main(args.fasta, args.qual, args.phred64)	


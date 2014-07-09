"""
barcode_FASTQ.py

as used to make mapping file for pMR002-SUL1pLibBC (5/2014)

Matt Rich, 5/2014

"""

def main(barcode1, barcode2, sequence, outfile):
	#open fastq files
#	b_in = open(barcode1, "rU")
#	c_in = open(barcode2, "rU")
#	s_in = open(sequence, "rU")
	f_out = open(outfile, "w")

	#for each read in the three fastq files
	for [b, c, s] in read_fastq_multi([barcode1, barcode2, sequence]):
		#find consensus barcode from b and c
		bc = barcode_consensus(b,c)
		if bc != None:
			print_fastq([ b[0]+"%%"+bc, s[1], s[2] ], f_out)

def barcode_consensus(bc1, bc2):
	#if barcodes are the same, return it
	if bc1[1] == bc2[1]:
		return bc1[1]
	elif len(bc1[1]) != len(bc2[1]):
		return None
	else:	
		#if barcodes aren't the same, take highest q-score
		quals = [ a > b for (a,b) in zip(bc1[2], bc2[2]) ]
		con_bc = []
		for x in range(len(bc1[1])):
			if quals[x]:
				con_bc.append(bc1[1][x])
			else:
				con_bc.append(bc2[1][x])
		return "".join(con_bc)
			
if __name__ == "__main__":
	
	from fastq_tools import read_fastq_multi, print_fastq
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('-1', action = 'store', type = 'string', dest = 'bc1', help = 'first barcode read (FASTQ)')
	parser.add_option('-2', action = 'store', type = 'string', dest = 'bc2', help = 'second barcode read (FASTQ)')
	parser.add_option('-s', action = 'store', type = 'string', dest = 'seq', help = 'variant sequence read (FASTQ)')
	(option, args) = parser.parse_args()
	
	main(option.bc1, option.bc2, option.seq, "out.fq")
		
		
		


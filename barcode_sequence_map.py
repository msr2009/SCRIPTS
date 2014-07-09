"""
barcode_sequence_map.py

from fastq containing sequences marked with a barcode in readID,
creates a consensus sequence for each barcode

Matt Rich, 5/2014
"""

def main(fq):
	bc_dict = {}
	
	for record in read_fastq(fq):
		#strip barcode out of ReadID
		bc = record[0].split("%%")[1]
		if bc in bc_dict:
			bc_dict[bc].append( (record[1], record[2]) )
		else:
			bc_dict[bc] = [ (record[1], record[2]) ]

	for x in bc_dict:
		#calculate consensus sequence for each read
		con = calculate_consensus(bc_dict[x])
		if con != None and "N" not in con:
			print "\t".join([x, con])
		else:
			print "\t".join([x, "NA"])

def calculate_consensus(read_list, threshold=0.8):
	if len(set(zip(*read_list)[0])) == 1:
		return read_list[0][0]
	else:
		seq = []
		for x in range( min([ len(x) for (x,y) in read_list ]) ):
			bases = {"A":0, "C":0, "T":0, "G":0, "N":0}
			#sum quality scores for base
			for r in read_list:
				if ord(r[1][x])-33 >= 20:
					bases[r[0][x]] += 1
			#if max base is less than threshold% of total score, trash the barcode
			if sum(bases.values()) == 0:
				return None
			if max(bases.values())/sum(bases.values()) < threshold:
				#print >> sys.stderr, read_list, max(bases.values())/sum(bases.values())
				return None	
			#choose highest-scoring base as correct
			else:
				seq.append( max(bases.iteritems(), key=operator.itemgetter(1))[0] )
		return "".join(seq)

if __name__ == "__main__":
	
	from fastq_tools import read_fastq, print_fastq
	from optparse import OptionParser
	import operator, sys

	parser = OptionParser()
	parser.add_option('--fq', action = 'store', type = 'string', dest = 'fq', help = "FASTQ file containing barcode-id'd sequences")
	(option, args) = parser.parse_args()
	
	main(option.fq)


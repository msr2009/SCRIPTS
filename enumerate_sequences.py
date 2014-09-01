"""
enumerate_seqeunces.py

makes FASTA file containing all single mutations of a sequence

Matt Rich, 8/2014
"""

def main(seq, start, stop, offset):
	bases = ["A","C","T","G"]
	#print wildtype
	print ">WT"
	print seq
	for i in range(start, stop):
		for b in bases:
			if b != seq[i]:
				new_seq = list(seq)
				new_seq[i] = b
				print ">n." + str(i+offset+1) + seq[i] + ">" + b
				print "".join(new_seq) 	

if __name__ == "__main__":
	
	from optparse import OptionParser
	from copy import copy
	
	parser = OptionParser()
	parser.add_option('--seq', action = 'store', type = 'string', dest = 'seq', help = "sequence to be mutated")
	parser.add_option('--start', action = 'store', type = 'int', dest = 'start', help = "position to begin mutatgenesis", default = 0)
	parser.add_option('--stop', action = 'store', type = 'int', dest = 'stop', help = "position to end mutagenesis")
	parser.add_option('--offset', action = 'store', type = 'int', dest = 'offset', help = "numbering offset", default = 0)
	(option, args) = parser.parse_args()
	
	if option.stop == None:
		option.stop = len(option.seq)

	main(option.seq.upper(), option.start, option.stop, option.offset)	


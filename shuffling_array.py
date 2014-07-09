"""
shuffling_array.py

takes alignment fasta as input, outputs array oligos

Matt Rich, 06/2014
"""

def main(aln, n, d, l, five, three):
	#calculate length of sequence in oligo
	o_len = l - len(five) - len(three)
	
	#compile alignment into a list	
	comp_seq = []
	#first sequence, populate list
	for s in SeqIO.parse(open(aln, "r"), "fasta"):

		#initialize comp_seq with longest sequence
		if len(comp_seq) == 0:
			comp_seq = [ [ str(s.seq[x:x+3]).upper() ] for x in range(0, len(s.seq), 3) ]
		#add rest of sequences to comp_seq	
		else:
			for x, y in zip(range(0, len(s.seq), 3), range(len(comp_seq))):
				comp_seq[y].append( str(s.seq[x:x+3]).upper())
		
	#after compiling all sequences, generate oligos by creating random full-length sequences
	oligos = set()

	while len(oligos) < n:
		#make random full-length sequence
		full = "".join([ choice(x) for x in comp_seq ])
		
		#remove dashes 
		full = full.replace("-", "")		
		print full		
		#create oligos by tiling across sequence
		for x in range(0, len(full), d):
			oligos.add(five + full[x:x+o_len] + three)
		
	#print file containing all oligos
	for o in oligos:
		print o

if __name__ == "__main__":
	
	from optparse import OptionParser
	from Bio import SeqIO
	from random import choice

	parser = OptionParser()
	parser.add_option('-a', '--alignment', action = 'store', type = 'string', dest = 'aln', help = "alignment (fasta)")
	parser.add_option('-n', action = 'store', type = 'int', dest = 'n', help = "number of oligos to output", default = 10000)
	parser.add_option('-l', '--length', action = 'store', type = 'int', dest = 'length', help = "oligo length (including 5' and 3' sequences)", default = 100)
	parser.add_option('-5', '--five-prime', action = 'store', type = 'string', dest = 'five', help = "5' adaptor sequence", default = "")
	parser.add_option('-3', '--three-prime', action = 'store', type = 'string', dest = 'three', help = "3' adaptor sequence", default = "")
	parser.add_option('-d', '--tile-dist', action = 'store', type = 'int', dest = 'tiledist', help = "tiling distance (default = 20bp)", default = 20)
	
	(option, args) = parser.parse_args()
	
	main(option.aln, option.n, option.tiledist, option.length, option.five, option.three)	


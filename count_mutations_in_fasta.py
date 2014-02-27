"""
count_mutations_in_fasta.py

compares two FASTA files, counts number of mutatations between them
"""

def main(fasta1, fasta2):
	f1 = SeqIO.to_dict(SeqIO.parse( open(fasta1,"rU"), "fasta"))
	f2 = SeqIO.to_dict(SeqIO.parse( open(fasta2,"rU"), "fasta"))
	
	mutations = 0
	length = 0	

	for c in f1:
		print c
		mutations += hammingDistance(f1[c], f2[c])
		length += len(f1[c])	
	
	print mutations, length, float(mutations)/float(length)

def hammingDistance(seq1, seq2):
	return len(seq1) - sum([ seq1[x] == seq2[x] for x in range(len(seq1))])

if __name__ == "__main__":
	from optparse import OptionParser
	from Bio import SeqIO
	
	parser = OptionParser()
	
	parser.add_option("--f1", action = "store", type = "string", dest = "fasta1", help = "first FASTA file")
	parser.add_option("--f2", action = "store", type = "string", dest = "fasta2", help = "second FASTA file")
	(option, args) = parser.parse_args()

	main(option.fasta1, option.fasta2)

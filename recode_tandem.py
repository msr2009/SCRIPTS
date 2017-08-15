"""
recode_tandem.py

recodes tandem proteins (like tdMCP) so that I can synthesize them in gblocks

alternates codons and changes to either optimal or non-optimal codon

Matt Rich, 2017
"""

def main(seq, codon_file):
	#read codon table
	codon_table = {}
	lookup = {}
	for line in open(codon_file, "r"):
		l = line.strip().split()
		lookup[l[1]] = l[0]
		if l[0] in codon_table:
			codon_table[l[0]][0].append(l[1])
			codon_table[l[0]][1].append(float(l[4]))
		else:
			codon_table[l[0]] = [[l[1]],[float(l[4])]]

	sorted_table = {}
	for a in codon_table:
		sorted_table[a] = [y for (x,y) in sorted(zip(codon_table[a][1], codon_table[a][0]))]

	new_seq = ""
	for i in range(0,len(seq), 3):
		if i % 2 == 0:
			new_seq += recode_codon(seq[i:i+3], sorted_table, lookup)
		else:
			new_seq += seq[i:i+3]

	print ">new_seq1"
	print new_seq

	new_seq = ""
	for i in range(0,len(seq), 3):
		if i % 2 == 1:
			new_seq += recode_codon(seq[i:i+3], sorted_table, lookup)
		else:
			new_seq += seq[i:i+3]

	print ">new_seq2"
	print new_seq

def recode_codon(old_codon, codons, lookup):
	if old_codon.upper() != "ATG" and old_codon.upper() != "TGG":
		used_codon = codons[lookup[old_codon]].index(old_codon)
		if used_codon == 0:
			return codons[lookup[old_codon]][1]
		else:
			return codons[lookup[old_codon]][0]
	else:
		return old_codon

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--seq', action = 'store', type = str, dest = 'seq', 
		help = "sequence for recoding")
	parser.add_argument('--codons', action = 'store', type = str, dest = 'codon', 
		help = "codon usage table")
	args = parser.parse_args()
	
	main(args.seq, args.codon)	


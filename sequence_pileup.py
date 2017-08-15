"""
sequence_pileup.py
 
from sequence (raw reads, one per line), makes pileup for each variant, outputs
sorted file, one position per line

Matt Rich, 02/2017
"""

def main(reads, protein):
	firstread = True
	bases = {"protein": ['A', 'C', 'D', 'E', 'F', 'G', 'H', 
						 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 
						 'R', 'S', 'T', 'V', 'W', 'Y', '*', 'X'],
			 "dna": ['A', 'C', 'G', 'T', 'N']}
		
	for line in open(reads, "r"):
		l = line.strip()
		if firstread:
			print "firstread"
			if protein == "protein":	
				dat = {i:{x:0 for x in bases["protein"]} for i in range(len(l))}
			else: 
				dat = {i:{x:0 for x in bases["dna"]} for i in range(len(l))}
			firstread = False
		for x in range(len(l)):
			dat[x][l[x]] += 1

	#print header
	if protein == "protein":
		print "\t".join(["pos"] + bases["protein"])
		for i in dat:
			print "\t".join([str(i)] + [ str(dat[i][j]) for j in bases["protein"]]) 
				
	else:
		print "\t".join(["pos"] + bases["dna"])
		for i in dat:
			print "\t".join([str(i)] + [ str(dat[i][j]) for j in bases["dna"]]) 

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--reads', '-r', action = 'store', type = str, dest = 'reads', 
		help = "file containing reads (one per line)")
	parser.add_argument('--protein', action = 'store_true', dest = 'protein', 
		help = "reads are protein sequence", default=False)
	
	args = parser.parse_args()
	
	if args.protein:
		main(args.reads, "protein")
	else:
		main(args.reads, "dna")


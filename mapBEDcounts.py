"""
mapBEDcounts.py

takes a BED file with "unlinked" mutations and counts occurrences of each mutation. Outputs in BED format.

"""

def main(bed, pileup):
	
	#open pileup file, store in dictionary
	piles = {}
	for line in open(pileup, "r"):
		l = line.strip().split('\t')
	 	piles[int(l[0])] = float(l[1])
	
	#open outfile
	f_out = open(bed+".counts", "w")
	#print header
	print >> f_out, "\t".join(['chr', 'start', 'stop', 'ref', 'obs', 'count', 'pop_ratio'])
	
	#read BED file in dictionary
	bed_counts = {}
	
	firstline = True
	seq = ""
	mut_type = ""
	
	for line in open(bed, 'r'):
		if "start" in line:
			continue
			
		l = line.strip().split('\t')
		
		pos = int(l[1])
		
		if firstline == True:
			seq = l[0]
			firstline = False
		
		if pos in bed_counts:
			if l[6] != "S":
				if l[6].startswith("I"):
					mut_type = "I-" + l[4]
				else:
					mut_type = l[6]
				if mut_type not in bed_counts[pos]:
					bed_counts[pos][mut_type] = 1
				else:
					bed_counts[pos][mut_type] += 1
			else:
				bed_counts[pos][l[4].upper()] += 1
		else:
			if l[6] != "S":
				if l[6].startswith("I"):
					mut_type = "I-" + l[4]
				else:
					mut_type = l[6]
				bed_counts[pos] = bed_counts[pos] = {"A":0, "C":0, "G":0, "T":0, mut_type:1, "ref":l[3]}
			else:		
				bed_counts[pos] = {"A":0, "C":0, "G":0, "T":0, "ref":l[3]}
				bed_counts[pos][l[4].upper()] += 1	
		
	#print out BED containing counts for all mutations
	for p in sorted(bed_counts.keys()):
		for b in sorted(bed_counts[p].keys()):
			if b != "ref":
				if b.startswith("I"):
					print >> f_out, '\t'.join([seq, str(p), str(p+len(b.split('-')[1])-1), bed_counts[p]["ref"], b, str(bed_counts[p][b]), str(bed_counts[p][b]/piles[p])])
				else:
					print >> f_out, '\t'.join([seq, str(p), str(p+1), bed_counts[p]["ref"], b, str(bed_counts[p][b]), str(bed_counts[p][b]/piles[p])])
				
if __name__ == "__main__":
	
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option("-b", "--bed", action = "store", type = "string", dest = "bed", help = "BED file containing mutations")
	parser.add_option("-p", "--pileup", action = "store", type = "string", dest = "pileup", help = "file containing pileup (read coverage at each base)")
	(option, args) = parser.parse_args()

	main(option.bed, option.pileup)
 
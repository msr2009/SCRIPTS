"""
mapBEDratios.py

Takes two BED files (preferably the output of annotate_orfsequences.py) and calculates ratios between files. 
Prints both ratios file and matrix for heatmaps
"""

def main( f1, f2, pad, minreads, matrix ):
	
	#read files
	counts1 = readBEDcounts(f1, pad, True, minreads)
	#pull out header from counts1
	h = counts1.pop("header")
	
	counts2 = readBEDcounts(f2, pad, False, minreads)
		
	#dictionary to store matrix information
	mat = {}
	#fill in dictionary with NA's by default
	for i in range(1, max([ (int(x.split(':')[0])/3)+1 for x in counts1.keys() ])):
		mat[i] = { "A":"NA", "C":"NA", "D":"NA", "E":"NA", "F":"NA", "G":"NA", "H":"NA", "I":"NA", "K":"NA", "L":"NA",
					"M":"NA", "N":"NA", "P":"NA", "Q":"NA", "R":"NA", "S":"NA", "T":"NA", "V":"NA", "W":"NA", "Y":"NA", "*":"NA" }
	
	#print header
	print "\t".join( h[0:5] + [f1.split('/')[-1]+"-count", f1.split('/')[-1]+"-frac", f2.split('/')[-1]+"-count", f2.split('/')[-1]+"-frac", "ratio" ] + h[7:] )
	
	#compare files
	for m in counts1:
		if m in counts2:
			pre = counts1[m][1][0:6]
			suf = counts1[m][1][7:]
			print '\t'.join( pre + [str(counts1[m][0]), counts2[m][1][5], str(counts2[m][0]), str(log(counts2[m][0]/counts1[m][0],2))] + suf )
			#add data to matrixdict
			if counts1[m][1][8] == "NS" or counts1[m][1][8] == "N":
				mat[int(counts1[m][1][7])][counts1[m][1][6][-1]] = str(log(counts2[m][0]/counts1[m][0],2))

	#open matrix outfile and print matrix
	f_mat = open(matrix, "w")
	print >> f_mat, '\t'.join(["position"] + sorted(mat[1].keys()))
	for pos in sorted(mat.keys()):
		print >> f_mat, '\t'.join([str(pos)] + [ mat[pos][x] for x in sorted(mat[pos].keys()) ])
	f_mat.close()
							
def readBEDcounts( bed, pad, header, minreads ):
	#read in files, keep track of total number of counts in each
	#key -> "start:obs" 
	#value -> [pop_freq, splitline] 
	counts = {}
	
	p = 0
	if pad == True:
		p = 1
			
	for line in open(bed, 'r'):
		#skip header 
		if line.startswith("chr"):
			if header == True:
				counts["header"] = line.strip().split('\t')
			continue
		
		#get counts for each mutation
		l = line.strip().split('\t')
		if int(l[5]) >= minreads:
			counts[l[1]+":"+l[4]] = [float(l[6])+p, l]

	return counts
		
	
if __name__ == "__main__":
	
	from math import log
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option("--f1", action = "store", type = "string", dest = "bed1", help = "BED file containing mutations (usually input)")
	parser.add_option("--f2", action = "store", type = "string", dest = "bed2", help = "BED file containing mutations (usually selected)")
	parser.add_option("--pad", action = "store_true", dest = "pad", help = "add 1 to all counts? default=False", default=False)
	parser.add_option("--min-reads", action = "store", dest = "minreads", type = "int", help = "read count threshold under which mutations are discarded. default=1", default=1)
	parser.add_option("--matrix", action = "store", dest = "matrix", type = "string", help = "filename for matrix output (default = 'matrix.txt'", default="matrix.txt")

	(option, args) = parser.parse_args()

	main(option.bed1, option.bed2, option.pad, option.minreads, option.matrix)
	
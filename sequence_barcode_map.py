"""
sequence_barcode_map.py

counts barcodes for given sequence, based on mapping file

Matt Rich, 05/2014
"""

def main(m):
	seqs = {}
	for line in open(m, "rU"):
		l = line.strip().split('\t')
		if l[1] != "NA":
			if l[1] in seqs:
				seqs[l[1]].append(l[0])
			else:
				seqs[l[1]] = [l[0]]

	for s in seqs:
		print "\t".join( [s, ",".join(seqs[s]), str(len(seqs[s]))] )

if __name__ == "__main__":
	
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('-m', action = 'store', type = 'string', dest = 'barcodemap', help = "tab-delimited barcode map")
	(option, args) = parser.parse_args()
	
	main(option.barcodemap)	

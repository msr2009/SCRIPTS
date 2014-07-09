"""
barcode_map_to_FASTA.py

converts tab-delimited barcode map to FASTA file containing all barcoded sequences

Matt Rich, 5/2014
"""

def main(m):
	m_out = open(m+".fa", "w")
	for line in open(m, "rU"):
		l = line.strip().split("\t")
		if l[1] != "NA":	
			print >> m_out, ">"+l[0]
			print >> m_out, l[1]

if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('-m', action = 'store', type = 'string', dest = 'barcodemap', help = "tab-delimited barcode map")
	(option, args) = parser.parse_args()
	
	main(option.barcodemap)


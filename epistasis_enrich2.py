"""
epistasis_enrich2.py

Matt Rich, 12/14
"""

def main(singles, doubles, mincounts):

	#read single mutant file and extract data
	s = {}
	for line in open(singles, "r"):
		if line.startswith("sequence") != True:	
			l = line.strip().split("\t")
			if int(l[1]) >= mincounts:
				s[l[0]] = [ l[-1], l[1] ] 
	
	#print a header for our output
	print "\t".join([ "sequence", "A", "B", "counts", "countsA", "countsB", "slope", "slopeA", "slopeB"])
	
	#read double mutant file and calculate epistasis
	for line in open(doubles, "r"):
		if line.startswith("sequence") != True:	
			l = line.strip().split("\t")
			if int(l[1]) >= mincounts:
				mutA, mutB = [ x.strip() for x in l[0].split(",") ]	
				try:
					out = [ l[0], mutA, mutB, l[1], s[mutA][1], s[mutB][1], l[-1], s[mutA][0], s[mutB][0] ]
					print "\t".join( [ str(x) for x in out ] )
				except KeyError:
					continue

if __name__ == "__main__":
	
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('--singles', action = 'store', type = 'string', dest = 'singles', help = "file containing single mutant data")
	parser.add_option('--doubles', action = 'store', type = 'string', dest = 'doubles', help = "file containing double mutant data")
	parser.add_option('--min-counts', action = 'store', type = 'int', dest = 'mincounts', help = "minimum input read counts", default = 1)
	(option, args) = parser.parse_args()
	
	main(option.singles, option.doubles, option.mincounts)	


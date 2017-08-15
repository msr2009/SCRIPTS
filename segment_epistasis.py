"""
segment_epistatis.py

finds all double mutants between two ranges of sites, and reports epistasis values for each

written for use with output from Enrich v2

Matt Rich, 08/2014
"""

def main(singles, doubles, r1, r2, epi, mincounts):
	#open singles file, read all singles
	s1 = {}
	s2 = {}
	for line in open(singles, "r"):
		l = line.strip().split('\t')
		if l[0] != "_wt" and l[0] != "sequence":
			#figure out position from HGVS_to_PosMut
			pos, ref, obs = HGVS_to_PosMut(l[0])
			#add to either s1/s2 if in either range
			if pos in r1:
				s1[l[0]] = float(l[-2])
			elif pos in r2: 
				s2[l[0]] = float(l[-2])
	
	#open doubles file, find mutants that contain sites in _both_ ranges
	for line in open(doubles, "r"):	
		l = line.strip().split("\t")
		if l[0] == "sequence":
			print "\t".join(l + ["epistasis"])
		elif l[0] != "_wt":
			#again, find position with HGVS_to_PosMut
			dpos = l[0].split(", ")
			#find mutants with sites in both ranges
			if dpos[0] in s1 and dpos[1] in s2:
				#calculate epistasis
				print "\t".join(l + [str(epistasis(s1[dpos[0]], s2[dpos[1]], float(l[-2]), epi))] )
			elif dpos[1] in s1 and dpos[0] in s2:
				#calculate epistasis
				print "\t".join(l + [str(epistasis(s1[dpos[1]], s2[dpos[0]], float(l[-2]), epi))] )	

def epistasis(s1, s2, d, method):
	if method == "additive":
		return d - s1 - s2
	elif method == "multiplicative":
		return d - s1 * s2

def HGVS_to_PosMut(hgvs):
	if hgvs != "_wt":
		posmut = re.match('n\.(-?\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs).groups()
		return [int(posmut[0]), posmut[1], posmut[2]]
	else:
		pass

if __name__ == "__main__":
	
	from optparse import OptionParser
	import re

	parser = OptionParser()
	parser.add_option('--singles', action = 'store', type = 'string', dest = 'singles', help = "single mutant data file")
	parser.add_option('--doubles', action = 'store', type = 'string', dest = 'doubles', help = "double mutant data file")
	parser.add_option('--r1', action = 'store', type = 'string', dest = 'r1', help = "position range 1")
	parser.add_option('--r2', action = 'store', type = 'string', dest = 'r2', help = "position range 2")
	parser.add_option('--method', action = 'store', type = 'string', dest = 'method', help = "epistasis type (additive, multiplicative)", default = "additive")
	parser.add_option('--min-counts', action = 'store', type = 'int', dest = 'mincounts', help = "minimum input counts for mutant", default=0)
	(option, args) = parser.parse_args()
	
	if option.r2 == None:
		option.r2 = option.r1

	r1 = [ int(x) for x in option.r1.split(",") ]
	r2 = [ int(x) for x in option.r2.split(",") ]
	
	if option.method not in ["additive", "multiplicative"]:
		raise TypeError("Unexpected epistasis type")

	main(option.singles, option.doubles, range(r1[0],r1[1]), range(r2[0],r2[1]), option.method, option.mincounts)	


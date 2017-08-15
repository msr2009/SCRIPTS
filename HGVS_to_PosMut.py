"""
comments

Matt Rich, DATE
"""
import re
import pandas as pd

def HGVS_to_PosMut(hgvs, pandas=False):
	if hgvs != "_wt":
		posmut = re.match('n\.(-?\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs).groups()
		if pandas:
			return pd.Series([int(posmut[0]), posmut[1], posmut[2]])
		else:
			return [int(posmut[0]), posmut[1], posmut[2]]
	else:
		pass

if __name__ == "__main__":
	
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = "enrich output file")
	(option, args) = parser.parse_args()
	
	for line in open(option.infile, "r"):
		l = line.strip().split('\t')	
		
		if line.startswith("sequence"):	#header line
			print "\t".join(l + ["pos", "ref", "obs"])
			continue
		
		posmut = HGVS_to_PosMut(l[0])
		print "\t".join(l + [str(posmut[0]), posmut[1], posmut[2]])			
	

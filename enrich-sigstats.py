"""
enrich-sigstats.py

calculates mann-whitney U test for a set of barcodes compared to wildtype

Matt Rich, 08/2014
"""

def main(f, q, d):
	#open file, calculate mann-whitney U test for each variant
	for line in open(f, "r"):
		l = line.strip().split('\t')
		#if line is header, print header + p-value
		if l[0] == "sequence":
			print "\t".join(l[0]+["p-value"])
		
		#then calculate mann-whitney for each variant and 
		#and print p-value
		

def get_query_distribution(f, mincounts):
	#open query enrich output (e.g., wildtype barcodes)
	#create np.array to store all slope values
	df = pd.read_csv(dfile, sep="\t", header=0)
	return np.array(df[df['count.0']>=mincounts]['slope'])
	
def get_variant_distribution(f, mincounts):
	#open barcodes for all variants, store all slope values
	#based on the variant sequence 
	d = {}
	for line in open(f, 'r'):
		if line.startswith("sequence") == True:
			continue	
		
		l = line.strip().split("\t")
		if l[1] >= mincounts:
			if l[-1] in d:
				np.append(d[l[-1]], d[l[-2]])
			else:
				d[l[-1]] = np.array(d[l[-2]])
	
	return d

if __name__ == "__main__":
	
	from optparse import OptionParser
	import pandas as pd
	import numpy as np
	from scipy.stats import mannwhitneyu, ttest_ind
	
	parser = OptionParser()
	parser.add_option('--OPT', action = 'store', type = 'string', dest = 'DEST', help = "HELP")
	(option, args) = parser.parse_args()
	
	main()	


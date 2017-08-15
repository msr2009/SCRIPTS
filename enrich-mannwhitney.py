"""
enrich-sigstats.py

calculates significance for a set of barcodes compared to wildtype

Matt Rich, 08/2014
"""

def main(f, q, d):
	#open file, calculate mann-whitney U test for each variant
	for line in open(f, "r"):
		l = line.strip().split('\t')
		#if line is header, print header + p-value
		if l[0] == "sequence":
			print "\t".join(l+["p-value"])
			continue
		#then calculate statistic for each variant and 
		#and print p-value
		if l[0] in d:
			print d[l[0]]
			print "\t".join( l + [str(mannwhitneyu(d[l[0]], q)[1])] )
		else:
			print "\t".join( l + [ 'NA' ] )

def get_query_distribution(f, mincounts):
	#open query enrich output (e.g., wildtype barcodes)
	#create np.array to store all slope values
	df = pd.read_csv(f, sep="\t", header=0)
	return np.array(df[df['count.0']>=mincounts]['slope'])
	
def get_variant_distribution(f, mincounts):
	#open barcodes for all variants, store all slope values
	#based on the variant sequence 
	d = {}
	for line in open(f, 'r'):
		if line.startswith("sequence") == True:
			continue	
		l = line.strip().split("\t")
		if float(l[1]) >= mincounts:
			if l[-1] in d:
				d[l[-1]] = np.append(d[l[-1]], np.array(l[-2]))
			else: 
				d[l[-1]] = np.array(l[-2])
	return d

if __name__ == "__main__":
	
	from optparse import OptionParser
	import pandas as pd
	import numpy as np
	from scipy.stats import mannwhitneyu, ttest_ind
	
	parser = OptionParser()
	parser.add_option('-f', action = 'store', type = 'string', dest = 'variants', help = "variants file to calculate significance")
	parser.add_option('-q', '--query', action = 'store', type = 'string', dest = 'query', help = "distribution to compare against")
	parser.add_option('-b', '--barcodes', action = 'store', type = 'string', dest = 'barcodes', help = "barcodes file")
	parser.add_option('--mincounts', action = 'store', type = 'int', dest = 'mincounts', help = "minimum input counts", default = 10)
	(option, args) = parser.parse_args()
	
	main(option.variants, get_query_distribution(option.query, option.mincounts), get_variant_distribution(option.barcodes, option.mincounts))	


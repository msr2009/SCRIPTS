"""
pandas_pwm.py

using pandas (generally as enrich-v2 output) to output a MEME file consisting of 
log-odds for a given motif

Matt Rich, 7/2014
"""

def main(dfile, datcol, mincounts, f, range, wt_score, wt_seq, missing):
	
	df = preprocess(pd.read_csv(dfile, sep="\t", header=0), datcol, mincounts)
	#open output file
	f_out = open(f, "w")
	#print a couple lines to f_out
	print >> f_out, "MEME version 4\n"
	print >> f_out, "ALPHABET= ACGT\n"
	print >> f_out, "strands: +\n"
	print >> f_out, " ".join(["MOTIF", f, datcol])	

	#make a pivot table of all the data 
	### KINDA HACKY -- THIS ASSUMES THERE'S ONE DATAPOINT FOR EACH POSITION
	data = df.pivot("pos", "mut", "dat")
	#make a series out the wtseq
	wt = pd.Series(list(wt_seq), index=data.index.values)		
	
	#add wt_score to each position's wt
	for i in data.index.values:
		data.loc[i][wt.loc[i]] = wt_score
	#then convert all NaNs to .01 (or something so there's no missing data
	data[data.isnull()] = missing
	
	#extract data range
	pwm = data.loc[range[0]:range[1]]
	print >> f_out, "log-odds matrix: "
	pwm =  pwm.apply(np.log2)*100
#	pwm = pwm.apply(np.round)
	#output pwm 
#	pwm.to_csv(f_out, sep="\t", na_rep="0", header=True, index_label="PO")
	pwm.to_csv(f_out, sep="\t", na_rep="0", header=False, index=False)

	f_out.close()

def preprocess(df, datcol, mincounts):
	#split HGVS into pos, mut
	new_df = df[df["count.0"] >= mincounts]["sequence"].transpose().apply(func=HGVS_to_PosMut)
	new_df.columns = ["pos", "ref", "mut"]
#	print new_df.shape
	new_df["dat"] = df[datcol]
	return new_df
	
def HGVS_to_PosMut(hgvs):
	if hgvs != "_wt":
		posmut = re.match('n\.(-?\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs).groups()
		return pd.Series([int(posmut[0]), posmut[1], posmut[2]])
	else:
		pass
	
if __name__ == "__main__":
	
	from optparse import OptionParser
	
	import pandas as pd
	import numpy as np
	from numpy.linalg import norm	
	import re
	
	parser = OptionParser()
	parser.add_option('--data', action = 'store', type = 'string', dest = 'data', help = "tab-delimited dataframe")
	parser.add_option('--dna', action = 'store_true', dest = 'dna', help = "is DNA sequence?", default = False)
	parser.add_option('--min-counts', action = 'store', type = 'int', dest = 'mincounts', help = 'minimum number of input counts', default = 0)
	parser.add_option('--datcol', action = 'store', type = 'string', dest = 'datcol', help = 'column to plot')
	parser.add_option('--offset', action = 'store', type = 'int', dest = 'offset', help = 'axis offset', default = 0)
	parser.add_option('--range', action = 'store', type = 'string', dest = 'range', help = 'comma-delimited range (start,end) for PWM extraction')
	parser.add_option('--wt-score', action = 'store', type = 'float', dest = 'wt', help = 'wildtype score', default = 1)
	parser.add_option('--wt-seq', action = 'store', type = 'string', dest = 'wt_seq', help = "wildtype sequence")
	parser.add_option('--missing', action = 'store', type = 'float', dest = 'missing', help = "value to use for missing data (before taking log)", default=1.0)
	(option, args) = parser.parse_args()
		
	r = [int(x) for x in option.range.split(',')]
	f_out = ".".join(option.data.split(".")[0:-1]) + "-".join(option.range.split(",")) + ".pwm"
	print "output printed to " + f_out
	main(option.data, option.datcol, option.mincounts, f_out, r, option.wt, option.wt_seq.upper(), option.missing )

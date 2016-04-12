"""
pandas_pwm.py

using pandas (generally as enrich-v2 output) to output a MEME file consisting of 
log-odds for a given motif

Matt Rich, 7/2014
"""

def main(dfile, datcol, mincounts, f, range, wt_score, wt_seq, missing, output, raw_out):
	
	df = preprocess(pd.read_csv(dfile, sep="\t", header=0), datcol, mincounts)

	#make a pivot table of all the data 
	### KINDA HACKY -- THIS ASSUMES THERE'S ONE DATAPOINT FOR EACH POSITION
	data = df.pivot("pos", "mut", "dat")

	#make a series out the wtseq
	wt = pd.Series(list(wt_seq), index=data.index.values)		
	#add wt_score to each position's wt
	for i in data.index.values:
		data.loc[i][wt.loc[i]] = wt_score
	#then convert all NaNs to missing data value
	data[data.isnull()] = missing
	#extract data range
	pwm = data.loc[range[0]+1:range[1]+1]
	print pwm
	
	#print raw data from pwm
	if raw_out != None:
		pwm.to_csv(raw_out, sep="\t", header=True, index_label="PO")

	#output as pwm or activity matrix
	if output == "meme":
		meme_out = open(f+".meme", "w")
		print >> meme_out, "MEME version 4\n"
		print >> meme_out, "ALPHABET= ACGT\n"
		print >> meme_out, "strands: +\n"
		print >> meme_out, " ".join(["MOTIF", f, datcol])	
		print >> meme_out, "log-odds matrix: "
		pwm2 =  pwm.multiply(100)
	#	pwm2 = pwm.apply(np.round)
		#output pwm 
		pwm2.to_csv(meme_out, sep="\t", na_rep="0", header=False, index=False)
		meme_out.close()
	
	if output == "pwm":
		pwm_out = open(f+".pwm", "w")	
		print >> pwm_out, "ID " + f + " " + datcol
		print >> pwm_out, "S. cerevisiae"
		pwm2 = pwm.apply(lambda x: 2**x)
		pwm2.to_csv(pwm_out, sep="\t", na_rep="0", header=True, index_label="PO", index=True)
		print >> pwm_out, "XX"
		print >> pwm_out, "//"
		pwm_out.close()	
		

def preprocess(df, datcol, mincounts):
	#split HGVS into pos, mut
	new_df = df[df["count.0"] >= mincounts]["sequence"].transpose().apply(func=HGVS_to_PosMut)
	new_df.columns = ["pos", "ref", "mut"]	
	if datcol == "slope40":
		print "extrapolating to 40 generations"
		new_df["dat"] = df["score"]*40.0+df["intercept"]
#		new_df["dat"][ new_df["dat"] < 0.01 ] = 0.01
	else:
		new_df["dat"] = df[datcol]
#	print new_df.loc[-193]
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
	parser.add_option('--output_format', action = 'store', type = 'string', dest = 'output', help = "output for meme or activity?", default = "meme")
	parser.add_option('--raw', action = 'store', type = 'string', dest = 'rawout', help = 'output raw data to file', default = None)
	(option, args) = parser.parse_args()
		
	r = [int(x) for x in option.range.split(',')]
	f_out = ".".join(option.data.split(".")[0:-1]) + "-".join(option.range.split(","))
	print "output printed to " + f_out

	main(option.data, option.datcol, option.mincounts, f_out, r, option.wt, option.wt_seq.upper(), option.missing, option.output, option.rawout )

"""
apply_drift.py

normalizes to drift and (optionally) wt slopes for enrich data

Matt Rich, 8/2014
"""

def main(expts, drifts, wt_norm, out):
	
	#read dataframes
	xdf = pd.read_csv(expts, header=0, sep="\t", index_col=0)

	dft = pd.read_csv(drifts, header=0, sep="\t", index_col=0)
	new_dft = pd.DataFrame(index=dft.index)
	new_dft["drift"] = dft["score"]

	#merge drift column into xdf
	new_xdf = pd.merge(xdf, new_dft, left_index=True, right_index=True, how="outer")

	#perform drift normalization
	new_xdf["drifted_score"] = new_xdf["score"] - new_xdf["drift"]

	#define wt_score
	wt_score = 0
	if wt_norm:
		try:
			wt_score = new_xdf.loc["_wt"]["drifted_score"]
		except KeyError:
			sys.stderr.write("_wt missing. No wildtype normalization performed.")
			pass
	
	#wt normalization
	new_xdf["norm_score"] = new_xdf["drifted_score"] - wt_score
	
	#print new dataframe
	new_xdf.sort("norm_score",ascending=False).to_csv(out, sep="\t", header=True, na_rep="NaN")

if __name__ == "__main__":
	
	import pandas as pd
	import sys
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('--drift', action = 'store', type = 'string', dest = 'drift', help = "variants file from drift control")
	parser.add_option('--expt', action = 'store', type = 'string', dest = 'expt', help = "variants file from experiment")
	parser.add_option('--wt_norm', action = 'store_true', dest = 'wt_norm', help = "normalize to wt slope (default=No)", default=False)
	parser.add_option('--out', action = 'store', type = 'string', dest = 'out', help = "output file for new dataframe")
	(option, args) = parser.parse_args()
	
	main(option.expt, option.drift, option.wt_norm, option.out)	


"""
enrich_score_by_position.py

Matt Rich, 07/14
"""

def main(data, datcol, out, countfilt):
	#read data into pandas dataframe
	df = pd.read_csv(data, header=0, sep="\t")
	
	#filter df
	df = df[df["count.0"] >= countfilt]

	#make dataframe of parsed HGVS
	df_new = df["sequence"].transpose().apply(func=HGVS_to_PosMut)
	df_new.columns = ["pos", "ref", "mut"]
	#extract data column of interest
	
	df_new["dat"] = df[datcol]

	#take mean of each mutation at a position
	df_piv = df_new.pivot('pos', 'mut', 'dat')
	#take mean of rows
	df_mean = df_piv.mean(1)
	#take sd of rows
	df_sd = df_piv.std(1)
	#plot
	fig, ax = plt.subplots()
#	ax.errorbar(df_piv.index, df_mean, yerr = df_sd, capsize = 0, color="r")
	ax.bar(df_piv.index, df_mean, color='r', align='center')
	ax.set_ylabel("Average wt-normalized slope")
	xtickloc = mpl.ticker.MultipleLocator(10)
	xtickform = mpl.ticker.ScalarFormatter()
	ax.xaxis.set_major_locator(xtickloc)
	ax.xaxis.set_major_formatter(xtickform)
	ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())

	plt.savefig("test.pdf")
	
def HGVS_to_PosMut(hgvs):
	if hgvs != "_wt":
		posmut = re.match('n\.(-?\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs)
		return pd.Series((int(posmut.groups()[0]), posmut.groups()[1], posmut.groups()[2]))

if __name__ == "__main__":
	
	from optparse import OptionParser
	import re
	import pandas as pd
	import numpy as np
	
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt

	parser = OptionParser()
	parser.add_option('--data', action = 'store', type = 'string', dest = 'data', help = "tab-delimited dataframe")
	parser.add_option('--out', '-o', action = 'store', type = 'string', dest = 'figout', help = 'output filename', default='fig.pdf')
	parser.add_option('--datcol', action = 'store', type = 'string', dest = 'datcol', help = 'column to plot')
	(option, args) = parser.parse_args()
	
	main(option.data, option.datcol, option.figout, 10)	


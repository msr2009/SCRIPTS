"""
comments

Matt Rich, DATE
"""

def main(df, wt_seq, dna, figout):
	
	#set the WT sequence of your region of interest
	WT_sequence=wt_seq
	if dna == False:
		aa = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'R', 'H', 'K', 'D', 'E', '*']
	else:
		aa = ['A', 'C', 'G', 'T']
	
	#This next line gives you the position of each WT amino acid in your "aa" list. We use this later to mark the WT positions in the heatmap. We reverse all our lists so that the heatmap is drawn top to bottom (default is bottom to top)
	WT_aa_index = [aa.index(i) for i in WT_sequence][::-1] 
	
	#Columns will denote position (first 100 positions here)
	column_labels = range(len(WT_sequence)) 
	row_labels = aa
	
	#make empty matrix to store heatmap values (filled with zeroes here)
	data = np.zeros((len(column_labels),len(row_labels))) 
	
	#fill in matrix by iterating over your dataframe
	for index,row in df.iterrows(): 
		amino = row.mut
		pos = int(row.pos)
		data[pos-1,aa.index(amino)] = row.dat
	
	#mask the 0 values so that we can color them grey later
	data = np.ma.masked_values(data, 0)[::-1] 
	#Now that we have the matrix with positions in columns and amino acids in rows, we can draw the actual heatmap 
	#set heatmap color scheme
	my_cmap = cm.get_cmap('RdBu_r') 
	#set color for masked array entries
	my_cmap.set_bad('grey') 
	
	fig,axarr = plt.subplots(1,1)
	fig.set_size_inches(20,42)
	heatmap = axarr.imshow(data, interpolation='none', aspect='auto', cmap=my_cmap) #vmin sets worst score, vmax sets best
	axarr.set_xticklabels(row_labels, minor=False)
	axarr.set_yticklabels(column_labels[::-1], minor=False)
	axarr.set_xticks(np.arange(data.shape[1]), minor=False)
	axarr.set_yticks(np.arange(data.shape[0]), minor=False)
	axarr.tick_params(axis='both', which='major', labelsize=26)
	axarr.invert_yaxis()
	axarr.set_frame_on(False)
	axarr.grid(False)
	#Mark the WT positions
	for i in column_labels:
		#hatch types: '-', '+', 'x', '\\', '/', '*', 'o', 'O', '.'
		axarr.add_patch(mpl.patches.Rectangle((WT_aa_index[i]-0.5,i-0.5), 1, 1, fill=False, color="white")) 
		#axarr.add_patch(mpl.patches.Wedge((WT_aa_index[i],i), .25, 0, 360, width=0.25, color="white")) #position, outer size, rotation, inner size
	fig.colorbar(heatmap)
#	fig.title("Score", fontsize=60)

	plt.savefig(figout)

def preprocess(df, datcol, takelog=True):
	takelog=False
	#split HGVS into pos, mut
	new_df = df["sequence"].transpose().apply(func=HGVS_to_PosMut)
	new_df.columns = ["pos", "ref", "mut"]
	if takelog:
		new_df["dat"] = df[datcol].apply(np.log2)
	else:
		new_df["dat"] = df[datcol]
	return new_df
	
def HGVS_to_PosMut(hgvs):
	if hgvs != "_wt":
		posmut = re.match('n\.(\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs)
		return pd.Series(posmut.groups())
	else:
		pass

if __name__ == "__main__":
	
	from optparse import OptionParser
	
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	
	import pandas as pd
	import numpy as np	
	import re
	
	parser = OptionParser()
	parser.add_option('--data', action = 'store', type = 'string', dest = 'data', help = "tab-delimited dataframe")
	parser.add_option('--dna', action = 'store_true', dest = 'dna', help = "is DNA sequence?", default = False)
	parser.add_option('--wt', action = 'store', type = 'string', dest = 'wt_seq', help = "wildtype sequence")
	parser.add_option('--out', '-o', action = 'store', type = 'string', dest = 'figout', help = 'output filename', default='fig.pdf')
	parser.add_option('--datcol', action = 'store', type = 'string', dest = 'datcol', help = 'column to plot')
	(option, args) = parser.parse_args()
	
##	print preprocess(pd.DataFrame())
	main(preprocess(pd.read_csv(option.data, sep="\t", header=0), option.datcol), option.wt_seq.upper(), option.dna, option.figout)	



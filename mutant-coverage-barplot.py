"""
mutant-coverage-barplot.py

Matt Rich, 7/2014
"""

def main(mat):
	#read mat, convert all sequences and store positions in list
	pos = {}
	for line in open(mat, "r"):
		hgvs = line.strip().split('\t')[0]
		if hgvs != "sequence" and hgvs != "_wt":
			#split sequence HGVS, extract positions
			for p in [ HGVS_to_PosMut(x)[0] for x in hgvs.split(",") ]:
				if p in pos:	
					pos[p] += 1
				else:
					pos[p] = 1
		
	fig, ax = plt.subplots()
	fig.set_size_inches(12,9)
	ax.bar(pos.keys(), pos.values(), color="blue", linewidth=0, edgecolor="blue", log=True)
	ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
	plt.xlim( (min(pos.keys()), max(pos.keys())) )
	plt.savefig(".".join(mat.split(".")[0:-1])+".bar.pdf")	

def HGVS_to_PosMut(hgvs):
	if hgvs != "_wt":
		posmut = re.match('\W?n\.(-?\d+)([ACDEFGHIJKLMNPQRSTUVWY]+)>([ACDEFGHIJKLMNPQRSTUVWY]+)', hgvs).groups()
		return [int(posmut[0]), posmut[1], posmut[2]]
	else:
		pass

def axis_labels(pos, offset):
	return str( pos-offset )

if __name__ == "__main__":
	
	from optparse import OptionParser

	import matplotlib
	matplotlib.use("Agg")
	import matplotlib.pyplot as plt
	
	import re
	
	parser = OptionParser()
	parser.add_option('-m', action = 'store', type = 'string', dest = 'mat', help = "enrich output matrix")
	(option, args) = parser.parse_args()
	
	main( option.mat ) 	


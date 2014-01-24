"""
matrix_boxplots.py

makes boxplots out of a matrix, based on a given set of row names

Matt Rich, 12/13
"""

def main(infile):
	import matplotlib.pyplot as plt
	
	#read file
	dat = []
	names = []
	for line in open(infile, 'r'):
		l = line.strip().split('\t')
		names.append(l[0])
		dat.append([float(x) for x in l[1:]])	
	#print dat
	
	
	fig, ax1 = plt.subplots(figsize=(10,6))
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	
	#make boxplots
	plt.figure()
	bp = plt.boxplot(dat)
	
	plt.xticks(range(len(names)), names)
	
	plt.savefig("plot1.pdf", format="pdf")
	
if __name__ == "__main__":
	
	from optparse import OptionParser
		
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'tab-delimited matrix of data')
	(option, args) = parser.parse_args()
	
	main(option.infile)	
	


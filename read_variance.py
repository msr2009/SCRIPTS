"""
read_variance.py

Given a ratios file as input, calculates mean, SD, and variance for variants binned by read counts.

Matt Rich 8/2013
"""

def main(ratios, counts, bins):
	
	#make dictionary to store ratios
	rdict = {}
	for b in bins:
		rdict[b] = []
		
	#make dictionary to store counts
	cdict = {}
	for line in open(counts,'r'):
		if line.startswith("seqID") != True:
			l = line.strip().split("\t")
			cdict[l[0]] = l[-1]#determineBin(int(l[-1]), bins)
	
	#populate dictionary of ratios by counts
	for line in open(ratios,'r'):
		if line.startswith("seqID") != True:
			l = line.strip().split('\t')
			#rdict[ cdict[l[0]] ].append(float(l[7]))
			print "\t".join([l[0], cdict[l[0]], l[7]])

	
	#output mean, sd, var
#	print "\t".join(["counts_bin", "n", "mean", "sd", "var"])	
#	for b in bins:
#		print "\t".join( [ str(b), str(len(rdict[b])), str(mean(rdict[b])), str(std(rdict[b])), str(var(rdict[b])) ] )			
		
def determineBin(count, b):
	b = b[::-1]
	if count >= max(b):
		return max(b)
	elif count <= min(b):
		return min(b)
	else:
		for x in range(len(b)):
			if count <= b[x] and count > b[x+1]:
				return b[x]
	
if __name__ == "__main__":
	from optparse import OptionParser
	from numpy import mean, std, var
	
	
	parser = OptionParser()
	parser.add_option('-r', '--ratios', action = 'store', type = 'string', dest = 'ratios', help = 'path to ratios file')
	parser.add_option('-c', '--counts', action = 'store', type = 'string', dest = 'counts', help = 'path to counts file')
	parser.add_option('--bins', action = 'store', type = 'string', dest = 'bins', 
			help = 'comma-delimited list of bins for read counts (default = 1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50)', default="1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50")
	(option, args) = parser.parse_args()

	main(option.ratios, option.counts, [int(x) for x in option.bins.split(",")])

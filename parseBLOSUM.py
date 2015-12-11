"""
comments

Matt Rich, DATE
"""

def parseBLOSUM(f):

	aa = []
	blosum_dict = {}	
	
	for line in open(f, "r"):
		if line.startswith("#"):
			continue
		elif line.strip().endswith("*"):
			aa = line.strip().split()	
			continue
		else:
			l = line.strip().split()
			blosum_dict[l[0]] = dict(zip(aa, [ int(x) for x in l[1:]]))
	
	return blosum_dict		

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--blosum', action = 'store', type = str, 
						dest = 'blosum', help = "file with blosum matrix")
	args = parser.parse_args()
	
	tmp = parseBLOSUM(args.blosum)	
	for x in tmp:
		print x, tmp[x]

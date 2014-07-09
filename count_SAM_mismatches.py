"""
count_SAM_mismatches.py

Matt Rich, 5/2014
"""

def main(sam):
	for line in open(sam, "rU"):
		if line.startswith("@") == True:
			continue
		l = line.strip().split('\t')
		if l[2] == "*":
			continue
		bc = l[0]
		m = 0
		d = 0
		for x in set(l):
			if x.startswith("XM:"):
				m = x.split(":")[2]
			elif x.startswith("XO:"):
				d = x.split(":")[2]
		print "\t".join([bc, m, d])		
		
if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('-s', action = 'store', type = 'string', dest = 'sam', help = "sam file")
	(option, args) = parser.parse_args()
	
	main(option.sam) 

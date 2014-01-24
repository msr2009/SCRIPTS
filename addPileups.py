"""
addPileups.py

Takes input of list of pileup files (from parseCross_matchOutput.py); adds across files for each position

Matt Rich, 9/2013
"""

def main(plist):
	all_pileups = {}
	
	for f in open(plist,"r").readlines():
		#read file into dictionary
		for line in open(f.strip(), 'r').readlines():
			l = line.strip().split('\t')
			try:
				all_pileups[int(l[0])] += int(l[1])
			except KeyError:
				all_pileups[int(l[0])] = int(l[1])
				
	#print combined pileup file
	for x in sorted(all_pileups.keys()):	
		print str(x) + "\t" + str(all_pileups[x])
		
		
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', action = 'store', type = 'string', dest = 'pileups', help = 'file containing list of pileups, one per line')

	(option, args) = parser.parse_args()

	main(option.pileups)

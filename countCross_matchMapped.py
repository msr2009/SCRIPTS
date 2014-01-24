"""
countCross_matchMapped.py

"""

def main(cmfile):
	cm = open(cmfile,'r')
	
	while True:
		line = cm.readline()
		if line == "":
			return
		
		l = line.strip().split('\t')
		counts =  str(len(l[1].split('-')[-1].split(',')))
		print "\t".join( l + [counts] )  

if __name__ == "__main__":
	from optparse import OptionParser	
	parser = OptionParser()
	parser.add_option('-f', action = 'store', type = 'string', dest = 'inputfile', help = 'mapped cross_match output')
	(option, args) = parser.parse_args()

	main(option.inputfile)		

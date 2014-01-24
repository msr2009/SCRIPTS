"""
Phd1ToQUAL.py

Takes input of a phd.1 file (ABI quality format) and converts it to a QUAL file for use with cross_match or other aligners

Matt Rich 4/2013
"""

def main(phd1):
	#flag for finding sequencing qualities
	in_read = False
	#list to hold quality values
	qual = []
	
	#loop through phd1 file
	p = open(phd1,'r')
	print ">" + phd1.split('_')[0] + ".i" 
	while True:
		line = p.readline()
		
		if line.strip() == "BEGIN_DNA":
			in_read = True
			continue
		elif line.strip() == "END_DNA":
			break
		
		if in_read == True:
			l = line.strip().split()
			qual.append(l[1])
		
	#print the output to std_out
	print " ".join(qual)
	p.close()
	
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', action = 'store', type = 'string', dest = 'phd1', help = 'phd.1 file')
	(option, args) = parser.parse_args()
	
	main(option.phd1)
	
			
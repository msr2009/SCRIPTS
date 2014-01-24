"""
getFileLength.py


"""

import sys

def getFileLength(f, lowercase):
	
	l = ""
	
	if f == "STDIN":
		l += sys.stdin.read()
	else:
		l += open(f,'r').read()
		
	full_line = "".join(l.split('\n'))
	
	if lowercase == True:
		full_line = list(full_line)
		for i in range(len(full_line)):
			if full_line[i].islower() == True:
				full_line[i] = ""	
		full_line = "".join(full_line)
	
	return full_line.strip()
	
if __name__ == "__main__":
	from optparse import OptionParser
		
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'inputfile', help = 'file to find length of. Use "STDIN" to read from stdin.')
	parser.add_option('--lowercase', action = 'store_true', dest = 'lowercase', help = "should I exclude lowercase characters? def=False", default=False)
	(option, args) = parser.parse_args()

	print len(getFileLength(option.inputfile, option.lowercase))
"""
changeFASTQLastChar.py

Script changes last character of FASTQ header line

Matt Rich 3/2013
"""

def main(f, char):
	#open file
	f_in = open(f, 'r')
	
	while True:
		read = [ f_in.readline().strip() for i in range(4) ]
		if read == ['','','','']:
			break
		
		read[0] = read[0][0:-1] + char
		print '\n'.join(read)
			
if __name__ == "__main__":
	from optparse import OptionParser	
	parser = OptionParser()
	parser.add_option('-f', action = 'store', type = 'string', dest = 'inputfile', help = 'FASTQ file')
	parser.add_option('-c', action = 'store', type = 'string', dest = 'newchar', help = 'new character (to replace specified one)')	
	(option, args) = parser.parse_args()

	main(option.inputfile, option.newchar)
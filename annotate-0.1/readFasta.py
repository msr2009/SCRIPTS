"""
readFASTA.py

function to read FASTA files into a dictionary. Replaces BioPython's SeqIO.parse(x,'fasta'),
although without some functionality.

Matt Rich 2012
"""

def readFASTA(f):
	f_out = {}
	k = ''
	val = ''
	
	#find first line that isn't a comment
	comment = True
	while comment == True:
		l = f.next()
		if l.startswith("#") != True:
			#initialize the first key
			k = l.strip()[1:]
			comment = False
	
	while True:
		try:
			l = f.next()
						
			if l.startswith('>'):
				f_out[k] = val
				k = l.strip()[1:]
				val = ''			
			else:
				val += l.strip()

		except StopIteration:
			val += l.strip()
			f_out[k] = val
			return f_out
						
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'fasta', help = 'FASTA file to be parsed')
	(option, args) = parser.parse_args()
	
	print readFASTA(open(option.fasta))
	


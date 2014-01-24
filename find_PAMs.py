"""
find_PAMs.py

Given a FASTA sequence, find all possible 20mers+PAM (NGG) on either strand.  Output list.

Matt Rich 5/13
"""

def main(fasta, flank5, flank3):
	#parse fasta file
	f = readFASTA(open(fasta,'r'))

	#loop through sequence, find all NGGs
	for i in range(20,len(f)):
		if f[i+1:i+3] == "GG":
			print flank5 + f[i-20:i] + flank3
	
	#loop through RC'd sequence, find all NGGs
	for i in range(20,len(complement(f[::-1]))):
		if f[i+1:i+3] == "GG":
			print flank5 + f[i-20:i] + flank3

def complement(base):
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'a':'T', 't':'A', 'c':'G', 'g':'C', 'n':'N'}
	base = list(base)
	return ''.join( [ comp[b] for b in base ] )	

# FASTA parser
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
			f_out = val
			return f_out

	
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'fasta', help = 'FASTA file of sequence')
	parser.add_option('-3', '--flank3', action = 'store', type = 'string', dest = 'flank3', help = "3' flanking sequence for output", default = "")
	parser.add_option('-5', '--flank5', action = 'store', type = 'string', dest = 'flank5', help = "5' flanking sequence for output", default = "")
	(option, args) = parser.parse_args()

	main(option.fasta, option.flank5, option.flank3)
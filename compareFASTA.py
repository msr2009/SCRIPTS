"""
compareFASTA.py

takes input of reference and observed FASTAs, outputs BED file with mismatches

Matt Rich 2/2014

"""

def main(ref, obs):
	
	r = readFASTA(open(ref, "rU"))
	o = readFASTA(open(obs, "rU"))

	for x in r:
		for y in o:
			for i in range(len(r[x])):
				#compare sequences
				if r[x][i] != o[y][i]:
					print "\t".join([y, str(i+1), str(i+2), r[x][i], o[y][i]])	
	
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
#			val += l.strip()
			f_out[k] = val
			return f_out
						
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('--ref', action = 'store', type = 'string', dest = 'reference', help = 'Reference FASTA file')
	parser.add_option('--obs', action = 'store', type = 'string', dest = 'observed', help = 'Observed FASTA (compared to reference) file')
	(option, args) = parser.parse_args()
	
	main(option.reference, option.observed)


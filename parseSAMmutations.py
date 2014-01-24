"""
parseSAMmutations.py

Takes input of SAM file of alignments; outputs file containing, for each read: 
0) read ID
1) read length
2) number of mismatches
3) number of indels
4) mutation locations
5) mutation identities

Matt Rich, 3/2013

"""

def main(sam, concordant):
	#for each line in SAM file
	for line in open(sam, 'r'):
		#skip header lines
		if line.startswith('@'):
			continue
		if concordant == True and 'YT:Z:CP' not in line:
			continue
		#otherwise, parse information from line
		l = line.strip().split('\t')
		aln_data = l[11:]
		read_id = l[0]
		read_length = getPosIden(l[5], l[9])
		read_start = l[3]
		try:
			n_mis = aln_data[[ aln_data[i].startswith('XM:') for i in range(len(aln_data)) ].index(True)].split(':')[2]
			n_ind = aln_data[[ aln_data[i].startswith('XG:') for i in range(len(aln_data)) ].index(True)].split(':')[2]
			print '\t'.join([read_id, str(read_length), n_mis, n_ind])	
		except ValueError:
			pass
	
		
def getPosIden(cigar, seq):
	#trim soft-clipped bases off seq
	try:
		seq = seq[int(re.match('^(\d+)S', cigar).groups()[0]):]
	except AttributeError:
		pass
	try:
		seq = seq[0:-1*int(re.match('.*[IDM](\d+)S$', cigar).groups()[0])]
	except AttributeError:
		pass
	return len(seq)
	
		


if __name__ == "__main__":
	from optparse import OptionParser	
	import re
	parser = OptionParser()
	parser.add_option('-f', action = 'store', type = 'string', dest = 'inputfile', help = 'SAM-format file containing alignments')
	parser.add_option('--require-concordant', action = 'store_true', dest = 'concordant', help = 'only output data from concordant read pairs?', default=False )	
	(option, args) = parser.parse_args()

	main(option.inputfile, option.concordant)		
	
	
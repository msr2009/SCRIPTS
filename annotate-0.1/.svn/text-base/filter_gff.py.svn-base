"""
filter_gff.py

filters gff3-formatted annotation files to not include lines whose descriptions don't begin with 'ID'
and specific annotations that are not wanted. For those annotations not beginning with 'ID' that are
wanted to be kept, changes to 'ID'

"""

def main(gff3):
	#set of annotations to exclude
	exclude_set = set(['chromosome', 'CDS', 'repeat_region', 'region', 'nucleotide_match', 'long_terminal_repeat', 'noncoding_exon', 'gene', 'pseudogene', 'intron'])
	#open output file
	f_out = open(gff3 + '.filtered', 'w')
	#loop through gff3 file
	for line in open(gff3, 'r'):
		if line.startswith('#') != True:
			l = line.strip().split('\t')
			#check that it's an annotation we want to keep
			try:
				if l[2] not in exclude_set:
					#if it is, check if the description starts with "ID"
					if l[8].startswith('ID') != True:
						d = list(l[8])
						x = d.index('=')
						d = ['I','D'] + d[x:]
						print >> f_out, '\t'.join(l[0:-1]) + '\t' + ''.join(d)
					else:
						print >> f_out, line.strip()
			except IndexError:
				return 0		
		else:
			print >> f_out, line.strip()
	f_out.close()	
			
if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'inputfile', help = 'gff3 file to be filtered')
	(option, args) = parser.parse_args()

	main(option.inputfile)
"""
read_fuser_overlaps.py

Find read_fuser overlaps, so that you don't have to!

Matt Rich, 4/2013
"""

def main(read1, read2, rc):
	from Bio import pairwise2
	import re
		
	#perform alignment with high gap penalties
	aln = ''
	if rc == False:	
		aln = pairwise2.align.localms(read1, complement(read2)[::-1], 2, -1, -5, -5)[0]
	else:
		aln = pairwise2.align.localms(read1, read2, 2, -1, -5, -5)[0]

	#forward read overlap start and reverse overlap end are given in alignment output
	r1_overlap_start = aln[3]
	r2_overlap_end = aln[4]
	
	#evaluate the number of dashes flanking each read in the alignment 
	#use an re
	dashes = re.compile("^(-*)[AaCcTtGgNn]+(-*)$")
	r1_groups = re.match(dashes, aln[0]).groups()
	r2_groups = re.match(dashes, aln[1]).groups()
	
	r1_overlap_start = 0 + len(r2_groups[0])
	r1_overlap_end = len(read1) - len(r2_groups[1]) - 1
	
	r2_overlap_start = 0 + len(r1_groups[0])
	r2_overlap_end = len(read2) - len(r1_groups[1])	- 1
	
	return [r1_overlap_start, r1_overlap_end, r2_overlap_start, r2_overlap_end]
	
def complement(base):
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'a':'T', 't':'A', 'c':'G', 'g':'C', 'n':'N'}
	base = list(base)
	return ''.join( [ comp[b] for b in base ] )
	
if __name__ == "__main__":
	
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-1', '--read1', action = 'store', type = 'string', dest = 'read1', help = 'Read 1')
	parser.add_option('-2', '--read2', action = 'store', type = 'string', dest = 'read2', help = 'Read 2')
	parser.add_option('--rc', action = 'store_true', dest = 'rc', help = 'Read2 already reverse-complemented (default = False)', default=False)	

	(option, args) = parser.parse_args()

	overlaps = main(option.read1, option.read2, option.rc)
	print "Read 1 overlap start: ", overlaps[0]
	print "Read 1 overlap end: ", overlaps[1]
	print "Read 2 overlap start: ", overlaps[2]
	print "Read 2 overlap end: ", overlaps[3]

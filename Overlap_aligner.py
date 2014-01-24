# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def main(r1, r2):
	from Bio import pairwise2
	from Bio.Seq import Seq
	
	# <codecell>
	
	forward_wt_read = Seq(r1)
	reverse_wt_read = Seq(r2).reverse_complement()
	
	
	alignment = pairwise2.align.globalms(forward_wt_read, reverse_wt_read, 2, -1, -3, -1)
	
	# <codecell>
	
	#print alignment[0][2:]
	#print alignment[0][0]
	#print alignment[0][1]
	
	# <codecell>
	
	forward_wt_aligned = alignment[0][0]
	reverse_wt_aligned = alignment[0][1]
		
	# <codecell>
	
	overlap = [i for i, x in enumerate(forward_wt_aligned) if reverse_wt_aligned[i] != "-" and forward_wt_aligned[i] != "-"]
	if forward_wt_aligned[0] == "-":
	    read1_overlap_start = 0
	    read1_overlap_end = len(overlap) - 1
	    read2_overlap_start = overlap[0]
	    read2_overlap_end = overlap[-1]
	else:
	    read1_overlap_start = overlap[0]
	    read1_overlap_end = overlap[-1]
	    read2_overlap_start = 0
	    read2_overlap_end = len(overlap) - 1
	
	# <codecell>
	
	return [read1_overlap_start, read1_overlap_end, read2_overlap_start, read2_overlap_end]
	
	
if __name__ == "__main__":
	
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-1', '--read1', action = 'store', type = 'string', dest = 'r1', help = 'Read 1')
	parser.add_option('-2', '--read2', action = 'store', type = 'string', dest = 'r2', help = 'Read 2')
	
	(option, args) = parser.parse_args()

	print main(option.r1, option.r2)



"""
remove_indel_reads.py

removes indel reads from a fastq file, and corresponding reads from other fastqs

Matt Rich, 6/2013
"""

def main(forward, reverse, index, start, stop, target, seq):
	#create outfiles for all indices (and one for headers that don't match any indices close enough)
	out_f = open(forward.strip(".fq")+".nodel.fq", 'w')
	out_r = open(reverse.strip(".fq")+".nodel.fq", 'w')
	out_i = open(index.strip(".fq")+".nodel.fq", 'w')
		
	#open fastq files		
	forw = open(forward, 'r')
	rev = open(reverse, 'r')
	ind = open(index, 'r')
		
	#start looping through files, reading in 4 lines at a time
	while True:
		
		i_read = [ ind.readline().strip() for i in range(4) ]
		if i_read[0] == '':
			break
		
		f_read = [ forw.readline().strip() for i in range(4) ]
		r_read = [ rev.readline().strip() for i in range(4) ]  
		
		#get sequence, truncate to necessary length
		s = ''
		if target == "forward":
			s = f_read[1][start:stop]
		elif target == "reverse":
			s = r_read[1][start:stop]
		else:
			s = i_read[1][start:stop]
 		
		#perform pairwise2 alignment to last 45 bases of the forward read
		aln = pairwise2.align.globalms(seq, s, 2, -1, -5, -5)
		#print aln[0]
		
		#if there isn't a gap, then print out the lines to output fastq files
		if "-" not in aln[0][1]:
			print >> out_f, '\n'.join(f_read)
			print >> out_r, '\n'.join(r_read)
			print >> out_i, '\n'.join(i_read)
		
	out_f.close()
	out_r.close()
	out_i.close()
			
if __name__ == '__main__':
	from optparse import OptionParser
	from Bio import pairwise2
	
	parser = OptionParser()
	parser.add_option('--forward', action = 'store', type = 'string', dest = 'forward', help = 'path to forward reads file')
	parser.add_option('--reverse', action = 'store', type = 'string', dest = 'reverse', help = 'path to reverse reads file')
	parser.add_option('--index', action = 'store', type = 'string', dest = 'index', help = 'path to index reads file')
	parser.add_option('--sequence', action = 'store', type = 'string', dest = 'seq', help = 'sequence to match for checking indels')
	parser.add_option('--target_read', action = 'store', type = 'string', dest = 'target', help = 'which read to check? (forward, index, reverse)')
	parser.add_option('--start', action = 'store', type = 'int', dest = 'start', default = None, help = 'starting position for alignment')
	parser.add_option('--end', action = 'store', type = 'int', dest = 'stop', default = None, help = 'ending position for alignment' )
	(option, args) = parser.parse_args()
	
	main(option.forward, option.reverse, option.index, option.start, option.stop, option.target, option.seq)	
	
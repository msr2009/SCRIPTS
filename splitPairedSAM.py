"""
splitPairedSAM.py

splits paired-end SAM file (like the output of Bowtie2) into two fastq files or two pairs of fasta and qual files

Matt Rich, 03/2013
"""


def main(paired_sam, fasta):
	
	#open sam file
	ps = open(paired_sam, 'r')

	#open output file(s)
	outname = '.'.join(paired_sam.split('.')[0:-1])
	if fasta == True:
		out_f = open(outname+'_F.fa','w')
		out_r = open(outname+'_R.fa','w')
		qual_f = open(outname+'_F.fa.qual','w')
		qual_r = open(outname+'_R.fq.qual','w')
	else:
		outf = open(outname+'_F.fq','w')
		outr = open(outname+'_R.fq','w')

	#loop through SAM file, splitting out reads
	while True:
		read1 = ps.readline().strip().split('\t')
		read2 = ps.readline().strip().split('\t')
		
		#get through the header
		if read1[0][0] == '@' or read1[0][0] == '#':
			if read2[0]['@'] == True or read2[0][0] == '#':
				continue
			else:
				read1 = read2
				read2 = ps.readline().strip().split('\t')
					
		#stop condition
		if read1 == [''] or read2 == ['']:
			try:
				qual_f.close()
				qual_r.close()
			except:
				pass

			out_f.close()
			out_r.close()
			break 

		#check to make sure the reads correspond to two pairs of the same read
		if read1[0] != read2[0]:
			raise ValueError('Read names in consecutive lines are not consistent. Exiting...')
			return 0

		if fasta == True:
			print >> out_f, '>'+read1[0] + '\n' + read1[8] 
			print >> out_r, '>'+read2[0] + '\n' + read2[8]
			print >> qual_f, '>'+read1[0] + '\n' + read1[9]
			print >> qual_r, '>'+read2[0] + '\n' + read2[9]
		else:
			print >> out_f, '@' + read1[0] + '\n' + read1[8] + '\n+\n' + read1[9]
			print >> out_r, '@' + read2[0] + '\n' + read2[8] + '\n+\n' + read2[9]

if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
        parser.add_option('--sam', action = 'store', type = 'string', dest = 'sam', help = 'SAM file to be split')
	parser.add_option('--fasta', action = 'store_true', dest = 'fasta', help = 'Output FASTA and QUAL files instead of FASTQ?', default = False)
	(option, args) = parser.parse_args()

	main(option.sam, option.fasta)

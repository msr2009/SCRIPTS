"""
splitFASTQforCross_match.py

Helper script that takes FASTQ file as input, output two files: one FASTA file with sequence, one QUAL file with phred-scaled quality scores

Matt Rich, 03/2013
"""

def main(fastq, phred, divreads, include):
	
	p = 64-phred
	
	#a few counters
	nfile = ''
	if divreads != float('inf'):
		nfile = 1
	
	nreads = 1
	
	#create set of unique ReadIDs
	uniqueIDs = set()
	if include != False:
		uniqueIDs = parseSAM(include)	
						
	#open fastq
	f = open(fastq, 'r')
	
	#create output files
	outfile = '.'.join(fastq.split('.')[0:-1])
	#make first outfiles
	if nfile == '':
		fasta_out = open(outfile+'.fa', 'w')
		qual_out = open(outfile+'.fa.qual', 'w')	
	else:
		fasta_out = open(outfile+'_'+str(nfile)+'.fa', 'w')
		qual_out = open(outfile+'_'+str(nfile)+'.fa.qual', 'w')	
	
	#read file, print out lines to files
	while True:
		
		if nreads % divreads == 0:
			#close the current files
			fasta_out.close()
			qual_out.close()
			
			#open the new files
			fasta_out = open(outfile+'_'+str(nfile)+'.fa', 'w')
			qual_out = open(outfile+'_'+str(nfile)+'.fa.qual', 'w')
			nfile += 1
		
		read = [ f.readline().strip() for i in range(4) ]
		
		if read == ['','','','']:
			break
		
		#only print reads that are unique, if we care
		if include != False and read[0][1:-2] in uniqueIDs:
		
			#print FASTA
			print >> fasta_out, ">"+read[0]
			print >> fasta_out, read[1]
			
			#print QUAL
			print >> qual_out, ">"+read[0]
			print >> qual_out, " ".join([ str(ord(x) - p) for x in list(read[3]) ])
		
			nreads += 1
			
		#otherwise, just print out all the reads
		elif include == False:
			#print FASTA
			print >> fasta_out, ">"+read[0]
			print >> fasta_out, read[1]
			
			#print QUAL
			print >> qual_out, ">"+read[0]
			print >> qual_out, " ".join([ str(ord(x) - p) for x in list(read[3]) ])
		
			nreads += 1
		
	#close the current (last) file
	fasta_out.close()
	qual_out.close()

def parseSAM(f):
	goodreads = set()

	for line in open(f,'r'):
		l = line.strip().split('\t')
		if l[0][0] == "@":
			goodreads.add(l[0][1:])
		else:
			goodreads.add(l[0])
	
	return goodreads
		
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--fastq', action = 'store', type = 'string', dest = 'fastq', help = 'FASTQ file to be split')
	parser.add_option('--phred', action = 'store', type = 'int', dest = 'phred', default = 33, help = 'phred correction (default=33)')
	parser.add_option('--divide', action = 'store', type = 'float', dest = 'divreads', default=float("inf"), help = 'should this file be divided into multiple smaller files? if so, how many reads per file?')
	parser.add_option('--include-reads', action = 'store', type = 'string', dest = 'include', default = False, help = 'File containing unique (non-duplicate) ReadIDs in the first column (can be SAM). Only these reads will be processed.')
	(option, args) = parser.parse_args()

	main(option.fastq, option.phred, option.divreads, option.include)

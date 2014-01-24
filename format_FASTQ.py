#/usr/bin/python

"""
format_FASTQ.py

Converts all non-ACTG characters in the sequence line (line 2) of a FASTQ file to N's. Other characters can be ignored by setting the -i option.

python format_FASTQ.py -f <FASTQFILE> -i <OTHERCHARACTERSTOBEIGNORED>

Matt Rich, Mar.2012
"""

def main(inputfile, ignorelist=''):
	#open new output file
	f_out = open(inputfile.split('.fq')[0]+'_N.fq', 'w')
	#counter for lines
	count = 0
	#add contents of ignorelist to ACTG
	ignore = 'ACTG' + ignorelist
	#read lines in from file to be converted
	for line in open(inputfile,'r'):
		#if it's the sequence line
		if count % 4 == 1:
			#change all non-ACTG's to N's
			print >> f_out, re.sub("[^%s]" % ignore, 'N', line.strip())
		#otherwise, print line as-is
		else:
			print >> f_out, line.strip()
		#increment counter
		count += 1
	
	#close file
	f_out.close()
		
if __name__ == "__main__":
	from optparse import OptionParser
	import sys, re
	
	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'inputfile', help = 'fasta file to be converted')
	parser.add_option('-i', '--ignore', action = 'store', type = 'string', dest = 'ignorelist', help = 'characters to be ignored (other than ACTG)')
	(option, args) = parser.parse_args()
	
	if option.ignorelist != None:
		main(option.inputfile, option.ignorelist)
	else:
		main(option.inputfile)
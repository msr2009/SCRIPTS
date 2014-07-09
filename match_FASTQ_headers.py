"""
match_fastq_headers.py

removes reads not found in one fastq file from all others

Matt Rich, 5/2014

"""

def main(fq_key, otherfqs):
	
	#first, read in the headers in fq_key
	k = set()
	for record in read_fastq(fq_key):
		k.add(record[0][0:-6])

	#then, read all the other fastqs in, printing only 
	#reads found in the key fq

	for f in otherfqs:	
		print f
		f_out = open(f+".match", "w")
		for record in read_fastq(f):
			if record[0][0:-6] in k:
				print_fastq(record, f_out)
		f_out.close()
			

if __name__ == "__main__":
	from fastq_tools import read_fastq, print_fastq
        from optparse import OptionParser

        parser = OptionParser()
        parser.add_option('--key', action = 'store', type = 'string', dest = 'key', help = 'one FASTQ to rule them all')
        parser.add_option('--fq', action = 'store', type = 'string', dest = 'fq', help = 'comma-delimited list containing all other fastqs')
        (option, args) = parser.parse_args()

        main(option.key, option.fq.split(","))	

"""
simulate_illumina_reads.py

takes input of fragment being sequenced, outputs fastq file containing reads

Matt Rich 5/2014

"""

def main(seq, flen, rlen, out, err_rate, avg_qual, num_reads):
		
	#write first read
	f_out = open(out + "_1.fq", "w")
	
	for i in range(num_reads):	
		r, q = simRead(seq, flen, err_rate, avg_qual)
		print >> f_out, "@simRead_" + str(i) + "_1"
		print >> f_out, r
		print >> f_out, "+"
		print >> f_out, q
	f_out.close()
	
	#write second read
	r_out = open(out + "_2.fq", "w")
	
        for i in range(num_reads):      
                r, q = simRead(complement(seq)[::-1], rlen, err_rate, avg_qual)
                print >> r_out, "@simRead_" + str(i) + "_2"
                print >> r_out, r
                print >> r_out, "+"
                print >> r_out, q
        r_out.close()


def complement(base):
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'a':'T', 't':'A', 'c':'G', 'g':'C', 'n':'N'}
	base = list(base)
	return ''.join( [ comp[b] for b in base ] )


def simRead(seq, length, err_rate, avg_qual):
	read = []
	qual = []
	seq = seq.upper()
	nucs = set("ACTG")
 	for l in range(length+1):
		read.append(seq[l])
		qual.append(chr(avg_qual+33))	
	return ( "".join(read), "".join(qual) )	


if __name__ == "__main__":
	
	from optparse import OptionParser
	from random import random, sample
	from copy import copy

	parser = OptionParser()
	
	parser.add_option("--seq", action = "store", type = "string", dest = "sequence", help = "sequence to be sequenced")
	parser.add_option("--flen", action = "store", type = "int", dest = "flen", help = "length of forward read")
	parser.add_option("--rlen", action = "store", type = "int", dest = "rlen", help = "length of reverse read")
	parser.add_option("--output", action = "store", type = "string", dest = "output", help = "prefix of output FASTQs")	

	(option, args) = parser.parse_args()

	main(option.sequence, option.flen, option.rlen, option.output, .01, 25, 1000)




"""
barcode_sequence_map-phrap.py

from fastq containing sequences marked with a barcode in readID,
creates a consensus sequence for each barcode using phrap

Matt Rich, 5/2014
"""

def main(fq, out, splitsize, tmpdir):
	bc_dict = {}
	
	for record in read_fastq(fq):
		#strip barcode out of ReadID
		bc = record[0].split("%%")[1]
		if bc in bc_dict:
			bc_dict[bc].append( (record[1], record[2], record[0]) )
		else:
			bc_dict[bc] = [ (record[1], record[2], record[0]) ]

	#if we're splitting...
	if splitsize != None:
		k = bc_dict.keys()
		#for each chunk of barcodes		
		for chunk in getChunks(k, splitsize):
			#print the barcodes as a FASTQ to a temp dir
			tmpf = open(tmpdir+chunk[0]+"_"+chunk[-1]+".fq", "w")
			for b in chunk:
				 [ print_fastq([x[2], x[0], x[1]], file=tmpf) for x in bc_dict[b] ]
			tmpf.close()
		return 1

	#for each barcode...	
	outfile = open(out, "w")
	for bc in bc_dict:
		#if barcode only has one read, check for correct length, Ns
		if len(bc_dict[bc]) == 1:
			if  "N" in bc_dict[bc][0][0]:
				print >> outfile, "\t".join([bc, "NA", "N"])
			else:
				print >> outfile, "\t".join([bc, bc_dict[bc][0][0], "*"])

		#if there are more than one read, send reads to phrap for assembly
		else:	
			#if the reads are all the same, we can just call the consensus
			if len(set(zip(*bc_dict[bc])[0])) == 1:
				print >> outfile, "\t".join([bc, bc_dict[bc][0][0], "*"])
				continue

			#print fasta and qual files
			fa_out = open(tmpdir+bc+".fa", "w")
			qual_out = open(tmpdir+bc+".fa.qual", "w")
			
			#sample 1000 reads if there are more than 1000
			sampled_reads = bc_dict[bc]

			if len(bc_dict[bc]) > 500:
				sampled_reads = random.sample(bc_dict[bc], 500)
				
			for read in sampled_reads:
				print >> fa_out, ">" + read[2]
				print >> fa_out, read[0]
				print >> qual_out, ">" + read[2]
				print >> qual_out, " ".join([ str(ord(y)-33) for y in read[1] ])
			
			fa_out.close()
			qual_out.close()
		
			call = "phrap " + tmpdir+bc +".fa -minmatch 20 -node_seg 16 &> phrap.out"
			phrap_out = subprocess.call(call, shell=True)

			#check how many contigs in phrap output
			if int(subprocess.check_output("grep -c '>' " + tmpdir+bc + ".fa.contigs", shell=True)) > 1:
				print >> outfile, "\t".join([bc, "NA", "con"])
			else:	
				#open file, read contig
				contig = ""
				for line in open(tmpdir+bc + ".fa.contigs", "r").readlines()[1:]:
					contig += line.strip()
				print >> outfile, "\t".join([bc, contig, "*"])
		
			#delete all bc phrap files
			subprocess.call("rm " + tmpdir+bc + ".fa*", shell=True)	
		
	outfile.close()


def getChunks(l, n):
	for i in xrange(0, len(l), n):
		yield l[i:i+n]


if __name__ == "__main__":
	
	from fastq_tools import read_fastq, print_fastq
	from optparse import OptionParser
	import subprocess, random

	parser = OptionParser()
	parser.add_option('--fq', action = 'store', type = 'string', dest = 'fq', help = "FASTQ file containing barcode-id'd sequences")
	parser.add_option('-o', action = "store", type = "string", dest = "outfile", help = "output file")
	parser.add_option('--split', action = "store", type = "int", dest = "splitsize", help = "if FQs should be split into barcode groups, how many BCs per FQ", default=None)
	parser.add_option('--tmp', action = "store", type = "string", dest = "tmpdir", help = "directory for storing split barcode fastqs", default=None)
	(option, args) = parser.parse_args()
	
	main(option.fq, option.outfile, option.splitsize, option.tmpdir)


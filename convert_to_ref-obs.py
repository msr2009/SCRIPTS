#/usr/bin/python

def main(bed, fasta_genome):
	#read fasta file of genome to compare against
	genome = []
	for record in SeqIO.parse( open(fasta_genome, 'r'), 'fasta' ):	
		genome.append(record.seq.tostring())
	
	#open output file
	f_out = open(bed+'.refobs', 'w')
		
	#read each snp from BED file, reformat
	for line in open(bed):
		l = line.strip().split('\t')
		#make list to store new line
		newline = l[0:4]
		#check which allele is reference
		c = chromosome_conversion(l[0][3:])
		if genome[c-1][int(l[1])] == l[4]:
			newline.append(l[4])
			newline.append(l[5])
		elif genome[c-1][int(l[1])] == l[5]:
			newline.append(l[5])
			newline.append(l[4])
		else:
			raise ValueError('Neither allele is reference!')
		newline += l[6:]
		print >> f_out, '\t'.join(newline)
	f_out.close()
		
def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':'Mito'}
	return chrom_conv[chrom_number]
	
if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	from Bio import SeqIO
	
	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'inputfile', help = 'fasta file of genome sequence')
	parser.add_option('-s', '--snps', action = 'store', type = 'string', dest = 'snps', help = 'BED file with snps')
	(option, args) = parser.parse_args()

	main(option.snps, option.inputfile)
"""
mutation_caller.py

takes input of two genomes - reference and experimental - and outputs list of SNP locations and identities

"""

def main(ref_file, exp_file):
	mutations = {}
	out_file = open('.'.join(exp_file.split('.')[0:-1])+'.mutations', 'w')
	#print header for output
	print >> out_file, '\t'.join(['chr', 'start', 'stop', 'ref', 'obs'])
	
	#open genome files, create 
	genome1 = []
	for record in SeqIO.parse(ref_file, 'fasta'):
		genome1.append(record.seq)
	genome2 = []
	for record in SeqIO.parse(exp_file, 'fasta'):
		genome2.append(record.seq)

	#loop through genomes, finding mutations
	for chrom in range(len(genome1)):
		#check that chromosomes are of equal length
		if len(genome1[chrom]) != len(genome2[chrom]):
			sys.exit('Chromosomes not of equal lengths - check order.')
		
		#find mutations
		if str(genome1[chrom]) != str(genome2[chrom]):
			m = {}
			# loop through chromosome, find mutations
			for pos in range(len(genome1[chrom])):
				if genome1[chrom][pos] != genome2[chrom][pos] and genome2[chrom][pos] != 'N':
					print >> out_file, '\t'.join([chromosome_conversion(chrom+1), str(pos+1), str(pos+1), genome1[chrom][pos], genome2[chrom][pos]])
	out_file.close()

def chromosome_conversion(chrom_number):
	chrom_conv = {1:'1', 2:'2', 3:'3', 4:'4', 5:'5', 6:'6', 7:'7', 8:'8', 9:'9', 10:'10', 
					11:'11', 12:'12', 13:'13', 14:'14', 15:'15', 16:'16', 17:'Mito'}
	return chrom_conv[chrom_number]

if __name__ == '__main__':
	from optparse import OptionParser
	from Bio import SeqIO
	import sys
	
	parser = OptionParser()
	parser.add_option('-r', '--reference', action = 'store', type = 'string', dest = 'reference', help = 'reference genome sequence (FASTA)')
	parser.add_option('-o', '--observed', action = 'store', type = 'string', dest = 'observed', help = 'observed genome sequence (FASTA)')
	(option, args) = parser.parse_args()

	main(option.reference, option.observed)
			
	








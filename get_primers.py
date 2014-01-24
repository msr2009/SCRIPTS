


def main(target_genes, features_file, genome_file, offset, flank=',', Tm=None):

	"""GET DATA"""
		
	#populate list of chromosomes in genome	
	genome = {}
	for record in SeqIO.parse( open(genome_file, 'r'), 'fasta' ):
		genome[chromosome_conversion(record.description.split(' [')[4].split('=')[1][0:-1])] = record.seq.tostring()
	
	#populate dictionary of genes and start/stop, chromosome
	#noncoding[ID] = [ chrom, start, stop, strand, name ]
	features = {}
	for line in open(features_file, 'r'):
		if line.startswith('#') != True:
			l = line.strip().split('\t')
			if l[2] == "gene":
				if l[8].split(';')[2].startswith('gene'):
					features[l[8].split(';')[0].split('=')[1]] = [ chromosome_conversion(l[0][3:]), int(l[3]), int(l[4]), l[6], l[8].split(';')[2].split('=')[1] ]
				else:
					features[l[8].split(';')[0].split('=')[1]] = [ chromosome_conversion(l[0][3:]), int(l[3]), int(l[4]), l[6], l[8].split(';')[0].split('=')[1] ]

	
	"""GET SEQUENCES"""
	
	#if using Tm calculation only
	if Tm != None:
		print '\t'.join(['ORF', 'Name', 'PrimerF', 'PrimerR', 'TM-F', 'TM-R'])
		for gene in open(target_genes, 'r'):
			seq = ''
			gene = gene.strip()
			#get info for gene
			[c, start, stop, strand, name] = features[gene][0:5]
			#retrieve sequence
			seq = genome[c][start-1-offset:stop+offset]
			if strand == '-':
				seq = complement(seq[::-1])
			
			#find 5' primer
			primer5 = ''
			for i in range(len(seq)):
				primer5 += seq[i]
				if len(primer5) >= 14:
					#check Tm of primer
					if calculate_Tm(primer5) > Tm:
						break
			
			#find 3' primer
			primer3 = ''
			#remove stop codon
			seq = complement(seq[0:-3][::-1])
			for i in range(len(seq)):
				primer3 += seq[i]
				if len(primer3) >= 14:
					#check Tm of primer
					if calculate_Tm(primer3) > Tm:
						break
			
			flankseq = flank.split(',')
			print '\t'.join([gene, name, primer5, primer3+flankseq[1], str(calculate_Tm(primer5)), str(calculate_Tm(primer3))])
	
	
	#otherwise, use Primer3
	else:
		#output file of primer3 input data
		p3 = open('primer3_input', 'w')
		
		#find sequences for genes in target list using genome sequence
		for gene in open(target_genes, 'r'):
			seq = ''
			gene = gene.strip()
			#get info for gene
			[c, start, stop, strand] = features[gene][0:4]
			#retrieve sequence
			seq = genome[c][start-1-offset:stop+offset]
			if strand == '-':
				seq = complement(seq[::-1])
			
			#print output for primer3 input
			print >> p3, 'SEQUENCE_ID='+gene
			print >> p3, 'SEQUENCE_TEMPLATE='+seq
			print >> p3, 'SEQUENCE_TARGET='+str(offset)+','+str(abs(start-stop)+1)
			print >> p3, 'PRIMER_TASK=generic'
			print >> p3, 'PRIMER_PICK_LEFT_PRIMER=1'
			print >> p3, 'PRIMER_PICK_INTERNAL_OLIGO=0'
			print >> p3, 'PRIMER_PICK_RIGHT_PRIMER=1'
			print >> p3, 'PRIMER_PRODUCT_SIZE_RANGE='+str(len(seq)-offset*2)+'-'+str(len(seq)-offset-1)
			print >> p3, 'PRIMER_PRODUCT_OPT_SIZE='+str(len(seq)-offset*2)
			print >> p3, 'PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0'
			print >> p3, 'PRIMER_THERMODYNAMIC_ALIGNMENT=0'
			print >> p3, '='
		
		p3.close()	
		cmd = 'primer3_core -output=primer3_output < primer3_input'
		os.system(cmd)
		
		#read in primer3_output, split by '=\n' to divide lines
		f = open('primer3_output','r')
		s = f.read().split('=\n')
		
		#separate flanking sequences
		flankseq = flank.split(',')
		
		#print header
		print '\t'.join(['ORF', 'PrimerF', 'PrimerR', 'TM-F', 'TM-R', 'Length'])
		
		#get data from each primer, print sequences
		for i in range(len(s)-1):
			[orf, seqF, seqR, tmF, tmR, length] = ['','','','','','']
			p = s[i].split('\n')
			for j in p:
				if j.startswith('SEQUENCE_ID='):
					orf = j.split('=')[1]
				if j.startswith('PRIMER_LEFT_0_SEQUENCE='):
					seqL = flankseq[0].lower() + j.split('=')[1]
				if j.startswith('PRIMER_RIGHT_0_SEQUENCE='):
					seqR = flankseq[1].lower() + j.split('=')[1]
				if j.startswith('PRIMER_LEFT_0_TM='):
					tmF = j.split('=')[1]
				if j.startswith('PRIMER_RIGHT_0_TM='):
					tmR = j.split('=')[1]
				if j.startswith('PRIMER_PAIR_0_PRODUCT_SIZE='):
					length = j.split('=')[1]
			
			if tmF == '':
				print '\t'.join([orf, 'NO PRIMERS FOUND'])
			else:		
				print '\t'.join([orf, seqL, seqR, tmF, tmR, length])
			

def calculate_Tm(s):
	sequence = s.upper()
	"""64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)"""
	return 64.9 + 41*(sequence.count('G') + sequence.count('C') - 16.4)/(len(sequence))	
			
def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':17, 'mitochondrion':17}
	if chrom_number.startswith('chr'):
		chrom_number = chrom_number[3:]
	
	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:			
		return chrom_conv[chrom_number]

def complement(base):
	compbase = []
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', '[':']', ']':'['}
	for i in range(len(base)):
		compbase.append(comp[base[i].upper()])
	return ''.join(compbase)
	
	
if __name__ == "__main__":
	from optparse import OptionParser
	import sys, os
	from Bio import SeqIO
	
	parser = OptionParser()
	parser.add_option('-l', '--target_list', action = 'store', type = 'string', dest = 'inputfile', help = 'file containing genes that need primers')
	parser.add_option('-f', '--features', action = 'store', type = 'string', dest = 'features', help = 'gff file containing genomic features')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'genome', help = 'fasta file containing genome sequence')
	parser.add_option('-o', '--offset', action = 'store', type = 'int', dest = 'offset', help = 'number of bases upstream and downstream of gene to include for priming sites')
	parser.add_option('--flank-5', action = 'store', type = 'string', dest = 'flank5', 
					help = "comma-delimited, 2-member list of flanking sequences to be added to 5'-end of forward AND reverse primers. E.g., 'ACGAGCGAGAG,GGACGTTTAGC'")
	parser.add_option('--tm', action = 'store', type = 'int', dest = 'melt', help = "optimal Tm. This option triggers a second workflow, in which primers are found by extending \ 
					from either the 5' or 3' end of the sequence, until the specified Tm is reached.")
	(option, args) = parser.parse_args()
	
	melt = None
	if option.melt != None:
		melt = option.melt

	if option.flank5 != None:
		main(option.inputfile, option.features, option.genome, option.offset, option.flank5, melt)
	else:
		main(option.inputfile, option.features, option.genome, option.offset, ',', melt)
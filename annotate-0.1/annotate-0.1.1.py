"""
yeast_annotation.py

takes input of BED formatted file with mutation locations (from mutation_caller.py) and annotates SNPs

python yeast_annotation.py -f <BED file containing mutations> -s <ORF sequences> -n <non-coding GFF file for annotations> -g <Genome sequence as FASTA>

Matt Rich, 12/19/2012

"""

def main(BED, orfs, noncoding_file, genome_file, upstream):

	"""GET DATA TOGETHER"""

	#populate list of chromosomes in genome	
	genome = {}
	temp_genome = readFASTA( open(genome_file, 'r') )
	for record in temp_genome:
		genome[chromosome_conversion(record.split(' [')[4].split('=')[1][0:-1])] = temp_genome[record]


	#gather data about genes:
	#orf id, start, stop (both 5'-3'), exons, introns, protein length, chr
	genes = {}
	temp_genes = readFASTA( open(orfs, 'r') )
	
	for record in temp_genes:
		start = 0
		stop = 0
		exons = map(lambda x: x.split('-'), record.split(', ')[1].split(' ')[-1].split(','))
		ch = chromosome_conversion(record.split('Chr ')[1].split(' ')[0])
		gene_id = record.split()[0]
		
		
		#find introns in genes that have them
		introns = []
		if gene_id.split('-')[0][-1] == 'C':
			exons = exons[::-1]
			start = int(exons[0][0])
			stop = int(exons[-1][1])
			if len(exons) > 1:
				for e in xrange(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		else:
			start = int(exons[0][0])
			stop = int(exons[-1][1])
			if len(exons) > 1:
				for e in xrange(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		
		genes[gene_id] = [gene_id] + [start, stop] + [exons] + [introns] + [len(temp_genes[record])/3.0] + [ch]
			
	#populate second dictionary of non-coding annotations
	#noncoding[ID] = [ ID, chrom, regiontype, start, stop ]
	noncoding = {}
	randints = range(1000000)
	for line in open(noncoding_file, 'r'):
		if line.startswith('#') != True:
			l = line.strip().split('\t')
			noncoding[randints.pop()] = [ l[8].split(';')[0].split('=')[1], chromosome_conversion(l[0][3:]), l[2], int(l[3]), int(l[4]) ]
		if line.startswith('##FASTA'):
			break

	"""ANNOTATE SNPS"""
	
	#open output file
	f_out = open(BED+'.annotated', 'w')
	
	#start reading BED file
	f_in = open(BED, 'r')

	#read header
	header = f_in.next()
	print >> f_out, header.strip() + '\t' + '\t'.join(['annotation', 'region', 'protein'])

	
	#chr, start, stop, ref, obs
	while True:
		try:
			indel = False
			line = f_in.next()
			l = line.strip().split('\t')
			
			snp_start = int(l[1])-1
			snp_end = int(l[2])-1
			snp_len = len(l[4])-len(l[3])
			if snp_len != 0: 
				indel = True
			
			chrom_conv = chromosome_conversion(l[0])
			
			#make copy of chromosome with mutation
			mut_chr = list(copy(genome[chrom_conv]))
			
			#print list(genome[chromosome_conversion(l[0])])[snp_start:snp_end], l[3]	#this line will print the reference base from the BED file and genome, to check for concordance.
			mut_chr[snp_start:snp_end] = l[4]
					
			annotation = False
			
			#loop through genes, trying to find one containing snp
			for g in genes:
				#if gene on correct chromosome
				if genes[g][6] == chrom_conv:
					#CRICK STRAND
					if genes[g][0].split('-')[0][-1] == 'C':
						#found gene containing snp
						if snp_start <= genes[g][1] and snp_start >= genes[g][2]:
							#check if snp is in an intron:
							for intron in genes[g][4]:
								if snp_start <= intron[0] and snp_start >= intron[1]:
									#found an intronic snp - check if it's near a splice site
									if snp_start in xrange(intron[0]-2, intron[0]+1) or snp_start in xrange(intron[1],intron[1]+3):
										print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
									#otherwise it's just intronic (boring!)
									else:
										print >> f_out, '\t'.join(l + ['intron', genes[g][0], 'NA'])
							
							#if the snp is within gene start-stop, but not in an intron, then it's in an exon
							#remove introns
							if len(genes[g][4]) > 0:
								#I'm not entirely certain how to implement this indel calling on genes with introns (since it necessitates adding/subtracting bases)
								#so just call indels like SNPs for those genes:
								
								#remove intronic sequences from gene
								mut_gene = []
								wt_gene = []
								
								for exon in genes[g][3]:
									#check first if the snp is near a splice site
									if snp_pos in xrange(int(exon[0])-2,int(exon[0])+1) or snp_pos in xrange(int(exon[1]),int(exon[1])+3):
										#if it's not in the start/stop regions (these aren't splice-sites)
										if snp_pos not in xrange(genes[g][1]-2,genes[g][1]+1) and snp_pos not in xrange(genes[g][2],genes[g][2]+3):
											print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
										
									mut_gene += list(complement(mut_chr[ int(exon[1])-1:int(exon[0]) ])[::-1])
									wt_gene += list(complement(genome[genes[g][6]][ int(exon[1])-1:int(exon[0]) ])[::-1])
															
								#loop through codons, find mismatch
								if len(l[3]) != len(l[4]):
									print >> f_out, '\t'.join(l + ['indel_in_gene_containing_introns'])
								
								for codon in xrange(0, len(mut_gene), 3):
									if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
										mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
										wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
										if mut_aa != wt_aa:
											if mut_aa == '*':
												print >> f_out, '\t'.join(l + ['nonsense', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]) 
											else:	
												print >> f_out, '\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
										else:
											print >> f_out, '\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
							
								annotation = True

							#otherwise, it's not a gene with introns, so call indels correctly	
							else:	
								#extract total sequence from start:stop
								mut_gene = list(complement(mut_chr[ int(genes[g][2])-1:int(genes[g][1])+snp_len ])[::-1]) 
								wt_gene = list(complement(genome[genes[g][6]][int(genes[g][2])-1:int(genes[g][1])])[::-1]) 
															
								for exon in genes[g][3]:
									#check first if the snp is near a splice site
									if snp_start in xrange(int(exon[0])-2,int(exon[0])+1) or snp_start in xrange(int(exon[1]),int(exon[1])+3):
										#if it's not in the start/stop regions (these aren't splice-sites)
										if snp_start not in xrange(genes[g][1]-2,genes[g][1]+1) and snp_start not in xrange(genes[g][2],genes[g][2]+3):
											print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
																								
								#loop through codons, find mismatch
								for codon in xrange(0, len(mut_gene), 3):
									if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
										mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
										wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
										#if it's a snp, print out as snp
										if indel == False:
											if mut_aa != wt_aa:
												if mut_aa == '*':
													print >> f_out, '\t'.join(l + ['nonsense', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]) 
												else:	
													print >> f_out, '\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
											else:
												print >> f_out, '\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
										
										#otherwise, print the first mismatch as an indel, then find the next stop
										else:
											#these would be in-frame
											if snp_len % 3 == 0:
												if snp_len > 0:
													print >> f_out, '\t'.join(l + ['in-frame-insertion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
											 	elif snp_len < 0:
													print >> f_out, '\t'.join(l + ['in-frame-deletion', genes[g][0], wt_aa + str(codon/3+1) + 'X'])
											else:
												#these are out of frame
												#translate sequence, find first stop
												trunc = translate_sequence(''.join(mut_gene))
												if '*' in trunc:
													if snp_len > 0:
														print >> f_out, '\t'.join(l + ['frameshift-insertion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa + "; Truncation:" + str(trunc.index('*')+1) + '/' + str(int(genes[g][5])) ])
													elif snp_len < 0:
														print >> f_out, '\t'.join(l + ['frameshift-deletion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa + "; Truncation:" + str(trunc.index('*')+1) + '/' + str(int(genes[g][5])) ])
												else:
													if snp_len > 0:
														print >> f_out, '\t'.join(l + ['frameshift-insertion', genes[g][0], str(codon/3+1) + mut_aa + "; Full-length"])
													elif snp_len < 0:
														print >> f_out, '\t'.join(l + ['frameshift-deletion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa + "; Full-length"])
										break
																				
								annotation = True
							
						#if snp_pos isn't in a gene, check if it's upstream of a gene
						elif snp_start in xrange(genes[g][1]+1,genes[g][1]+upstream+1):
							print >> f_out, '\t'.join(l + ["5'-upstream", genes[g][0], str(genes[g][1]-snp_start)])
							annotation = True	
	
					#WATSON STRAND
					else:
						#found gene containing snp
						if snp_start >= genes[g][1] and snp_start <= genes[g][2]:
	
							#check if snp is in an intron:
							for intron in genes[g][4]:
								if snp_start >= intron[0] and snp_start <= intron[1]:
									#found an intronic snp - check if it's splice site
									if snp_start in xrange(intron[0],intron[0]+3) or snp_start in xrange(intron[1]-2,intron[1]+1):
										print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
									#otherwise, it's just intronic
									else:
										print >> f_out, '\t'.join(l + ['intron', genes[g][0], 'NA'])
									annotation = True	
							
							#if the snp is within the gene start-stop, but not in an intron, then it's in an exon

							
							#remove introns
							if len(genes[g][4]) > 0:
								#if the snp is within the gene start-stop, but not in an intron, then it's in an exon
								#remove intronic sequences
								mut_gene = []
								wt_gene = []
								
								for exon in genes[g][3]:
									#first, check if in a splice-site
									if snp_pos in xrange(int(exon[0]),int(exon[0])+3) or snp_pos in xrange(int(exon[1])-2,int(exon[1])+1):
										#if it's not the start/stop
										if snp_pos not in xrange(genes[g][1],genes[g][1]+3) and snp_pos not in xrange(genes[g][2]-2,genes[g][2]+1):
											print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
										
									mut_gene += list(mut_chr[ int(exon[0])-1:int(exon[1]) ])
									wt_gene += list(genome[genes[g][6]][int(exon[0])-1:int(exon[1])])
		
								#loop through codons, find mismatch								
								if len(l[3]) != len(l[4]):
									print >> f_out, '\t'.join(l + ['indel_in_gene_containing_introns'])

								for codon in xrange(0, len(mut_gene), 3):
									if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
										mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
										wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
										#check if synonymous, non-synonymous
										if mut_aa != wt_aa:
											if mut_aa == '*':
												print >> f_out, '\t'.join(l + ['nonsense', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]) 
											else:	
												print >> f_out, '\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
										else:
											print >> f_out, '\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
										annotation = True
		
								annotation = True

							#otherwise, this is a gene without introns, so call indels correctly
							else:
								#remove intronic sequences
								mut_gene = list(mut_chr[ int(genes[g][1])-1:int(genes[g][2])+snp_len ])
								wt_gene = list(genome[genes[g][6]][int(genes[g][1])-1:int(genes[g][2])])
								
								for exon in genes[g][3]:
									#first, check if in a splice-site
									if snp_start in xrange(int(exon[0]),int(exon[0])+3) or snp_start in xrange(int(exon[1])-2,int(exon[1])+1):
										#if it's not the start/stop
										if snp_start not in xrange(genes[g][1],genes[g][1]+3) and snp_start not in xrange(genes[g][2]-2,genes[g][2]+1):
											print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
										
									#mut_gene += list(mut_chr[ int(exon[0])-1:int(exon[1]) ])
									#wt_gene += list(genome[genes[g][6]][int(exon[0])-1:int(exon[1])])
		
								#loop through codons, find mismatch
								for codon in xrange(0, len(mut_gene), 3):
									if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
										mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
										wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
										#check for mismatch, report first one and stop
										#check if synonymous, non-synonymous
										if indel == False:
											if mut_aa != wt_aa:
												if mut_aa == '*':
													print >> f_out, '\t'.join(l + ['nonsense', genes[g][0], wt_aa + str(codon/3+1) + mut_aa]) 
												else:	
													print >> f_out, '\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
											else:
												print >> f_out, '\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
										
										#otherwise, print the first mismatch as an indel, then find the next stop
										else:
											#these would be in-frame
											if snp_len % 3 == 0:
												if snp_len > 0:
													print >> f_out, '\t'.join(l + ['in-frame-insertion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
											 	elif snp_len < 0:
													print >> f_out, '\t'.join(l + ['in-frame-deletion', genes[g][0], wt_aa + str(codon/3+1) + 'X'])
											else:
												#these are out of frame
												#translate sequence, find first stop
												trunc = translate_sequence(''.join(mut_gene))
												if '*' in trunc:
													if snp_len > 0:
														print >> f_out, '\t'.join(l + ['frameshift-insertion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa + "; Truncation:" + str(trunc.index('*')+1) + '/' + str(int(genes[g][5])) ])
													elif snp_len < 0:
														print >> f_out, '\t'.join(l + ['frameshift-deletion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa + "; Truncation:" + str(trunc.index('*')+1) + '/' + str(int(genes[g][5])) ])
												else:
													if snp_len > 0:
														print >> f_out, '\t'.join(l + ['frameshift-insertion', genes[g][0], str(codon/3+1) + mut_aa + "; Full-length"])
													elif snp_len < 0:
														print >> f_out, '\t'.join(l + ['frameshift-deletion', genes[g][0], wt_aa + str(codon/3+1) + mut_aa + "; Full-length"])
										break
		
								annotation = True
												
						#if snp isn't in a gene, check if it's in the upstream region of a gene
						elif snp_start in xrange(genes[g][1]-upstream,genes[g][1]):
							print >> f_out, '\t'.join(l + ["5'-upstream", genes[g][0], str(snp_start-genes[g][1])])
							annotation = True						
		
			#after checking all genes for genic or non-coding snps, check non-coding elements
			if annotation == False:
				#if snp isn't upstream of a gene, check if it's in a non-coding region of the genome
				for n in noncoding:
					#if correct chromosome
					if noncoding[n][1] == chrom_conv:
						if snp_start >= noncoding[n][3] and snp_start <= noncoding[n][4]:
							print >> f_out, '\t'.join(l + [noncoding[n][2], noncoding[n][0], 'NA'])
							annotation = True		
				#if annotation is still false after check non-coding elements, then the snp is just intergenic
				if annotation == False:
					print >> f_out, '\t'.join(l + ['intergenic', 'NA', 'NA'])
		
		except StopIteration:
			f_out.close()
			return
									
									
	"""OTHER USEFUL FUNCTIONS"""

def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':17, 'mitochondrion':17,
				'to':17, 'on':17, 'M':17}
	
	chrom_number = chrom_number.split('chr')[-1]
	
	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:			
		return chrom_conv[chrom_number]

def complement(base):
	compbase = []
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	for i in xrange(len(base)):
		try:
			compbase.append(comp[base[i].upper()])
		except KeyError:
			compbase.append('N')
	return ''.join(compbase)

def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
             'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
             'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
             'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
             'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
             'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
             'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
             'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
             'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
             'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
             'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
             'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
             'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
             'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
             'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
             'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' }
	try:
		return lookup[codon.lower()]
	except KeyError:
		return 'X'

# translate DNA -> amino acid
def translate_sequence(seq):
	translated_seq = ''
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq

# FASTA parser
def readFASTA(f):
	f_out = {}
	k = ''
	val = ''
	
	#find first line that isn't a comment
	comment = True
	while comment == True:
		l = f.next()
		if l.startswith("#") != True:
			#initialize the first key
			k = l.strip()[1:]
			comment = False
	
	while True:
		try:
			l = f.next()
						
			if l.startswith('>'):
				f_out[k] = val
				k = l.strip()[1:]
				val = ''			
			else:
				val += l.strip()

		except StopIteration:
			val += l.strip()
			f_out[k] = val
			return f_out


if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	from copy import copy
	from random import randint
	
	parser = OptionParser()
	parser.add_option('-m', '--mutations', action = 'store', type = 'string', dest = 'inputfile', help = 'BED file with mutations')
	parser.add_option('-c', '--coding', action = 'store', type = 'string', dest = 'sequences', help = 'fasta file of coding sequences')
	parser.add_option('-n', '--non-coding', action = 'store', type = 'string', dest = 'noncoding', help = 'gff file containing non-coding regions')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'genome', help = 'fasta file containing genome sequence')
	parser.add_option('--upstream', action = 'store', type = 'int', dest = 'upstream', help = "Distance (in bp) upstream of gene start codons to include in 5'-upstream annotation (default=500)", default=500)
	
	(option, args) = parser.parse_args()

	main(option.inputfile, option.sequences, option.noncoding, option.genome, option.upstream)

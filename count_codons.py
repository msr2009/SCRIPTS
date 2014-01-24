"""
count_codons.py

Matt Rich, 4/2013
"""


def main(gene_list, codonpair_list, orf_fasta): 
	
	count_codons = { x: [0,0,0,0] for x in codonpair_list }

	#read fasta file, save genes in gene list
	genes = {}
	for f in SeqIO.parse(open(orf_fasta,'r'), 'fasta'):
		geneid = f.description.split()[0]
		if geneid in gene_list:
			genes[geneid] = [f.seq.tostring(), True]
		else:
			genes[geneid] = [f.seq.tostring(), False]
			
	
	for g in genes:
		#get codons
		codons = [ genes[g][0][i:i+3] for i in range(0,len(genes[g][0])-3,3) ] 
		#loop through codons, looking for codon pairs of interest
		for p in codonpair_list:
			pair = [p[0:3], p[3:6]]
			for c in range(len(codons)-1):
				if codons[c] == pair[0] and codons[c+1] == pair[1]:
					#if the codon pair is in the gene, count it
					count_codons[p][1] += 1
					#and if the gene is one of the subset, count that
					if genes[g][1] == True:
						count_codons[p][0] += 1
				if codons[c] == pair[1] and codons[c+1] == pair[0]:
					#count the reverse codon pairs
					count_codons[p][3] += 1
					if genes[g][1] == True:
						count_codons[p][2] += 1
	
	print '\t'.join(['codon_pair','subset_genes_forward','total_genes_forward','subset_genes_reverse','total_genes_reverse'])
	for c in count_codons:
		print '\t'.join([c, str(count_codons[c][0]), str(count_codons[c][1]), str(count_codons[c][2]), str(count_codons[c][3])])
		
					
if __name__ == "__main__":
	from Bio import SeqIO
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('--codons', action = 'store', type = 'string', dest = 'codons', help = 'comma-delimited list of codon-pairs')
	parser.add_option('--orfs', action = 'store', type = 'string', dest = 'orfs', help = 'fasta file of orfs')
	parser.add_option('--genes', action = 'store', type = 'string', dest = 'genes', help = 'file containing subset of genes to discriminate against, one per line')
	
	(option, args) = parser.parse_args()
	
	gene_list = []
	for line in open(option.genes,'r'):
		gene_list.append(line.strip())
	
	codon_list = []	
	if option.codons == "all":
		all_codons = []
		bases = ['A','C','T','G']
		for j in bases:
			for k in bases:
				for l in bases:
					all_codons.append(j+k+l)
		for m in all_codons:
			for n in all_codons:
				codon_list.append(m+n)
	else:
		codon_list = option.codons.split(',')
		
	main(gene_list, codon_list, option.orfs)

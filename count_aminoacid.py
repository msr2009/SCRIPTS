#/usr/bin/python

"""
count_aminoacid.py

Takes input of fasta file containing coding sequences, outputs amino acid frequency.

--translate asks if the file has translated sequences.  If the FASTA file is of protein sequences, then -t should 
be 'T', otherwise -t = 'F'.

python count_aminoacid.py --input INPUTFASTA --translate SEQSTRANSLATED?
 
"""

def main(inputfile, translate):
	
	#initialize dictionary to store amino acid counts
	aminoacids = {}
	
	#count amino acids
	#translate sequences
	if translate == False:
		for record in SeqIO.parse(inputfile, 'fasta'):
			for aa in translate_sequence(str(record.seq)):
				if aa in aminoacids:
					aminoacids[aa] += 1
				else:
					aminoacids[aa] = 1
	#don't translate sequences
	else:
		for record in SeqIO.parse(inputfile, 'fasta'):
			for aa in str(record.seq):
				if aa in aminoacids:
					aminoacids[aa] += 1
				else:
					aminoacids[aa] = 1
	
	#print counts
	for i in sorted(aminoacids.keys()):
		print '\t'.join([i, str(aminoacids[i])])


#lookup table for translation script
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
	return lookup[codon.lower()]

# translate DNA -> amino acid
def translate_sequence(seq):
	translated_seq = ''
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq

if __name__ == '__main__':

	from optparse import OptionParser
	from Bio import SeqIO
	
	parser = OptionParser()
	parser.add_option('-i', '--input', action = 'store', type = 'string', dest = 'inputfasta', help = 'FASTA file of coding sequences')
	parser.add_option('-t', '--translate', action = 'store', type = 'string', dest = 'translate', help = 'Are the sequences translated? T/F')
	(option, args) = parser.parse_args()
	
	translate = True
	if option.translate in set(['F', 'f', 'False', 'false', 'no']):
		translate = False
	
	main(option.inputfasta, translate)
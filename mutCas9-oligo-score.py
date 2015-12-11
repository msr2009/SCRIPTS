"""


comments

Matt Rich, DATE
"""

from fasta_iter import fasta_iter
from parseBLOSUM import parseBLOSUM
import numpy as np
import warnings 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



def main(grna, seq, blosum, mdist, fasta_name, forceC, pre_out):
	#open an output file
	dat_out = open(pre_out+".cda1.dat", "w")
	
	targets = []
	score_min = [float('Inf')]
	score_max = [-float('Inf')]

	if grna == None:
		#first find all targets
		targets = PAMfinder(seq)
	else:
		#read gRNAs from file
		for line in open(grna, "r"):
			#skip any headers
			if line.startswith("#") or line.startswith("chr"):
				continue
			#then read the position data
			l = line.strip().split("\t")
			targets.append([l[1], l[2], None, l[3]])

	for t in targets:
		blos_scores = [] #store all our blosum scores here
		#identify the residue and position of each site in gRNA
		strand = 1
	
		for p in range(len(t[2])):
			if t[3] == "+":
				res = (t[0]+p)/3
				pos = (t[0]+p)%3
				aa = seq[res*3:res*3+3]
				if forceC and t[2][p] != "C":
					blos_scores.append([np.nan, np.nan, np.nan, np.nan])				
				else:
					blos_scores.append(blosAA(aa, pos, blosum))				
					
			if t[3] == "-":
				res = (t[0]-1-p)/3 #<<<have to decrement these by 1
				pos = (t[0]-1-p)%3
				aa = seq[res*3:res*3+3]
				strand = -1
				if forceC and t[2][p] != "G":
					blos_scores.append([np.nan, np.nan, np.nan, np.nan])				
				else:	
					blos_scores.append(blosAA(aa, pos, blosum))
				

		#transpose blos_scores to same orientation as mutation distributions
		blos_scores = np.array(blos_scores).transpose()
		
		#normalize blos_scores based on empirical mutation distribution
		if mdist:
			if t[3] == "+":
				blos_scores = np.multiply(blos_scores, mdist[0])
			else:
				blos_scores = np.multiply(blos_scores, mdist[1])
	
		if np.nanmin(np.nansum(blos_scores, axis=0)) < score_min:
			score_min = np.nanmin(np.nansum(blos_scores, axis=0))
		if np.nanmax(np.nansum(blos_scores, axis=0)) > score_max:
			score_max = np.nanmax(np.nansum(blos_scores, axis=0))
	
		printAndPlot(blos_scores, t, fasta_name, strand, dat_out)

	plt.axis([0, len(seq), score_min-.5, score_max+.5])
	plt.xlabel("Coding sequence (bp)")
	plt.ylabel("Normalized BLOSUM64 score")
	plt.title(fasta_name)

	plt.savefig(pre_out+".cda1.pdf")
	plt.clf()
	dat_out.close()
			
#print and plot scores
def printAndPlot(mat, t, name, strand, outfile):
		mat_sum = np.nansum(mat,axis=0)
		blos_sums =	",".join([ str(x) for x in mat_sum ]) 
		if not np.isnan(np.nansum(mat)):
			print >> outfile, "\t".join([str(x) for x in [name] + t + \
										[blos_sums, np.nansum(mat)]])			
			#also add points to plot
			if strand == 1:
				plt.scatter(range(t[0], t[1], strand), mat_sum, c="r")
			elif strand == -1:
				plt.scatter(range(t[0], t[1], strand), mat_sum, c="b")
	

# lookup table for codon translation
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

#reverse complement a sequence
def reverse_complement(seq):
	complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
	return "".join([complement[x.upper()] for x in seq[::-1]])

#find all possible guide target sequences
def PAMfinder(seq, pams=set(["AGG","CGG","GGG","TGG"])):
	targets = []
	seq2 = reverse_complement(seq)
	for i in range(20,len(seq)):
		if seq[i:i+3] in pams:
			targets.append([i-20, i, seq[i-20:i], "+"])
		if seq2[i:i+3] in pams:
			targets.append([len(seq)-i+20, len(seq)-i,
							seq2[i-20:i], "-"])
	return targets

#read mutation distribution
def readMutDist(f):
	if f:
		m = []
		n = []
		for line in open(f, "rU").readlines()[1:]: 
			m.append([ float(x) for x in line.strip().split()[1:] ])
		for i in [3,2,1,0]:
			n.append(m[i])
		return (np.array(m), np.array(n))
	else:
		return None

#mutagenize amino acid
def mutAA(aa, pos, base):
	aa = list(aa)
	aa[pos] = base 		
	return "".join(aa) 

#lookup blosum score for all mutations at a site
def blosAA(aa, pos, blos):
	bases = ["A", "C", "G", "T"]
	bs = [np.nan, np.nan, np.nan, np.nan]
	wt_aa = translate_sequence(aa)
	mut_aa = [ translate_sequence(mutAA(aa, pos, b)) for b in bases ]
	for t in range(4):	
		if aa[pos] != bases[t]:
			bs[t] = blos[wt_aa][mut_aa[t]]
	return bs					

#	return [ blos[wt_aa][translate_sequence(mutAA(aa, pos, b))] for b in bases ] 

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	#take input as either list of pre-computed gRNAs or greedily find all PAMs
	parser.add_argument('--grna', action = 'store', type = str, dest = 'grna',
		help = "optional BED file of pre-designed set of gRNAs", default = None)
	parser.add_argument('--seq', action = 'store', type = str, dest = 'sequence', 
		help = "FASTA file containing sequence to design targets for")
	parser.add_argument('--blosum', action = 'store', type = str, 
		dest = 'blosum', help = "blosum matrix to use (from NCBI)")
	parser.add_argument('-m', '--mutdist', action = 'store', type = str, 
		dest = 'md', help = "distribution of mutations gRNA", 
		default = None)
	parser.add_argument('--forceC', action = 'store_true', 
		dest = 'forceC', help = "only calculate for mutation from cytosine",
		default = False)
	parser.add_argument('--out', action = 'store', type = str, 
		dest = 'outfile', help = "prefix for output (defaults to first word \
							of FASTA entry")
	args = parser.parse_args()

	outfile = ""

	with warnings.catch_warnings():
		warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
		warnings.filterwarnings('ignore', '', FutureWarning)

		for s in fasta_iter(args.sequence):
			fasta_name = s[0].strip().split()[0]
			print fasta_name
			if args.outfile != None:
				main(args.grna, s[1], parseBLOSUM(args.blosum), 
					readMutDist(args.md), fasta_name, 
					args.forceC, args.outfile)	
			else:
				main(args.grna, s[1], parseBLOSUM(args.blosum), 
					readMutDist(args.md), fasta_name, 
					args.forceC, fasta_name)	
			


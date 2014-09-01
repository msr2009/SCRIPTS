"""
uORFs.py	

Matt Rich, 8/14
"""

def main(seq, offset):
	#print header for output file
	print "\t".join(["start", "seq", "length", "uORF_seq", "note"])
	#loop through sequence	
	for i in range(len(seq)-3):
		#find site where you could mutate to ATG
		if uORFs(seq[i:i+3]):
			out = [str(i+offset), seq[i:i+3]]
			#translate what you can and find shortest uORF
			trans = translateSequence(seq[i:])
			uORF = trans.split("*")[0]
			uORF = "M" + uORF[1:]
			if uORF.endswith(">"):	#orf is out of frame
				out += ["NA", "NA", "out-of-frame"]
			elif len(trans.split("*")) == 1:	#in frame fusion
				out += [str(len(uORF)), uORF, "in-frame"]
			else:
				out += [str(len(uORF)), uORF, "NA"]
			#print output
			print "\t".join(out)

def translateSequence(seq):
	translated_seq = ""
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	if len(seq) % 3 != 0:
		translated_seq += ">"
	
	return translated_seq
 
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

def uORFs(mer):
	s = set(["CTG", "GTG", "TTG", "AAG", "ACG", "AGG", "ATT", "ATA", "ATC"])
	if mer in s:
		return True
	else:
		return False

if __name__ == "__main__":
	
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('--seq', action = 'store', type = 'string', dest = 'seq', help = "sequence to search") 
	parser.add_option('--offset', action = 'store', type = 'int', dest = 'offset', help = "sequence offset", default=0) 
	(option, args) = parser.parse_args()
	
	main(option.seq.upper(), option.offset)	


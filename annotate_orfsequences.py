"""
annotate_orfsequences.py

compares two files -- a wildtype sequence in FASTA format, and a BED file containing mutations of the other sequence

Matt Rich 5/2013

"""

def main(wildtype, variants, start, stop, typecol):
		
	#open wildtype sequence, save it.
	wt_seq = readFASTA(open(wildtype,'r')).values()[0]
	
	#open output file
	f_out = open(variants + '.annotated', 'w')
	
	#for each variant sequence:
	for line in open(variants,'r'):
		
		l = line.strip().split('\t')
		#this should be, at least in the case of the parseCross_matchOutput.py output:
		#[WTsequenceID, start, stop, ref, obs, qual, type (S,D-n,I-n)]
		
		#print and skip the header
		if "start" in line:
			print >> f_out, "\t".join(l + ['annotation', 'aa_position', 'aa_type'])
			continue
		
		#if the position is outside of the start-stop range, go to the next mutation
		if int(l[1]) < start or int(l[1]) > stop:
			continue
		
		#which codon is it?
		codon = (int(l[1])-1)/3
		wt_aa = lookup_codon(wt_seq[3*codon:3*codon+3])
		
		#if it's a substitution
		if l[typecol] == "S" or l[typecol] in set(["A","C","T","G"]):
			mut_codon = list(copy(wt_seq[3*codon:3*codon+3]))
			mut_codon[ int(l[1])%3-1 ] = l[4]
			mut_aa = lookup_codon(''.join(mut_codon))
			
			muttype = ""
			if mut_aa == wt_aa:
				muttype = "S"
			elif mut_aa != wt_aa:
				if mut_aa == "*":
					muttype = "N"
				else:
					muttype = "NS"
			
			mutinfo = [wt_aa+str(codon+1)+mut_aa, str(codon+1), muttype]
			
			print >> f_out, '\t'.join(l + mutinfo)
		
		#if the mutation is an deletion
		elif l[typecol].startswith("D"):
			#figure out where it starts, where the next stop is
			mut_seq = list(copy(wt_seq))
			#remove base(s)
			del mut_seq[int(l[1])-1:int(l[1])+int(l[typecol].split('-')[1])-1]
			t_mut_seq = translate_sequence("".join(mut_seq))			
			try:
				nextstop = t_mut_seq.index("*")
				print >> f_out, '\t'.join(l + [wt_aa+str(codon+1)+t_mut_seq[codon+1] + ", next-stop:"+str(nextstop+1), str(codon+1), "F"])				
			except ValueError:
				print >> f_out, '\t'.join(l + [wt_aa+str(codon+1)+t_mut_seq[codon+1] + ", Full-length", str(codon+1), "F"])
			
		elif l[typecol].startswith("I"):
			mut_seq = list(copy(wt_seq))
			#add bases
			if l[4].startswith("I"):
				mut_seq[int(l[1])-1] = l[4].split('-')[1]
			else:
				mut_seq[int(l[1])-1] = l[4]
			t_mut_seq = translate_sequence("".join(mut_seq))			
			try:
				nextstop = t_mut_seq.index("*")
				print >> f_out, '\t'.join(l + [wt_aa+str(codon+1)+t_mut_seq[codon+1]+", next-stop:"+str(nextstop+1), str(codon+1), "F"])	
			except ValueError:	
				print >> f_out, '\t'.join(l + [wt_aa+str(codon+1)+t_mut_seq[codon+1] + ", Full-length", str(codon+1), "F"])
	
	f_out.close()
				
			
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
			f_out[k] = val
			break
			
	return f_out

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
	return ''.join([ lookup_codon(seq[b:b+3]) for b in xrange(0, len(seq)-3, 3) ])

			
if __name__ == "__main__":
	
	from optparse import OptionParser
	from copy import copy
	
	parser = OptionParser()
	parser.add_option('-w', '--wt', action = 'store', type = 'string', dest = 'wt', help = 'FASTA file of wildtype sequence')
	parser.add_option('-v', '--variants', action = 'store', type = 'string', dest = 'var', help = 'BED file of variant sequences to compare to wildtype')
	parser.add_option('--pos-range', action = 'store', type = 'string', dest = 'pos', help = 'Comma-delimited start and stop (start,stop) positions. All mutations outside of these positions are ignored.')
	parser.add_option('--counts-input', action = 'store_true', dest = 'countsinput', help = 'input is from mapBEDcounts', default=False)
	(option, args) = parser.parse_args()
	
	start = float('-inf')
	stop = float('inf')
	if option.pos != None:
		[start, stop] = [ int(x) for x in option.pos.split(',')]
	
	typecol = 6
	if option.countsinput == True:
		typecol = 4
		
	main(option.wt, option.var, start, stop, typecol)

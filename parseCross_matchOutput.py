"""
parseCross_matchOutput.py

Parses cross_match alignment mismatch output, can make comparisons to original sequence

Matt Rich, 03/2013
"""

def main(cm_output, ref_fasta, offset, printBED, cutoff, noheader, verbose):
	#read fasta reference sequence, make dictionary
	ref_seq = SimpleFastaParser(open(ref_fasta, 'r'))
	
	#read cross_match output
	cm = open(cm_output, 'r')
	
	parsing = False
	pos = 0
	discrepancies = []
	mutation = ""
	read = ""
	length = 0
	disc_type = ""
	mapped_ref = ""
	comp = False
	header = True
	if noheader == True:
		header = False
	
	pileup = [0]*(len(ref_seq[ref_seq.keys()[0]]))

	while True:	
		line = cm.readline()	
		if line.startswith('Maximal single base matches') == True:
			while line.strip() != '':
				line = cm.readline()
			parsing = True
			continue

		#need to print out the last entry
		if 'matching entries (first file)' in line:			
			if printBED == "counts":
				print discrepancies[0][0]
			
			elif printBED == "bed":	
				for d in discrepancies[1:]:
					#if d is deletion:
					#deleted base is one base BEFORE position given by cross_match
					if d[0].startswith('D') == True:
						print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+int(d[0].split('-')[1])), d[3], d[2], d[4], d[0]])
						#print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+int(d[0].split('-')[1])), ref_seq[discrepancies[0][2]][ d[1]-offset-1:d[1]-offset + int(d[0].split('-')[1])+1 ], d[2], d[4], d[0]])
						
					#if d is insertion
					elif d[0].startswith('I') == True:
						print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+int(d[0].split('-')[1])), d[3], d[3] + d[2], d[4], d[0]])
					
					#if d is substitution
					elif d[0] == 'S':
						print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+len(d[2])), d[3], d[2], d[4], d[0]])
			break	
						
		if parsing == True:

			if line.strip() != '':
				l = line.strip().split()					
				#if it's the first line:
				try:
					int(l[0])
					#if the first line has an int at it's first position, then it's a new read					
					if discrepancies != []:
						#print out all the stored values, then reinitialize
						#readID, mutations (e.g., 45,128,533-A,C,X), types, length, total mutations, mapped_ref
						if printBED == "counts": 
							print discrepancies[0]

						elif printBED == "bed":
							if header == True:
								print '\t'.join(['chr', 'start', 'end', 'ref', 'obs', 'qual', 'type'])
								header = False
							
							for d in discrepancies[1:]:
								#if d is deletion:
								#deleted base is one base BEFORE position given by cross_match
								if d[0].startswith('D') == True:
									print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+int(d[0].split('-')[1])), d[3], d[2], d[4], d[0]])
									#print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+int(d[0].split('-')[1])), ref_seq[discrepancies[0][2]][ d[1]-offset-1:d[1]-offset + int(d[0].split('-')[1])+1 ], d[2], d[4], d[0]])
									
								#if d is insertion
								elif d[0].startswith('I') == True:
									print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+int(d[0].split('-')[1])), d[3], d[3] + d[2], d[4], d[0]])
								
								#if d is substitution
								elif d[0] == 'S':
									print '\t'.join([discrepancies[0][2], str(d[1]-offset), str(d[1]-offset+len(d[2])), d[3], d[2], d[4], d[0]])
						
						else:
							mut_string, type_string = mutLocs(discrepancies[1:], offset)
							print '\t'.join([discrepancies[0][0], mut_string, type_string, discrepancies[0][1], str(len(discrepancies[1:])), discrepancies[0][2]])							
												
						discrepancies = []					
						
					#read new data					
					read = l[4]
					length = l[6]
					mapped_ref = l[-4]
					if 'C' in l:
						comp = True
					else:
						comp = False
					
					#append information about mapped read
					discrepancies.append( [read, length, mapped_ref] )
				
					#add to pileup
					#have to get alignment endpoints first
					p_indices = l[-3:]
					#remove value inside parentheses
					for i in range(3):
						if p_indices[i].startswith("(") == True:
							del p_indices[i]
							break
					
					pileup = incrementPileup(pileup, int(p_indices[0]), int(p_indices[1])) 
				
				#otherwise, make a list of the relevant information from the discrepancy lines
				except ValueError:
					pos = int(l[3])
					disc_type = l[0]
					quality = l[2].split('(')[1][0:-1]
					if int(quality) < cutoff:
						continue
						
					if comp == True:
						mutation = complement(l[2].split('(')[0])
					else:
						mutation = l[2].split('(')[0]
					if disc_type != "E":
						if disc_type == 'D' or disc_type == 'I':
							disc_type += "-1"
						try:
							discrepancies.append( [disc_type, pos, mutation, ref_seq[discrepancies[0][2]][pos-1], quality] )						
						except KeyError:
							pass
	
	#print out pileup into {FILENAME}.pileup
	pileup_out = open(cm_output.rstrip(".crossmatch")+"a.pileup", "w")
	for p in range(len(pileup)):
		print >> pileup_out, "\t".join([str(p-offset), str(pileup[p])])
	pileup_out.close()

# increment pileup list
def incrementPileup(pileup, start, stop):
	for i in range(min([start,stop])-1, max([start,stop])):
		pileup[i] += 1
	return pileup

# create mutation location string
def mutLocs(mut_lists, offset):
	muts = []
	locs = []
	types = []
	for i in range(len(mut_lists)):
		
		locs.append(str(int(mut_lists[i][1])-offset))
		t = mut_lists[i][0]
		types.append(t)
		if t[0] == 'D':
			if len(t) == 1:
				muts.append('X')
			else:
				muts.append(''.join(['X' for i in range(len(t)-2)]))
		elif t[0] == 'I':
			muts.append('i' + mut_lists[i][2])
		else:
			muts.append(mut_lists[i][2])
	return ','.join(locs) + '-' + ','.join(muts), ','.join(types)

#this is biopython's code
def SimpleFastaParser(handle):
    """Generator function to iterator over Fasta records (as string tuples).

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> for values in SimpleFastaParser(open("Fasta/dups.fasta")):
    ...     print(values)
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    """
    f_out = {}
    
    #Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = handle.readline()

        #Remove trailing whitespace, and any internal spaces
        #(and any embedded \r which are possible in mangled files
        #when not opened in universal read lines mode)
        f_out[title] = "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return f_out # StopIteration

def complement(seq):
	complementary_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
	return ''.join([complementary_dict[base] for base in seq])


if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--fasta', action = 'store', type = 'string', dest = 'fasta', help = 'FASTA file of reference sequence')
	parser.add_option('-c', '--cross_match_output', action = 'store', type = 'string', dest = 'cross_match', help = 'cross_match output file (must include discrepancies information)')
	parser.add_option('--offset', action = 'store', type = 'int', dest = 'offset', help = 'offset used in calculating position (default = 0)', default = 0)
	parser.add_option('--output', action = 'store_true', dest = 'output_BED', help = 'how to output parsed data? options:bed, counts', default = "bed")
	parser.add_option('--cutoff', action = 'store', dest = 'qual_cutoff', type = 'int', help = 'lower quality cutoff to count mismatch as real', default = 0)
	parser.add_option('--no_header', action = 'store_true', dest = 'noheader', help = 'do not print header in BED file', default=False)
	parser.add_option('-v', '--verbose', action = 'store_true', dest = 'verbose', help = 'print status every so often', default=False)

	(option, args) = parser.parse_args()

	main(option.cross_match, option.fasta, option.offset, option.output_BED, option.qual_cutoff, option.noheader, option.verbose)

import math
import datetime
import re

''' define a series of functions for standard operations on genetic sequences: '''

# define a reverse complementation dictionary:
revcomp_dict = {
	'A':'T',
	'C':'G',
	'G':'C',
	'T':'A'}

# define a dictionary of the IUPAC code:
iupac_dict = {
	'A':'[A]',
	'C':'[C]',
	'G':'[G]',
	'T':'[T]',
	'M':'[AC]',
	'R':'[AG]',
	'W':'[AT]',
	'S':'[CG]',
	'Y':'[CT]',
	'K':'[GT]',
	'V':'[ACG]',
	'H':'[ACT]',
	'D':'[AGT]',
	'B':'[CGT]',
	'X':'[ACGT]',
	'N':'[ACGT]'}


# define a function that yields the reverse complement of a DNA sequence:
def revcomp(sequence):
	if sequence == '':
		return sequence
	else:
		sequence = sequence.upper()
		sequence = list(sequence)
		sequence.reverse()
		output = ""
		for char in sequence:
			if char == "A":
				output = output + "T"
			if char == "C":
				output = output + "G"
			if char == "G":
				output = output + "C"
			if char == "T":
				output = output + "A"
			if char == "N":
				output = output + "N"
		return output
	
	
# define a function that returns the indexes of a string and its rev.comp in a sequence:
def find(string, sequence):
	indexes = []
	string = string.upper()
	sequence = sequence.upper()
	# add sense instances:
	k = 0
	segments = sequence.split(string)
	segments.pop(len(segments)-1)
	for segment in segments:
		k += len(segment)
		if len(indexes) == 0:
			indexes.append(k)
		else:
			if indexes.count(k) == 0:
				indexes.append(k)
		k += len(string)
		
	# add antisense instances (rev.comp):
	k = 0	
	segments = sequence.split(gs.revcomp(string))
	segments.pop(len(segments)-1)
	for segment in segments:
		k += len(segment)
		if len(indexes) == 0:
			indexes.append(k)
		else:
			if indexes.count(k) == 0:
				indexes.append(k)
		k += len(string)
	
	return indexes


# define a function that returns the indexes of a string and its rev.comp in a sequence dictionary:
def dictfind(string, sequences):
	dictionary = {}
	string = string.upper()
	for key in sequences:
		sequence = sequences[key].upper()
		indexes = gs.find(string, sequence)
		if not indexes == []:
			dictionary[key] = indexes
	return dictionary
	
	
# define a function to count the incidences of a string and its rev.comp within a sequence:
def count(string, sequence):
	counts = 0
	string = string.upper()
	sequence = sequence.upper()
	counts += sequence.count(string)
	counts += sequence.count(gs.revcomp(self,string))
	return counts


# define a function to count the incidences of a string and its rev.comp within a sequence dictionary:
def dictcount(string, dictionary):
	counts = 0
	string = string.upper()
	for key in dictionary:
		counts += dictionary[key].upper().count(string)
		counts += dictionary[key].upper().count(gs.revcomp(self,string))
	return counts


# define a function to build a dictionary of restriction enzymes: 
def build_rebase(rebase_bionet):
	re_dict = {}
	start = False
	for line in rebase_bionet:
		if line.find('Rich Roberts') >= 0:
			start = True
			line = rebase_bionet.next()
		if start == True and len(line) > 10:
			buffer = line.split()
			re_dict[buffer[0]] = buffer[-1].replace('^', '')
	return re_dict


# define a function to find IUPAC defined sites in a sequence: 
def find_sites(input, set, sequence):
	iupac_dict = {
	'A':'[A]',
	'C':'[C]',
	'G':'[G]',
	'T':'[T]',
	'M':'[AC]',
	'R':'[AG]',
	'W':'[AT]',
	'S':'[CG]',
	'Y':'[CT]',
	'K':'[GT]',
	'V':'[ACG]',
	'H':'[ACT]',
	'D':'[AGT]',
	'B':'[CGT]',
	'X':'[ACGT]',
	'N':'[ACGT]'}
	
	site = set[input]
	pattern = ''
	positions = []
	for i in site:
		pattern += iupac_dict[i]
	searchpattern = re.compile(pattern)
  	sites = searchpattern.findall(sequence)
	temppos = searchpattern.finditer(sequence)
	for i in temppos:
		begin, end = i.span()
		positions.append(begin)
	return sites, positions

# specify a dictionary of amino acid single-letter and three-letter codes:
aa_1code_3code_dict = {
	"A":"Ala",
	"C":"Cys",
	"D":"Asp",
	"E":"Glu",
	"F":"Phe",
	"G":"Gly",
	"H":"His",
	"I":"Ile",
	"K":"Lys",
	"L":"Leu",
	"M":"Met",
	"N":"Asn",
	"P":"Pro",
	"Q":"Gln",
	"R":"Arg",
	"S":"Ser",
	"T":"Thr",
	"V":"Val",
	"W":"Trp",
	"Y":"Tyr",
}

# specify a dictionary of amino acid single-letter and three-letter codes:
aa_1code_name_dict = {
	"A":"Alanine",
	"C":"Cysteine",
	"D":"Aspartic Acid",
	"E":"Glutamic Acid",
	"F":"Phenylalanine",
	"G":"Glycine",
	"H":"Histidine",
	"I":"Isoleucine",
	"K":"Lysine",
	"L":"Leucine",
	"M":"Methionine",
	"N":"Asparagine",
	"P":"Proline",
	"Q":"Glutamine",
	"R":"Arginine",
	"S":"Serine",
	"T":"Threonine",
	"V":"Valine",
	"W":"Tryptophan",
	"Y":"Tyrosine",
}

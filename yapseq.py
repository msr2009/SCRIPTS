import math, gs


''' define a series of functions of general applicability to processing yapseq data: '''

aa_charged_positive = ['R','H','K']
aa_charged_negative = ['D','E']
aa_polar_uncharged = ['S','T','N','Q']
aa_hydrophobic = ['A','I','L','M','F','W','Y','V']
aa_special = ['C','G','P']

aas = aa_charged_positive + aa_charged_negative + aa_polar_uncharged + aa_hydrophobic + aa_special

aa_small = ['A','C','D','G','N','P','S','T','V']
aa_large = list(set(aas).difference(set(aa_small)))

''' define a function to build a dictionary of YapSeq input lines '''
def build(input_file):
	output_dict = {}
	input_key = open(input_file)
	header = input_key.readline()
	input_line = input_key.readline()
	while input_line:
		input_items = input_line.rstrip('\n').split('\t')
		lineID = input_items[0]
		if not lineID in output_dict:
			output_dict[lineID] = input_items[1:]
		else:
			print "ERROR: ID is not unique to line in file"
		input_line = input_key.readline()
	return output_dict

''' define a function to build a dictionary of the mutation frequency spectrum per position '''
def build_unlinked_dict(input_file, mode=""):
	output_dict = {}
	input_lines = open(input_file).readlines()

	items = input_lines.pop(0).rstrip().split('\t')[1:]	
	for input_line in input_lines:
		if mode == "NA-fix":
			input_line = input_line.replace("NA", "0")
		input_values = input_line.rstrip('\n').split('\t')
		pos = int(input_values.pop(0))
		output_dict[pos] = {}
		k = 0
		for item in items:
			output_dict[pos][item] = float(input_values[k])
			k += 1
	return output_dict
	
''' define a function to build a dictionary of the mutation frequency spectrum per position '''
def build_frequency_dict(input_file, mode=""):
	output_dict = {}
	input_lines = open(input_file).readlines()

	for pos in range(0, len(input_lines)-1):
		output_dict[pos] = {}
	
	items = input_lines.pop(0).rstrip('\n').split('\t')
	items.pop(0)
	
	for input_line in input_lines:
		if mode == "NA-fix":
			input_line = input_line.replace("NA", "0")
		input_values = input_line.rstrip('\n').split('\t')
		pos = int(input_values.pop(0))
		k = 0
		for item in items:
			output_dict[pos][item] = float(input_values[k])
			k += 1
	
	return output_dict


''' define a function to build a tally of sequences observed from filtered reads'''
def build_tally_dict(input_file):
	properties_dict = {}
	tally_dict = {}
	input_key = open(input_file)
	header = input_key.readline()
	input_line = input_key.readline()
	unique_data = {}
	unique_seqs = 0
	reads = 0
	
	while input_line:
		readID, sequence, match_count, mutation_count, mutation_location, mutation_identity, max_mutation_run = input_line.rstrip('\n').split('\t')
		reads += 1
		if not mutation_location in tally_dict:
			tally_dict[mutation_location] = {}
			properties_dict[mutation_location] = {}
		if mutation_identity in tally_dict[mutation_location]:
			tally_dict[mutation_location][mutation_identity] += 1
		else:
			tally_dict[mutation_location][mutation_identity] = 1
			properties_dict[mutation_location][mutation_identity] = [sequence, match_count, mutation_count, mutation_location, mutation_identity, max_mutation_run]
			unique_seqs += 1
		unique_data[reads] = unique_seqs
		input_line = input_key.readline()	
	return tally_dict, properties_dict, unique_data
	

''' define a function that builds a dictionary representation of a column of data '''
def build_value_dict(input_file, id_index, value_index, mode):
	value_dict = {}
	input_key = open(input_file)
	header = input_key.readline()
	input_line = input_key.readline()
	
	if mode == "str":
		def format(x):
			return x
	elif mode == "float":
		def format(x):
			return float(x)		
	elif mode == "int":
		def format(x):
			return int(x)	
			
	while input_line:
		input_items = input_line.rstrip('\n').split('\t')
		value_dict[input_items[id_index]] = format(input_items[value_index])
		input_line = input_key.readline()
	return value_dict


''' define a function that builds a dictionary representation of 'mapCounts.py' output '''
def build_norm_dict(input_file):
	return build_value_dict(input_file, 0, 7, "float")
	

''' define a function that builds a dictionary representation of 'mapCounts.py' output '''
def build_count_dict(input_file):
	return build_value_dict(input_file, 0, 7, "int")
	

''' define a function that builds a dictionary representation of log2ratios from 'mapRatios.py' output '''
def build_log2ratio_dict(input_file):
	return build_value_dict(input_file, 0, 7, "float")
		

''' define a function that normalizes the counts in a counts dictionary '''
def norm_count_dict(count_dict):
	norm_dict = {}
	total = 0
	for seqID in count_dict:
		total += count_dict[seqID]
	for seqID in count_dict:
		norm_dict[seqID] = float(count_dict[seqID])/total
	return norm_dict

''' define a function that returns the ratio between two frequency dictionaries '''
def gen_ratio_dict(A_dict, B_dict):
	ratio_dict = {}
	for seqID in set(A_dict).intersection(set(B_dict)):
		ratio_dict[seqID] = float(A_dict[seqID])/B_dict[seqID]
	return ratio_dict

''' define a function that estimates the ratio between two count dictionaries, uses fake 0.5 reads '''
def fix_ratio_dict(A_dict, B_dict, A_sum, B_sum):
	ratio_dict = {}
	for seqID in A_dict:
		A_norm = float(A_dict[seqID])/A_sum
		if seqID in B_dict:
			B_norm = float(B_dict[seqID])/B_sum
		else:
			B_norm = float(0.5)/B_sum
		ratio_dict[seqID] = float(A_norm)/B_norm
	return ratio_dict
	
''' define a function that returns the log2-ratio between two frequency dictionaries '''
def gen_log2ratio_dict(A_dict, B_dict):
	ratio_dict = {}
	for seqID in set(A_dict).intersection(set(B_dict)):
		ratio_dict[seqID] = math.log(float(A_dict[seqID])/B_dict[seqID],2)
	return ratio_dict

''' define a function that print amino acid codes for use with R: '''
def printAAtable(filename, stop=""):
	f = open(filename, "w")
	print >>f, "\t".join(["one.letter","three.letter","name"])
	names = ""
	codes = ""
	for aa in sorted(gs.aa_1code_3code_dict.keys()):
		print >>f, "\t".join([aa, gs.aa_1code_3code_dict[aa], gs.aa_1code_name_dict[aa].replace(" ","_")])
		names += '"' + gs.aa_1code_name_dict[aa] + '",'
		codes += '"' + gs.aa_1code_3code_dict[aa] + '",'
	if not stop == "":
		print >>f, "\t".join([stop,"Stop","Stop"])
		names += '"Stop",'
		codes += '"Stop",'	
	print names
	print codes
	f.close()
	
def printReference(infile):
	f1 = open("wildtype_1code.txt", "w")
	f3 = open("wildtype_3code.txt", "w")
	fn = open("wildtype_names.txt", "w")
	reference = open("../Python/input/" + infile).read().rstrip()
	print >>f1, "\t".join(["position", "variable", "value"])
	print >>f3, "\t".join(["position", "variable", "value"])
	print >>fn, "\t".join(["position", "variable", "value"])
	k = 0
	for aa in reference:
		print >>f1, "\t".join([str(k), aa, "0"])	
		print >>f3, "\t".join([str(k), gs.aa_1code_3code_dict[aa], "0"])
		print >>fn, "\t".join([str(k), gs.aa_1code_name_dict[aa].replace(" ","_"), "0"])
		k += 1
		
	f1.close()
	f3.close()
	fn.close()
import math
import datetime
import re

''' define a series of functions of general applicability: '''

# define a function to build a dictionary from a 'bed' file: 
def build(input_lines, id_column, value_column, mode='num'):
	dictionary = {}
	for line in input_lines:
		data = line.rstrip('\n').split('\t')
		feature, value = data[id_column], data[value_column]
		if mode == 'num':
			value = float(value)
		if dictionary.has_key(feature):
			print 'Error: feature already in dictionary!'
		dictionary[feature] = value
	return dictionary
	
''' define a function that builds a dictionary representation of a column of data '''
def build_file(input_file, id_column, value_column, mode):
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
		value_dict[input_items[id_column]] = format(input_items[value_column])
		input_line = input_key.readline()
	return value_dict

# combine items in two different lists of equal length into a list of item pairs:
def ziplists(list1, list2):
	pairs = []
	k = 0
	if not len(list1) == len(list2):
		print "Error: Lists are of unequal length!" 
	while k < len(list1):
		x1, x2 = list1[k], list2[k]
		pairs.append([x1,x2])
		k += 1
	return pairs
	
# print items in list to a specified output file:
def printlines(output_lines, output_file):
	f = open(output_file,'w')
	for output_line in output_lines:
		print >>f, output_line 
	f.close()

# print text to a specified output file:
def printtext(output_text, output_file):
	f = open(output_file,'w')
	print >>f, output_text 
	f.close()

# define a function to find the most common item in a list:	
def getmostcommon(ilist, iset):
	kmax = 0
	ks = []
	if len(iset) > 1:
		for item in iset:
			if ilist.count(item) > kmax:
				kmax = ilist.count(item)
				imax = item
		return imax
	elif len(iset) == 1:
		return list(iset)[0]
	else:
		print 'ERROR!'

# define a function to recover all items in a list of tab-delimited lines of items:
def tabrecover(input_list, index):
	output_list = []
	for line in input_list:
		items = line.split('\t')
		output_list.append(items[index])
	return output_list

# define a function to translate a dictionary from dict[id#1] --> dict[id#2]:
def translate(input_dict, translation_dict):
	output_dict = {}
	for key in input_dict:
		if key in translation_dict:
			output_dict[translation_dict[key]] = input_dict[key]
	return output_dict

# define a function to remove empty strings from list:
def clean(input_list):
	output_list = []
	for item in input_list:
		if not item == '':
			output_list.append(item)
	return output_list
	
""" Returns the keys of dictionary sorted by their values """
def valuesort(input_dict):
	items=input_dict.items()
	backitems=[[v[1],v[0]] for v in items]
	backitems.sort()
	return [ backitems[i][1] for i in range(0,len(backitems)) ]

""" Returns the keys of dictionary sorted by their values at list index """
def valuesort_index(input_dict, index, mode):
	items=input_dict.items()
	if mode == 'make.numeric':
		backitems=[[float(v[1][index]),v[0]] for v in items]
	else:
		backitems=[[v[1][index],v[0]] for v in items]
	backitems.sort()
	return [ backitems[i][1] for i in range(0,len(backitems)) ]

# define a function to reduce a dictionary to a set of keys:
def dictreduce(input_dict, input_list):
	output_dict = {}
	for key in input_list:
		output_dict[key] = input_dict[key]
	return output_dict
	
# define a function to reduce a list to a subset:
def listreduce(input_list, input_items, index):
	output_list = []
	for line in input_list:
		items = line.split('\t')
		if items[index] in input_items:
			output_list.append(line)
	return output_list

# define a function that determines whether two features overlap:
def testoverlap(qstart, qend, tstart, tend):
	if (qstart >= tstart and qstart <= tend) or (qend <= tend and qend >= tstart):
		return True
	else:
		return False
		
# define a function that returns the coordinate parts from a 'chr:start-end' coordinate string:
def getcoord(coord):
	chrm, poss = coord.split(':')
	start, end = map(int, poss.split('-'))
	return [chrm, start, end]

# define a function that counts the occurrences of elements in a list in a sequence:
def basecount(sequence, bases):
	x = 0
	for base in bases:
		x += sequence.count(base)
	return x
	
# define a function to remove numeric characters from a string: '''
def stripnumbers(string):
	for x in range(0,10):
		string = string.replace(str(x),"")
	return string

# define a function to invert a dictionary:
def dictinvert(input_dict):
	output_dict = {}
	for key in input_dict:
		print input_dict[key]
		output_dict[input_dict[key]] = key
	return output_dict
"""
compile_summaries.py

Compiles relevant information from the *.summary files of a folder, makes table

"""

def main(f_in):
	
	#get files
	summaries = []
	for i in os.listdir(f_in):
		if ".summary" in i and i != '.summary':
			summaries.append(i)
	
	mutations = ['total', 'coding-synonymous', 'coding-nonsynonymous', 'nonsense', "5'-upstream",  'intergenic']
	#get summary data
	all_data = {}
	for i in sorted(summaries):
		all_data[i.split('.')[0]] = getMutations(i, mutations)
	
	#print out table
	print '\t' + '\t'.join(sorted(all_data.keys()))
	
	for i in range(len(mutations)):
		print mutations[i] + '\t',
		for j in sorted(all_data.keys()):
			print str(all_data[j][mutations[i]]) + '\t',
		print '\n'
	
def getMutations(f_summary, m):	
	muts = {}
	#mutation types we want
	for i in m:
		muts[i] = 0
	
	for line in open(f_summary, 'r'):	
		l = line.strip().split('\t')
		if l[0] in muts:
			muts[l[0]] = l[1]
	return muts

if __name__ == '__main__':
			
	from optparse import OptionParser
	import os
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'f_input', help = 'folder containing summary files')
	(option, args) = parser.parse_args()
	
	main(option.f_input)

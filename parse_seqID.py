#/usr/bin/python

"""
parse_seqID.py

searches through counts or ratios files, prints lines to stdout that contain a given mutation or mutations

python parse_seqID.py -f FILE -m MUTATIONS

Matt Rich, 11/2011
"""

def main():
	
	#parse mutations to find
	mutations = set(option.mutations.strip().split(','))
	
	#open file
	line_count = 0
	for line in open(option.infile, 'r'):
		if line_count == 0:
			line_count += 1
			print line.strip()
			continue
		#parse seqID
		seq_muts = []
		sID = map(lambda x: x.split(','), line.strip().split('\t')[0].split('-'))
		#create list of mutations in sequence
		for i in range(len(sID[0])):
			seq_muts.append(sID[0][i]+sID[1][i])
		
		#compare seq_muts to mutations
		#print sequences that contain only mutations
		if len(seq_muts) == len(mutations):
			counter = 0
			for m in seq_muts:	
				if m in mutations:
					counter += 1
			if counter == len(mutations):
				print line.strip()
		#print all sequences that contain mutations		
		elif len(seq_muts) >= len(mutations) and option.contains == True:
			counter = 0
			for m in seq_muts:	
				if m in mutations:
					counter += 1
			if counter == len(mutations):
				print line.strip()
		line_count += 1

if __name__ == '__main__':
	
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'path to file')
	parser.add_option('-m', '--mode', action = 'store', type = 'string', dest = 'mutations', help = "comma-delimited list of mutations to find, e.g., '2P,16W,34A'")
	parser.add_option('-c', '--contains', action = 'store_true', dest = 'contains', default = False, help = "output all sequences containing mutation(s)? Default to output lines with only inputted mutations.")
	(option, args) = parser.parse_args()

	main()
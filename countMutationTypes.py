"""
countMutationTypes.py

counts transitions, transversion, and indels in parsed cross-match output BED file
"""

def main(file_list, score_indels, frequencies):

	data = {}
	#split list of files, load each one separately
	for f in file_list.split(','):
		data[f] = {"I":0, "D":0, "A>T":0, "A>C":0, "A>G":0, "T>A":0, "T>C":0, "T>G":0, "G>A":0, "G>C":0, "G>T":0, "T>A":0, "T>G":0, "T>C":0, "C>A":0, "C>T":0, "C>G":0}
		#read file, parse mutation types
		for line in open(f, 'r'):
			#skip header
			if line.startswith('chr') == True:
				continue
			
			l = line.strip().split('\t')
			if l[6][0] == "I" or l[6][0] == "D":
				if score_indels == True:
					data[f][l[6][0]] += 1
			else:	
				data[f][l[3].upper() + '>' + l[4].upper()] += 1
	
	print '\t'.join(["bedfile"] + sorted(data[f].keys()))			
	for d in sorted(data.keys()):
		line = [d]
		total = float(sum(data[d].values()))
		for x in sorted(data[d].keys()):
			if frequencies == False:
				line.append( str(data[d][x]) )		
			else:
				line.append( str(data[d][x]/total) )
		print '\t'.join(line)
				
				
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-b', '--bed', action = 'store', type = 'string', dest = 'bedfile', help = 'BED file containing all mutations (can be comma-delimited list of multiple files)')
	parser.add_option('--score_indels', action = 'store_true', dest = 'score_indels', help = 'Score indels, as well as substitutions', default = False)
	parser.add_option('--frequencies', action = 'store_true', dest = 'frequencies', help = 'Output frequencies instead of counts', default = False)
	(option, args) = parser.parse_args()

	main(option.bedfile, option.score_indels, option.frequencies)

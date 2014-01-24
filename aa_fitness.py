#/usr/bin/python

"""
aa_fitness.py

grabs fitness values for sequences with [0-N] amino acids within a region

python aa_fitness.py -f FITNESSFILE -o OUTPUTFOLDER -s STARTPOS -e ENDPOS -a AMINOACID(S)TOCOUNT -m MAXFLANKINGMUTATIONS

"""

def main():		
	
	#create list of lists to store data 
	diff = option.pos_end - option.pos_start + 1
	aas = {}
	for i in range(diff):
		aas[i] = []
	
	aa = []
	for a in sorted(option.aminoacid.strip().split(',')):
		aa.append(a)

	print option.infile

	#print out the list at the end	
	f_out = open(option.f_output_folder + option.infile.split('/')[-1] + '_fitness_' + str(option.pos_start) + '-' + str(option.pos_end) + '_' + ''.join(aa) + '_max' + str(option.max_mut), 'w')
	
	#print header
	print >> f_out, option.aminoacid.strip() + '\tfitness'
	
	#loop through file, put fitness scores in proper bins
	for line in open(option.infile, 'r'):
		if line.startswith('seqID') != True:
			l = line.strip().split('\t')
			
			#exclude variants with mutations outside region	
			if l[4] == 'NA':
				c = 0
				for a in aa:
					c += list(l[1])[option.pos_start:option.pos_end+1].count(a)
				print >> f_out, str(c) + '\t' + l[7]	
			else:
				muts = map(lambda x: int(x), l[4].split(','))
				outside_muts = set(muts).difference(set(range(option.pos_start, option.pos_end+1)))
				if len(outside_muts) > option.max_mut:
					continue	
				
			#count aa's and add to bucket
			c = 0
			for a in aa:
				c += list(l[1])[option.pos_start:option.pos_end+1].count(a)			

			print >> f_out, str(c) + '\t' + l[7] 

		
if __name__ == '__main__':
	from optparse import OptionParser
	import sys

	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'path to fitness file (Enrich ratio file)')
	parser.add_option('-o', '--output', action = 'store', type = 'string', dest = 'f_output_folder', help = 'path to output file')
	parser.add_option('--start', '-s', action = 'store', type = 'int', dest = 'pos_start', help = 'start point for counting')
	parser.add_option('--end', '-e', action = 'store', type = 'int', dest = 'pos_end', help = 'end point for counting')
	parser.add_option('--aa', '-a', action = 'store', type = 'string', dest = 'aminoacid', help = 'amino acid to count (comma-delimited for multiple amino acids)')
	parser.add_option('--max_mutation', '-m', action = 'store', type = 'int', dest = 'max_mut', help = 'maximum number of mutations outside defined sequence')
	(option, args) = parser.parse_args()

	main()

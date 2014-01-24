"""
batch_check_viability.py

after annotating genomes, checks viability and reports stats for all genomes in a folder.

python batch_check_viability.py -f FOLDER

"""

def main(folder, essentials_file, synthetic_file):
	#create dictionary to store data
	#data[file] = [viable, nssnp, ssnp, ss, intron, 5', inter-ele, intergenic]
	data = {}
	
	#get list of genome files
	f = os.listdir(folder)
	f = map(lambda x: folder + x, f)
	#loop through list of files,
	for i in f:
		if i.endswith('.annotated') == True:
			#check viability of genome
			filename = i.split('/')[-1]
			data[filename] = [genome_viability.main(i, essentials_file, synthetic_file)]
			#count stats of functional elements hit
			elements = []
			element_set = set('coding-synonymous', 'coding-nonsynonymous', 'intron', "5'-upstream", 'splice-site', 'intergenic')
			for line in open(i, 'r'):
				if line.startswith('chr') == True:
					continue
				l = line.strip().split('\t')
				if l[5] in element_set:
					elements.append(l[5])
				else:
					elements.append('intergenic-element')
			#append to dictionary entry
			data[filename].append(elements.count('coding-nonsynonymous'))
			data[filename].append(elements.count('coding-synonymous'))			
			data[filename].append(elements.count('splice-site'))
			data[filename].append(elements.count('intron'))
			data[filename].append(elements.count("5'-upstream"))
			data[filename].append(elements.count('intergenic-element'))
			data[filename].append(elements.count('intergenic'))
	
	#after getting all data parsed, output results
	f_out = open(folder + folder.split('/')[-2] + '_viability_stats.txt', 'w')
	#print header
	print >> f_out, '\t'.join(['file', 'total_snvs', 'coding-nonsynonymous', 'coding-synonymous', 'intron', "5'-upstream", 'intergenic-element', 'intergenic', 'viable?'])
	
	#first, calculate some overall results (going to make a double loop, but I want the totals at the top. oh well.)
	total = [0,0,0,0,0,0,0,0]
	for i in data:
		if data[i][0] == True:
			total[0] += 1
		for j in range(1,len(data[i])):
			total[j] += data[i][j]
	nfiles = float(len(f))
	print >> f_out, '\t'.join( ['total', str(int(nfiles)), str(total[1]/nfiles), str(total[2]/nfiles), str(total[3]/nfiles), str(total[4]/nfiles), 
								str(total[5]/nfiles), str(total[6]/nfiles), str(total[7]/nfiles), str(total[0]/nfiles) ] )
	
	#then, loop through again, printing out the data for each genome
	for i in data:
		print >> f_out, '\t'.join( [i.split('/')[-1], str(sum(data[i][1:])), str(data[i][1]), str(data[i][2]), str(data[i][3]),
								str(data[i][4]), str(data[i][5]), str(data[i][6]), str(data[i][7]), str(data[i][0]) ] )
	
	f_out.close()

			
if __name__ == "__main__":
	from optparse import OptionParser
	import sys, genome_viability
	
	parser = OptionParser()
	parser.add_option('-f', '--folder', action = 'store', type = 'string', dest = 'inputfile', help = 'folder with annotated genomes')
	parser.add_option('-e', '--essential', action = 'store', type = 'string', dest = 'essentials', help = 'list of essential genes')
	parser.add_option('-s', '--synthetic', action = 'store', type = 'string', dest = 'synth', help = 'synthetic lethal dataset (constanzo, 2010; lenient cutoff)')
	parser.add_option('-m', '--matrix', action = 'store', type = 'string', dest = 'synth', help = 'substitution matrix for determining mutation effect (default=BLOSUM80')
	(option, args) = parser.parse_args()

	viable = main(option.inputfile, option.essentials, option.synth)

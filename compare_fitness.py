#/usr/bin/python

"""
compare_fitness.py

creates list of shared variants in two Enrich ratios files, outputs paired fitness values for each

python compare_fitness.py -f1 FILE1 --f2 FILE2 -o OUTPUTFILENAME

"""

def main():
	#dictionaries to store fitness data (readID: fitness)
	f1_dict = {}
	f2_dict = {}
	#read fitness values
	for line in open(option.f1, 'r'):
		if line.startswith('seqID') != True:
			l = line.strip().split('\t')
			f1_dict[l[0]] = l[7]
	for line in open(option.f2, 'r'):
		if line.startswith('seqID') != True:
			l = line.strip().split('\t')
			f2_dict[l[0]] = l[7]
	
		
	f_out_name = option.f1.split('ratios_')[0] + option.f1.split('ratios_')[1].split('_qc_')[0] + '-' + option.f2.split('ratios_')[1].split('_qc_')[0] + '_fitness'
	if ".m1" in f1_dict:
		f_out_name += '.m1'
	if ".m2" in f1_dict:
		f_out_name += '.m2'
	f_out = open(f_out_name, 'w')
	
	
	#print header line
	print >> f_out, '\t'.join(['seqID', option.f1.split('ratios_')[1].split('_qc_')[0], option.f2.split('ratios_')[1].split('_qc_')[0]])
	#print shared variant fitnesses
	for i in f1_dict:
		if i in f2_dict:
			print >> f_out, '\t'.join([i, f1_dict[i], f2_dict[i]])
	f_out.close()
	
if __name__ == '__main__':

	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('--f1', action = 'store', type = 'string', dest = 'f1', help = 'path to file1')
	parser.add_option('--f2', action = 'store', type = 'string', dest = 'f2', help = 'path to file2')
	parser.add_option('-o', action = 'store', type = 'string', dest = 'f_out', help = 'output file')
	(option, args) = parser.parse_args()

	main()
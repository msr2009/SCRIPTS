"""
annotation_summary.py

creates a summary of the annotations for a given BED file.
"""

def main(f_in, mut_types = ''):
	f_out_handle = '.'.join(f_in.strip().split('.')[:-2]) + '.summary'
	f_out = open(f_out_handle, 'w')
	
	m = set(mut_types.split(','))
	mut_outfiles = {}
	if len(m) != 0:
		for i in m:
			mut_outfiles[i] = open('.'.join(f_in.strip().split('.')[:-2]) + '.' + i, 'w')
	
	annotations = {}
	count_chrpos = set()
	
	for line in open(f_in,'r'):
		if line.startswith('chr\t') == True:
			continue
		l = line.strip().split('\t')
		#count annotations
		if l[5] not in annotations:
			annotations[l[5]] = 1
		else:
			annotations[l[5]] += 1
		#write information to mutations_type file if type matches
		if l[5] in m:
			print >> mut_outfiles[l[5]], '\t'.join(l[6:8])	
		
		count_chrpos.add('-'.join(l[0:2]))	

	for a in sorted(annotations.keys()):
		print >> f_out, '\t'.join([a, str(annotations[a])])
	print >> f_out, '\t'.join(['total', str(len(count_chrpos))])
	
	f_out.close()
	for i in mut_outfiles:
		mut_outfiles[i].close()

if __name__ == '__main__':
			
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'f_input', help = 'annotations file')
	parser.add_option('-m', '--mutation_type', action = 'store', type = 'string', dest = 'mutation_type', help = 'comma-delimited list of mutation types to provide more information for')
	(option, args) = parser.parse_args()
	
	mut_types = ''
	if option.mutation_type != None:
		mut_types = option.mutation_type
	
	main(option.f_input, mut_types)

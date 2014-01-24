"""
combine_enrich_replicates.py

script takes multiple enrich counts or ratios files, creates a single file containing all variants and counts

Matt Rich
6/2013
"""

def main(filelist, pad, header):

	all_data = []
	column = 0
	all_keys = set()
	
	for f in filelist:
		if f.split('/')[-1].startswith('counts'):
			column = 8
		else:
			column = 7
		
		#make a dictionary to store data
		f_data = {}
		
		#open file, add variant, data to dict
		for line in open(f, 'r'):
			if line.startswith('seqID'):
				continue
			l = line.strip().split('\t')
			f_data[l[0]] = l[column]
		
		#add f_data keys to list of all keys
		all_keys.update(f_data.keys())
		
		#add f_data to all_data
		all_data.append(f_data)
	
	#after gathering all data, print out each variant, padding with zeroes where needed
	#print header first
	print '\t'.join(["seqID"]+header)
	
	for v in all_keys:
		outline = [v]			
		for f in range(len(filelist)):
			if v in all_data[f]:
				outline.append(all_data[f][v])
			else:
				if pad == True:
					outline.append("0")
				else:
					outline.append("NA")
		
		if "NA" not in outline:
			print '\t'.join(outline)
			
		
if __name__ == '__main__':
	from optparse import OptionParser
				
	parser = OptionParser()
	parser.add_option('-f', '--files', action = 'store', type = 'string', dest = 'files', help = 'comma-delimited list of counts/ratios files')
	parser.add_option('--pad_zeroes', action = 'store_true', dest = 'pad', help = 'should I pad missing values with zeroes? default removes missing variants', default=False)
	parser.add_option('--header', action = 'store', dest = 'header', type = 'string', help = 'comma-delimited list of names to use as header. Must correspond to --files.')
	(option, args) = parser.parse_args()
	
	if option.header != None:
		main(option.files.split(','), option.pad, option.header.split(','))
	else:
		main(option.files.split(','), option.pad, option.files.split(','))
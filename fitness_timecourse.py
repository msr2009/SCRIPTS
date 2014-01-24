#/usr/bin/python

"""
mutation_timecourse.py

takes input of multiple ratios files, creates output of file with columns for each variant's fitness for each ratio file



"""

def main():

	data = []
	
	#split infiles list
	infile = option.infiles.split(',')
	#open each file, create dictionary for each variant/fitness
	for f in infile:
		d = {}
		for line in open(f, 'r'):
			l = line.strip().split('\t')
			d[l[0]] = l[7]
		data.append(d)
	#print header line
	print '\t'.join( ['variant', 'input'] + map(str, range(1, len(infile)+1)) )

	#make new dictionary combining data:
	#for each variant - missing data == "NA"

	combined_data = {}
	for i in range(len(data)):
		for v in data[i]:
			#check if v is already in combined data
			if v in combined_data:
				#if it is, check for missing data
				if len(data[i][v]) == i+1:	
					combined_data[v].append(data[i][v])
				else:
					for j in range(i-len(combined_data[v])+1):
						combined_data[v].append('NA')
					combined_data[v].append(data[i][v])
			else:
				combined_data[v] = [0]
				for j in range(len(combined_data[v])-1):
					combined_data[v].append('NA')
				combined_data[v].append(data[i][v])	

	for v in combined_data:
		if len(combined_data[v]) != len(data)+1:
			for i in range(len(data)+1-len(combined_data[v])):
				combined_data[v].append('NA')
		print '\t'.join( map(lambda x: str(x), [v] + combined_data[v] ) )

if __name__ == '__main__':

	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infiles', help = 'comma-delimited list of ratios files')
	parser.add_option('-m', '--include_missing', action = 'store_true', dest = 'missing', help = 'include missing data in output?', default = False)
	(option, args) = parser.parse_args()

	main()

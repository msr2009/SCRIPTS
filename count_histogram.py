"""
count_histogram.py

creates raw histogram data for quality scores. can calculate for multiple columns of file at once, outputting into a single file. 
Also, can filter based on the maximum number of mutations per read.


"""

def main():
	#dictionary to store hists
	hists = []
	cols = []
	colnames = []
	firstline = 1

	#if columns to counts inputted
	if option.ncol != None:
		cols = map(lambda x: int(x), option.ncol.split(','))
		for i in cols:
			hists.append({})
	
	#open file, read lines
	for line in open(option.infile, 'r'):
		line = line.strip().split()
		#skip wildtype lines
		if '100' in line:
			continue
		#skip lines that have more mutations than maxmut
		if True in set([len(x.split(',')) > option.maxmut for x in line]):
			continue
		
		#on the first line, create appropriate number of hist dictionaries in [hists] if not inputted already
		if firstline == 1:
			#if len(cols) == 0 because ncol not inputted
			if len(cols) == 0:
				for i in range(len(line)):
					cols.append(i)
					colnames.append(line[i])
					hists.append[{}]
			else:
				for i in cols:
					colnames.append(line[i])
			firstline = 0
		
		#count data from each line
		elif firstline == 0:
			for c in range(len(cols)): 
				#count each column depending on mode
				tmp_dat = map(lambda x: int(float(x)), line[cols[c]].split(','))
				if option.mode == 'all':
					for score in tmp_dat:	
						if score in hists[c]:
							hists[c][score] += 1
						else:
							hists[c][score] = 1
				elif option.mode == 'min':
					if min(tmp_dat) in hists[c]:
						hists[c][min(tmp_dat)] += 1
					else: 
						hists[c][min(tmp_dat)] = 1
				elif option.mode == 'max':
					if max(tmp_dat) in hists[c]:
						hists[c][max(tmp_dat)] += 1
					else:
						hists[c][max(tmp_dat)] = 1
		
	#print output header
	f_out = open(option.infile+'-histogram-' + option.mode, 'w')
	print >> f_out, 'bin' + '\t' + '\t'.join(colnames)
	#find max score
	ma = 0
	for i in range(len(hists)):
		if max(hists[i]) > ma:
			ma = max(hists[i])	
	#find min score
	mi = 100
	for i in range(len(hists)): 
		if min(hists[i]) < mi:
			mi = min(hists[i])
	
	#print histograms
	for i in range(mi,ma+1):
		tmp_line = [str(i)]
		for j in range(len(hists)):
			if i in hists[j]:
				tmp_line.append(str(hists[j][i]))
			else:
				tmp_line.append('0')
		print >> f_out, '\t'.join(tmp_line)
	f_out.close()
		
if __name__ == '__main__':

	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'path to file')
	parser.add_option('-m', '--mode', action = 'store', type = 'string', dest = 'mode', help = 'output mode: all,min,max')
	parser.add_option('-c', '--columns', action = 'store', type = 'string', dest = 'ncol', help = 'comma-delimited list of which columns to count (default=all in file)')
	parser.add_option('--maxmut', action = 'store', type = 'int', dest = 'maxmut', help = 'max number of mutations in each read to use in counting')
	(option, args) = parser.parse_args()

	main()
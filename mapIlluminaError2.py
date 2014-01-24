#!/usr/bin/env python
#This script outputs a histogram of Sanger-scaled quality scores at each position from a raw Illumina .fq file

import optparse, time, pdb, sys, time, numpy, os

# define a function to adjust Illumina Shendure pipeline scores Sanger-scaled (i.e. ASCII 33-126) scores, see DF-127
def adjustscore(x):
	return x-33
	
def main():
	
	print time.asctime(time.localtime())
	
	parser = optparse.OptionParser()
	parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path from script to the data directory')
	parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'Illumina (.fq) file')
	(option, args) = parser.parse_args()
	
	#open the files for input and output
	f_infile = open((option.path + 'raw/' + option.infile), 'U')		
	f_output = open((option.path + 'raw/' + option.infile.rstrip('.fq') + '_quality_histogram'), 'w')
	
	print >> f_output, '\t'.join(['position', 'quality', 'counts'])
	
	# initialize some variables
	lcount = 0
	read_length = 0
	single_line = True
	
	read_quality = []
	#initialize the dictionary to hold the histogram
	data = dict()

	'''
	for keys in data:
		print keys
		print data[keys]
	'''
	
	while True:
		line = f_infile.readline().rstrip()

		#check to see if EOF has arrived
		if len(line) == 0:
			for k in data:
				for j in data[k]:
					print >> f_output, '\t'.join(map(str, [k, j, data[k][j]]))			
			f_output.close()
			
			#for keys in data:
				#print '\t'.join(map(str, data[k].values()))
			#os.system('rscript ../R/Hist_plotter.r ' + option.path + ' ' + option.infile + '_histogram_col' + str(option.colnum) + '_output')
			
			print 'GRACEFUL EXIT: EOF'
			print time.asctime(time.localtime())
			break
			
		#find the max number of mutations in a run
		elif lcount % 4 == 3:
			read_quality = map(adjustscore, map(ord, line))
			#print read_quality
			#print len(read_quality)
			for i in range(read_length):
				data[i][read_quality[i]] += 1				
			
		elif lcount % 4 == 1: 
			if single_line == True:
				read_length = len(line)
				# create dictionary to hold data
				for i in range(read_length):
					data[i] = dict() 
					for j in range(93):
						data[i][j] = 0
				# flip single_line so this doesn't happen again
				single_line = False
		
		lcount += 1
							
if __name__ == '__main__':
	        	main()
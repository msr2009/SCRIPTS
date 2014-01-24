"""
parse_Biotek_growthcurves.py

Script parses the bullshit output from our Biotek plate reader, getting rid of all the extraneous 
bullshit information and changes the times to something reasonable, instead of bullshit HH:MM:SS.

Matt Rich 04/2013
"""

def main(bullshit_infile, mapping_file, combine_replicates, theonlydatawewant, smooth, out_rates, dt):
	
	curves = {}	
	
	#open the shitty file and make it less shitty:
	#find the data that we want (which isn't the only shit in this file)
	stupid_transposed_data = [x.split('\t') for x in "".join(open(bullshit_infile,'r').read().split(theonlydatawewant)[1:3]).split("Results")[0].strip().split('\r\n')]
	#transpose it so that all the data for each well is on the same row (not sure why this isn't the case anyway)
	better_data = zip(*stupid_transposed_data)

	#open the mapping file, and store that information
	mapping = {}
	for line in open(mapping_file,'r'):
		
		if line.startswith('Row') == True:
			continue
			
		l = line.strip().split('\t')
		mapping[l[0]+l[1]] = l[2] + ":" + l[3]

	#loop through each curve, and store those data together, based on their well mapping
	for c in better_data:
		if c[0] == "Time":
			curves["Time"] = [ convertHHMMSStoDecimal(x) for x in c[1:] ] 
		elif c[0] == "T\xb0 ":
			continue
		else:
			strain_cond = mapping[c[0]]
			data = list([float(x) for x in c[1:]])
			if smooth == True:
				data = smoothData(data)
				
			if strain_cond in curves:
				curves[strain_cond].append(data)
			else:
				curves[strain_cond] = [data]
	
	#figure out what the time interval is, convert to hours
	if dt == 0.0:
		intervals = set( [ float(curves["Time"][i+1]) - float(curves["Time"][i]) for i in range(len(curves["Time"])-1) ] )
		if len(intervals) == 1:
			dt = intervals.pop()/60.0
		else:
			print "Time intervals are not consistent, setting dt to 30min (0.5hr)"
			dt = 0.5
	
	#open file to output rates, if out_rates == true
	if out_rates == True:
		if combine_replicates == True:
			f_rates = open(bullshit_infile.rstrip(".txt")+".rates-means.txt",'w')	
		else:
			f_rates = open(bullshit_infile.rstrip(".txt")+".rates.txt",'w')	
		
	#open output file(s) and write output
	if combine_replicates == False:
		#don't do any math to combine replicates
		f_out = open(bullshit_infile.rstrip(".txt")+".growthcurves.txt",'w')
		#print header
		print >> f_out, '\t'.join(['Strain', 'Condition'] + curves["Time"])
		#print data
		for c in curves:
			if c != "Time":
				sc = c.split(':')
				for d in curves[c]:
					if out_rates == True:
						print >> f_rates, '\t'.join(sc + [str(calculateGrowthRate(d, dt))])
					print >> f_out, '\t'.join(sc + [str(x) for x in d])
		f_out.close()
	
	else:
		#do the maths
		f_avg = open(bullshit_infile.rstrip(".txt")+".growthcurves-means.txt",'w')
		f_sd = open(bullshit_infile.rstrip(".txt")+".growthcurves-stdev.txt",'w')
		
		#print header
		print >> f_avg, '\t'.join(['Strain', 'Condition'] + curves["Time"])
		print >> f_sd, '\t'.join(['Strain', 'Condition'] + curves["Time"])

		#print data
		for c in curves:
			if c != "Time":
				sc = c.split(':')
				#calculate mean, std using numpy arrays
				print >> f_avg, '\t'.join( sc + [str(x) for x in list(np.mean(np.array(curves[c]), axis=0))] )
				print >> f_sd, '\t'.join( sc + [str(x) for x in list(np.std(np.array(curves[c]), axis=0))] )
				if out_rates == True:
					print >> f_rates, '\t'.join(sc + [str(calculateGrowthRate(list(np.mean(np.array(curves[c]), axis=0)), dt))])
		f_avg.close()
		f_sd.close()
		
def convertHHMMSStoDecimal(HHMMSS):
	[H,M,S] = [float(x) for x in HHMMSS.split(":")]
	return str(H*60 + M + S/60)

def smoothData(vector):
	a = vector[0:2]
	b = vector[-2:]
	sv = [ np.mean(vector[i-2:i+3]) for i in range(len(vector))[2:-2] ]
	return a + sv + b 
	
def calculateGrowthRate(vector, dt):
	vector = [ log(x) for x in vector ]
	#calculate slopes 
	slopes = [ (vector[x+1]-vector[x])/dt for x in range(8,len(vector)-1) ]
	#sort values in slopes, then average the 3rd-8th highest
	return sum(sorted(slopes)[-8:-2])/6.0
	
	
if __name__ == '__main__':
	from optparse import OptionParser
	import numpy as np
	from math import log
	
	parser = OptionParser()
	parser.add_option('-f', '--bullshit', action = 'store', type = 'string', dest = 'infile', help = 'file containing data, in bullshit format')
	parser.add_option('-m', '--mapping', action = 'store', type = 'string', dest = 'mapping', help = 'file that contains information mapping wells to samples (row,column,sample,condition)')
	parser.add_option('--combine_replicates', action = 'store_true', dest = 'combine', help = 'should the script combine (outputting mean and SD) replicates with same name and condition?', default=False)
	parser.add_option('-d', '--data_name', action = 'store', type = 'string', dest = 'data_name', help = 'name of data (default = "Read 2:660")', default="Read 2:660")
	parser.add_option('--smooth', action = 'store_true', dest = 'smooth', help = 'should the data be smoothed? (each point becomes average of 7 points)', default=False)
	parser.add_option('--rates', action = 'store_true', dest = 'rates', help = 'output file containing growth rates?', default=False)
	parser.add_option('--dt', action = 'store', type = 'float', dest = 'deltaT', help = 'time between measurements (in hours)', default=0)	

	(option, args) = parser.parse_args()
	
	main(option.infile, option.mapping, option.combine, option.data_name, option.smooth, option.rates, option.deltaT)	

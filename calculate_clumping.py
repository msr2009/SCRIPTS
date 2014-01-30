"""
calculate_clumping.py

calculates clumping statistic from CellProfiler output
this is essentially just a copy of calculate_invasion.py, but for clumping 

Matt Rich 1/2014
"""

def main(num, exc, decile, noprint):
	outer = {}
	
	#read numerator file
	for line in open(num, "rU").readlines()[1:]:
		l = line.strip().split('\t')
		if l[2] + l[3] not in exc:
			outer[":".join([l[0],l[2],l[3]])] = [l[6],l[9]]
	
	#print header
	if noprint != True:
		print "\t".join(["plate", "row", "column", "score", "mean_intensity"]) 

	#do maths, print output
	wells_dat = {}
	dat_list = []
	
	for well in outer:
		if outer[well][0] == "nan":
			if noprint != True:	
				print "\t".join(well.split(":") + ["NA", "NA"])
		else:
			if noprint != True:
				print "\t".join(well.split(":") + [str(x) for x in [outer[well][1], outer[well][0]]])
			wells_dat[well] = outer[well][1]	
			dat_list.append(outer[well][1])
	
	if decile != None:
		#open some output files
		deciles = [10,20,30,40,50,60,70,80,90,100]
		out_files = {}
		for d in deciles:
			out_files[d] = open(decile+"."+str(d)+".txt", "w")
			#print header
			print >> out_files[d], "\t".join(["plate", "row", "column", "score", "percentage"])
		#loop through the scores, outputting them to the appropriate files
		for w in wells_dat:
			p = percentileofscore(dat_list, wells_dat[w])
			#print w, wells_dat[w], p
			
			if p < deciles[0]:
				print >> out_files[deciles[0]], "\t".join(w.split(":") + [str(wells_dat[w]), str(p)])
			else:
				for x in range(len(deciles)-1):
					if p > deciles[x] and p <= deciles[x+1]:
						print >> out_files[deciles[x+1]], "\t".join(w.split(":") + [str(wells_dat[w]), str(p)])	
						break
		
		#close all out_files
		for f in out_files:
			out_files[f].close()
		
	
if __name__ == "__main__":
	from optparse import OptionParser
	from scipy.stats import percentileofscore
	
	parser = OptionParser()
	parser.add_option("-d", "--data", type = "string", action = "store", dest = "data", help = "CellProfiler data output file")
	parser.add_option("--exclude_wells", type = "string", action = "store", dest = "exclude", help = "comma-delimited list of wells to ignore from every plate", default = "")
	parser.add_option("--output_deciles", type = "string", action = "store", dest = "deciles", help = "file prefix for decile output", default = None)
	parser.add_option("--no_print", action = "store_true", dest = "noprint", help = "don't print output", default = False)
	(option, args) = parser.parse_args()

	main(option.data, set(option.exclude.split(",")), option.deciles, option.noprint)

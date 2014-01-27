"""
calculate_invasion.py

calculates invasion statistic from CellProfiler output

Matt Rich 1/2014
"""

def main(num, den, exc, decile, noprint):
	inner = {}
	outer = {}
	
	#read numerator file
	for line in open(num, "rU").readlines()[1:]:
		l = line.strip().split('\t')
		if l[2] + l[3] not in exc:
			inner[":".join([l[0],l[2],l[3]])] = l[4:]
	
	#read denominator file
	for line in open(den, "rU").readlines()[1:]:
		l = line.strip().split('\t')
		if l[2] + l[3] not in exc:
			outer[":".join([l[0],l[2],l[3]])] = l[4:]
	
	#print header
	if noprint != True:
		print "\t".join(["plate", "row", "column", "score", "pre_background", "post_background", "pre_intensity", "post_intensity"]) 

	#do maths, print output
	wells_dat = {}
	dat_list = []
	
	for well in outer:
		if outer[well][0] == "nan":
			if noprint != True:	
				print "\t".join(well.split(":") + ["NA", "NA", "NA", "NA", "NA"])
		else:
			pre_background = float(outer[well][1])
			post_background = float(outer[well][0])
			pre_intensity = float(inner[well][3])
			post_intensity = float(inner[well][2])
			dat = (post_intensity-post_background)/(pre_intensity-pre_background)
			if noprint != True:
				print "\t".join(well.split(":") + [str(x) for x in [dat, pre_background, post_background, pre_intensity, post_intensity]])
			wells_dat[well] = dat	
			dat_list.append(dat)
	print len(dat_list)
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
	parser.add_option("-i", "--inner", type = "string", action = "store", dest = "inner", help = "inner circle intensities")
	parser.add_option("-o", "--outer", type = "string", action = "store", dest = "outer", help = "outer circle intensities")
	parser.add_option("--exclude_wells", type = "string", action = "store", dest = "exclude", help = "comma-delimited list of wells to ignore from every plate", default = "")
	parser.add_option("--output_deciles", type = "string", action = "store", dest = "deciles", help = "file prefix for decile output", default = None)
	parser.add_option("--no_print", action = "store_true", dest = "noprint", help = "don't print output", default = False)
	(option, args) = parser.parse_args()

	main(option.inner, option.outer, set(option.exclude.split(",")), option.deciles, option.noprint)

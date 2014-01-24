"""
parse_Biotek_plate.py

Matt Rich 01/2014
"""

def main(m, p, all_out, compare_strains):
	#first, make mapping dictionary
	plate_map = {}
	for line in open(m, "r"):
		l = line.strip().split('\t')
		if line.startswith("#") == True:
			#print "ignoring" + l[0]+l[1]+l[2]
			continue
		plate_map[l[0] + l[1] + l[2]] = l[3] + ":" + l[4]
	
#	for x in sorted(plate_map.keys()):
#		print x, plate_map[x]	
	
	#some counters
	plate_counter = 1
	row_counter = 0
	rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
	
	data = {}
	
	#read plate files one at a time
	for plate in p:
		
		#then read rows one at a time
		for line in open(plate, "r"):
			
			r = line.strip().split("\t")	
			#read data from wells in each row
			for col in range(len(r)):
				
				well_name = str(plate_counter)+rows[row_counter]+str(col+1)
				
				#add to data if already seen once
				try:
					if plate_map[well_name] in data:
						data[plate_map[well_name]].append(float(r[col]))
					else:
						data[plate_map[well_name]] = [float(r[col])]
				except KeyError:
					pass

			row_counter += 1 
		plate_counter += 1
		row_counter = 0
	
#	for x in data:
#		print x, data[x]
	
	#print all data to file
	f_out = open(all_out, "w")
	#print header
	print "\t".join(["strain", "mean", "sd"] + compare_strains)
	for strain in data:
		for x in data[strain]:
			print >> f_out, "\t".join(strain.split(":") + [str(x)])
		#do some maths
		print "\t".join([ strain, str(mean(data[strain])), str(std(data[strain])) ] + [ str(ttest_ind(data[strain], data[comp])[1]) for comp in compare_strains ] )
		
	f_out.close()		

if __name__ == "__main__":
	from optparse import OptionParser
  	from scipy.stats import ttest_ind
	from numpy import mean, std	
	parser = OptionParser()
	parser.add_option("-m", "--mapping", action = "store", type = "string", dest = "mapping", help = "tab-delimited mapping file for plates")
	parser.add_option("-p", "--plates", action = "store", type = "string", dest = "plates", help = "comma-delimited list of plates")
	parser.add_option("--output_all", action = "store", type = "string", dest = "allout", help = "print all data, you know, for boxplots and stuff")
	parser.add_option("--compare_strains", action = "store", type = "string", dest = "compstrains", help = "strains to calculate t-test against", default="")
	(option, args) = parser.parse_args()

	main(option.mapping, option.plates.split(","), option.allout, option.compstrains.split(","))
    

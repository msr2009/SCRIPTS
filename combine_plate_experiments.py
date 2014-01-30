"""
combine_plate_experiments.py

combines data from 96-well plate-based experiments

Matt Rich 1/2014
"""

def main(infiles):
	#store data in dat
	dat = {}
	header = []
	previous = 0
	#loop through list of files containing data
	for f in infiles:	
		for line in open(f, "rU"):
			#this is the header
			if line.startswith("plate"):
				if "plate" in dat:
					dat["plate"] = dat["plate"] + line.strip().split("\t")[3:]
				else:
					dat["plate"] = line.strip().split("\t")
				continue
			
			l = line.strip().split("\t")
			well = ":".join(l[0:3])
			if well in dat:
				dat[well] = dat[well]+l[3:]
			else:
				dat[well] = ["NA"]*previous + l[3:]

		#keep track of how many columns are in this file
		previous = len(l[3:])	
		
	#print header
	print "\t".join(dat["plate"])
	dat.pop("plate")

	#print dat
	for d in dat:
		print "\t".join(d.split(":") + dat[d])

if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option("-f", "--filelist", action = "store", type = "string", dest = "filelist", help = "comma-delimited list of files to be combined by well")
	(option, args) = parser.parse_args()

	main(option.filelist.split(","))

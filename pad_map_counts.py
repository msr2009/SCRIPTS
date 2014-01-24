"""
pad_map_counts.py

Takes input and selected map_counts output files, and pads all variants with 1 count, creating variants as necessary in selected file.

"""

def main(infile, selfile):
	
	indict, in_counts = makeCountsDict(infile)
	seldict, sel_counts = makeCountsDict(selfile)
	
	#open files for output
	in_out = open(infile+".padded", "w")
	sel_out = open(selfile+".padded", "w")
	#write headers to outfiles
	print >> in_out, "\t".join(indict["seqID"])
	print >> sel_out, "\t".join(seldict["seqID"])
	#delete headers from dicts
	del indict["seqID"]
	del seldict["seqID"]
	#calculate total number of counts after padding
	in_counts += len(indict)
	#sel_counts gets added 1 for each variant in seldict, plus 1 for each variant in indict that ISN'T in seldict
	sel_counts += len(seldict) + len(set(indict.keys()).difference(set(seldict.keys())))
	
	#pad all variants with 1 count
	for variant in indict:
		#pad input variant with a count, print to in_out
		print >> in_out, "\t".join(indict[variant][:-2] + [str((int(indict[variant][-1])+1)/in_counts), str(int(indict[variant][-1])+1)])
			
		#if variant in selected library, pad it and print to sel_out
		if variant in seldict:
			print >> sel_out, "\t".join(seldict[variant][:-2] + [str((int(seldict[variant][-1])+1)/sel_counts), str(int(seldict[variant][-1])+1)])
		#if it isn't, make a new entry in sel_out with one count
		else:	
			print >> sel_out, "\t".join(indict[variant][:-2] + [str(1/sel_counts), str(1)]) 

	in_out.close()
	sel_out.close()
	
def makeCountsDict(f):
	d = {}
	total_counts = 0
	for line in open(f, "r"):
		l = line.strip().split("\t")
		d[l[0]] = l
		if l[0] != "seqID":
			total_counts += int(l[-1])
	return d, float(total_counts)

if __name__ == '__main__':
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option("--input", action = "store", type = "string", dest = "inputfile", help = "input counts file")
	parser.add_option("--selected", action = "store", type = "string", dest = "selectedfile", help = "selected counts file")
	(option, args) = parser.parse_args()
	
	main(option.inputfile, option.selectedfile)	
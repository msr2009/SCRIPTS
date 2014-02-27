"""
parseVCF.py

Matt Rich 2/2014

"""

def main(vcf):
	for line in open(vcf, "rU"):
		#skip header lines
		if line.startswith("#"): continue	
		#otherwise, make dict out of data
		l = line.strip().split("\t")
		out = [l[1], l[3], l[4]]
		dat = { key: value for (key, value) in [x.split("=") for x in l[7].split(";")] }
		dp4 = dat["DP4"].split(",")	
		out.append( int(dp4[0])+int(dp4[1]) )
		out.append( int( dp4[2])+int(dp4[3]) )
		print "\t".join( [str(x) for x in out] )

if __name__ == "__main__":
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-v", "--vcf", action = "store", type = "string", dest = "vcf", help = "VCF file containing mutations")
	(option, args) = parser.parse_args()

	main(option.vcf)	

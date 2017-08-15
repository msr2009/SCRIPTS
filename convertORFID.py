"""
convertORFID.py

a quickie to add a column to a file with the standard name of an ORF to a file

Matt Rich, 1/2015
"""

def main(fin, fconv, col, sep, header):
	
	#read in conversion file
	c = {}
	for line in open(fconv, "r"):
		l = line.strip().split("\t")
		if len(l) != 1:
			c[l[0]] = l[1]
		else:
			c[l[0]] = l[0]

	#read file for conversion, add column to end with standard orf name	
	for line in open(fin, "r"):
		l = line.strip().split(sep)
		#if there's a header add the column name there first
		if header:
			print sep.join(l + ["standardORF"])
			header = False
			continue
		#otherwise, do the conversion and print out the standard name
		try:
			print sep.join(l + [ c[l[col]] ])
		except KeyError:
			print sep.join(l + [ l[col] ]) 

if __name__ == "__main__":
	
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'fin', help = "file for conversion")
	parser.add_option('--conversion', action = 'store', type = 'string', dest = 'conv', help = "file with standard and systematic orfs", default = "/net/fields/vol2/home/mattrich/YEAST/orfID-genename.txt")
	parser.add_option('--sep', action = 'store', type = 'string', dest = 'sep', help = "separator", default = "\t")
	parser.add_option('--col', action = 'store', type = 'int', dest = 'col', help = "column containing systematic orf name", default = 0)
	parser.add_option('--no-header', action = 'store_false', dest = "header", help = "file does not have a header", default = True)

	(option, args) = parser.parse_args()
	
	main(option.fin, option.conv, option.col, option.sep, option.header)	


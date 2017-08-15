"""
meme2tab.py

creates single tab file for every motif found in a meme file

Matt Rich, 10/2014
"""

def main(meme):

	f = open(meme,"r").read()
#	r = re.compile("MOTIF (.*)\n\n.*\n?\n?letter-probability matrix: .*\n([\.\d\W\n]*)[URL\w\./:]*\n")
	r = re.compile("MOTIF (.*)\n")
	p = re.compile("letter-probability matrix: .*\n([\.\d\W\n]*)\n")
	motifs = r.findall(f)
	pwms = p.findall(f)

	for x in range(len(motifs)):
		print ">" + ".".join(motifs[x].lstrip("MOTIF ").split()) + "." + str(x)
		dat = [ y.strip().split() for y in pwms[x].split("\n") ]
		for y in dat:
			print "\t".join([str(round(float(j)*1000)) for j in y])		

if __name__ == "__main__":
	
	from optparse import OptionParser
	import re


	parser = OptionParser()
	parser.add_option('--meme', action = 'store', type = 'string', dest = 'meme', help = "meme file to parse")
	(option, args) = parser.parse_args()
	
	main(option.meme)	


"""
parseGFF.py 

parses GFF file to get information needed for yeast_annotation.py
"""

#parse GFF file into genes and noncoding annotated elements
def parseGFF(gff, exclude):
	#dictionaries to store genes or noncoding elements
	genes = {}
	noncoding = {}
	counter = 0

	#read GFF file, parse information
	for line in open(gff, 'r'):
		#skip header lines
		if line.startswith('#') == True:
			continue
		elif line.startswith('>') == True:
			break
			
		l = line.strip().split('\t')
		
		#check if annotation is in exclude list, and skip if it is
		if l[2] in exclude:
			continue
		
		#list to hold annotation information found in l[8]
		ann_list = {y[0]: y[1] for y in [x.split('=') for x in l[8].split(';')]}

		#get chromosome
		ch = l[0][3:]
		#check annotation
		#if it's a gene, then add the start, stop, and chromosome to the genes dictionary
		if l[2] == "gene":
			genes[ann_list['ID']] = [ann_list['ID'], int(l[3]), int(l[4]), [  ], [  ], 0, chromosome_conversion(l[0][3:]), l[6] ]		
		#if there's a CDS or intron associated, then add those as well (CDS=exon)	
		elif l[2] == "CDS" and ann_list['Parent'] in genes:
			genes[ann_list['Parent']][3].append([int(l[3]),int(l[4])])
		elif l[2] == "intron" and ann_list['Parent'] in genes:
			genes[ann_list['Parent']][4].append([int(l[3]),int(l[4])])
		
		#otherwise, it'll be a noncoding region (or tRNA, etc...)
		else:
			noncoding[counter] = [l[8].split(';')[0].split('=')[1], chromosome_conversion(l[0][3:]), l[2], int(l[3]), int(l[4])]
			counter += 1
	
	#before returning, we have to flip the positions for any gene that's on the Crick strand
	for g in genes:
		if genes[g][-1] == "-":
			genes[g] = [ genes[g][0], genes[g][2], genes[g][1], [ x[::-1] for x in genes[g][3] ], [ x[::-1] for x in genes[g][4] ], (sum([ abs(x-y) for x,y in genes[g][3] ])+1)/3.0, genes[g][6], genes[g][7] ]
		else:
			genes[g][5] = (sum([ abs(x-y) for x,y in genes[g][3] ])+1)/3.0
			
	return (genes, noncoding)	
				
def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':17, 'mitochondrion':17,
				'to':17, 'on':17, 'M':17}
	
	chrom_number = chrom_number.split('chr')[-1]
	
	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:			
		return chrom_conv[chrom_number]		 


if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option("-g", "--gff", action = 'store', dest = 'gff', type = 'string', help = "path to GFF file")
	parser.add_option("--exclude", action = 'store', dest = 'exclude', type = 'string', help = "comma-delimited list of annotations to exclude")
	(option, args) = parser.parse_args()
	
	exclude = set(['chromosome', 'repeat_region', 'region', 'nucleotide_match', 'long_terminal_repeat', 'noncoding_exon', 'pseudogene', 'telomere'])
	if option.exclude != None:
		exclude = set(option.exclude.split(','))

	parseGFF(option.gff, exclude)
	
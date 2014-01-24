"""
annotate_genomes.py

wrapper script performs annotation on single yeast genome

1) mutation_caller.py
2) yeast_annotation.py

can be run by only calling 'python annotate_genomes.py -m GENOMEFILE' if using default genome, annotation files

"""

def main(genome, wt, orf, nc, only):
	
	#call mutation_caller
	print "Calling mutations for " + genome
	mut_call = "python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/mutation_caller.py -r " + wt + " -o " + genome
	subprocess.call( mut_call, shell=True )
	
	#call yeast_annotation on mutation_caller output
	genome2 = '.'.join(genome.split('.')[0:-1])+'.mutations'
	print "Annotating mutations for " + genome2
	ann_call = "python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/yeast_annotation.py -f " + genome2 + " -s " + orf + " -n " + nc + " -g " + wt
	subprocess.call( ann_call, shell=True )
	
	#call annotation_summary on yeast_annotation output
	genome3 = '.'.join(genome.split('.')[0:-1])+'.mutations.annotated'
	print "Creating summary for " + genome3
	sum_call = "python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/yeast_annotation.py -f " + genome3
	subprocess.call( sum_call, shell=True )

if __name__ == '__main__':
			
	from optparse import OptionParser
	import sys, math, os, re, subprocess
	
	parser = OptionParser()
	parser.add_option('-m', '--mutant', action = 'store', type = 'string', dest = 'f_input', help = 'path to mutated genome fasta file')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'wtgenome', help = 'path to wildtype genome')
	parser.add_option('-o', '--coding', action = 'store', type = 'string', dest = 'orf', help = 'ORF annotation file')	
	parser.add_option('-n', '--non-coding', action = 'store', type = 'string', dest = 'noncoding', help = 'non-coding annotation file')
	parser.add_option('--annotate_only', action = 'store_true', dest = 'annotate_only', help = 'only annotate mutations (BED format only)', default=False)
	(option, args) = parser.parse_args()
	
	if option.wtgenome != None:
		wildtype = option.wtgenome
	else:
		wildtype = '/net/fields/vol1/home/mattrich/YEAST/S288C_reference_sequence_R64-1-1_20110203.fsa'
	
	if option.orf != None:
		orf_annotation = option.orf
	else:
		orf_annotation = '/net/fields/vol1/home/mattrich/YEAST/orf_coding_all_R64-1-1_20110203.fasta'
		
	if option.noncoding != None:
		nc_annotation = option.noncoding
	else:
		nc_annotation = '/net/fields/vol1/home/mattrich/YEAST/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered'
	
	main(option.f_input, wildtype, orf_annotation, nc_annotation, option.annotate_only)

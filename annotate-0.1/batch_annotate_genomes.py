"""
batch_annotate_genomes.py

wrapper script performs annotations on all genomes found in a folder

1) mutation_caller.py
2) yeast_annotation.py

can be run by only calling 'python batch_annotate_genomes.py -f GENOMEFILEFOLDER' if using default genome, annotation files

"""

def main(folder, wt, orf, nc, clus=False):
	
	#get list of genome files
	f = os.listdir(folder)
	f = map(lambda x: folder + x, f)
	#loop through list, performing annotations
	for i in f:
		if i.endswith('annotated') == False and i.endswith('mutations') == False:
			#call annotate_genomes.py
			if clus==False:
				print "=============================================================="
				print i
				annotate_call = "python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/annotate_genomes.py -m " + i
				subprocess.call( annotate_call, shell=True )
			else:
				print "calling mutations and annotating " + i + " on SGE"
				annotate_call = "qsub -V -S /bin/bash /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/batch_annotate_genomes_SGE.sh " + i
				subprocess.call( annotate_call, shell=True )


if __name__ == '__main__':
			
	from optparse import OptionParser
	import sys, math, os, re, subprocess
	
	parser = OptionParser()
	parser.add_option('-f', '--folder', action = 'store', type = 'string', dest = 'f_input', help = 'folder containing files to be annotated')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'wtgenome', help = 'path to wildtype genome')
	parser.add_option('-o', '--coding', action = 'store', type = 'string', dest = 'orf', help = 'ORF annotation file')	
	parser.add_option('-n', '--non-coding', action = 'store', type = 'string', dest = 'noncoding', help = 'non-coding annotation file')
	parser.add_option('--cluster', action = 'store_true', dest='cluster', help = 'run in parallel on SGE cluster?', default=False)
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
	
	main(option.f_input, wildtype, orf_annotation, nc_annotation, option.cluster)
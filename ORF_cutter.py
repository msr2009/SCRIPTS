"""
ORF_cutter.py
blah blah blah
Script that searches sequences for restriction enzyme sites

python ORF_cutter.py --orfs ORFS_FASTA --sites COMMA-DELIMITED LIST OF SEQUENCES
"""

from optparse import OptionParser
import sys, re
from Bio import SeqIO
from Bio.Restriction import *

def cutORFs(orfs_file, sites_list, mode):
		
	#create dictionary of ORF sequences
	genes = {}
	for record in SeqIO.parse( open(orfs_file, 'r'), 'fasta' ):
		genes[record.id.split(' ')[0]] = record.seq.tostring()
	
	cuts = {}	
	if mode == 'single':
		for i in sites_list.split(','):
			cuts[i.upper()] = set()
		#loop through ORFs, checking if they contain sites
		for orf in genes:
			for s in cuts:
				if s in genes[orf]:
					cuts[s].add(orf)
	
	elif mode == 'all':
		for i in sites_list.split(','):
			for j in sites_list.split(','):
				if j != i:
					cuts[i+','+j] = set()
		#loop through ORFs, checking if they contain sites for either enzyme
		for orf in genes:
			for s in cuts:
				if s.split(',')[0].upper() not in genes[orf] or s.split(',')[1].upper() not in genes[orf]:
					cuts[s].add(orf)
					
	#make list of all NEB restriction enzymes
	rb = RestrictionBatch(first=[], suppliers='N')
		
	for s in cuts:
	#see if sequences are known restriction enzymes
		enzyme = False
		for re in rb:
			if re.site == s:
				enzyme = re
				break
					
		if enzyme != False:
			print enzyme, s, len(cuts[s]), (len(genes)-len(cuts[s]))/float(len(genes)), ','.join(list(cuts[s]))
		else:
			print s, s, len(cuts[s]), (len(genes)-len(cuts[s]))/float(len(genes)), ','.join(list(cuts[s]))
			
if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option('-o', '--orfs', action = 'store', type = 'string', dest = 'orfs', help = 'fasta file of ORF sequences')
	parser.add_option('-s', '--sites', action = 'store', type = 'string', dest = 'sites', help = 'comma-delimited list of ')
	parser.add_option('--mode', action = 'store', type = 'string', dest = 'mode', help = 'mode for analyzing enzymes: single, pair, all')
	(option, args) = parser.parse_args()
	
	#check if sequences or enzyme names were given
	sites = []
	sites_RB = RestrictionBatch(first=[])
	for s in option.sites.split(','):
		try:
			sites_RB.add(eval(s))
		except NameError:
			sites.append(s)
		
	if len(sites_RB) != 0:
		cutORFs(option.orfs, ','.join([re.site for re in sites_RB]), option.mode)
	else:
		cutORFs(option.orfs, ','.join(sites), option.mode)

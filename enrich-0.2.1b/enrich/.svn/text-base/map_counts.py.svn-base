#!/usr/bin/env python
'''
map_counts: The map_counts module identifies unique variants in read_aligner output files and produces an output file that contains a list of those unique sequences as well as the number of times they appear normalized by the total number of counts in the input file (i.e. the frequency).
'''
###foo

import sys, time, optparse # import general libraries
import enrich_util # import project-specific libraries

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

def main(project_directory, input_file, qc_filename, dna_offset, grid = 'L', min_counts = 1):

	if grid != 'L': #print logging information, if being called from the cluster
		print time.asctime(time.localtime())
		print project_directory, input_file, min_counts
	
	if "_DNA_" in input_file:
		filter_chars = "N"
	elif "_PRO_" in input_file:
		filter_chars = "X"
	
	try:
		min_counts = int(min_counts)
	except:
		print 'Error: minimum counts argument is not integer.  Setting to 1.'
		min_counts = 1
	
	try:
		dna_offset = int(dna_offset)
	except:
		print 'Error: dna_offset cannot be cast as integer. Setting to 0.'
		dna_offset = 0
	
	try:
		#build a dictionary of the sequences observed and counts:
		tally_dict, properties_dict = enrich_util.build_tally_dict(project_directory + 'data/tmp/' + input_file)
		#a couple counters for stats:
		#where mutations are in sequence
		mutation_position = []
		#how many mutations in sequence
		mutation_number = []
		#total number of mutations for each mutation-count
		mutation_tally = []
		#variation in mutated base
		base_identity = []
		
		total_variants = 0
		
		#build a dictionary of the sequence codes and counts:
		count_dict = {}
		for mutation_location in sorted(tally_dict):
			for mutation_identity in sorted(tally_dict[mutation_location]):
				seqID = mutation_location + '-' + mutation_identity
				count_dict[seqID] = tally_dict[mutation_location][mutation_identity]
				if mutation_position == []:
					mutation_position = [ 0 for i in range(len(properties_dict[mutation_location][mutation_identity][0])) ]
					mutation_number = [ 0 for i in range(len(properties_dict[mutation_location][mutation_identity][0])) ]
					mutation_tally = [ 0 for i in range(len(properties_dict[mutation_location][mutation_identity][0])) ]
					if "_DNA_" in input_file:
						base_identity = [ {'A':0, 'C':0, 'T':0, 'G':0} for i in range(len(properties_dict[mutation_location][mutation_identity][0])) ]
		seqIDs = enrich_util.valuesort(count_dict)

		#build a dictionary of the normalized tally:
		norm_dict = enrich_util.norm_count_dict(count_dict)

	except:
		print 'Error: could not build dictionary of counts from input file'
		return(1)
		
	try:
		#print sequences and counts:
		f = open(project_directory + 'data/output/' + 'counts_' + input_file,'w')
		f_1 = open(project_directory + 'data/output/' + 'counts_' + input_file + '.m1', 'w')
		f_2 = open(project_directory + 'data/output/' + 'counts_' + input_file + '.m2', 'w')
	
		print >>f, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )
		print >>f_1, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )
		print >>f_2, '\t'.join(["seqID","sequence","match_count","mutation_count","mutation_location","mutation_identity","max_mutation_run","sequence_frequency","sequence_count"] )
		
		seqIDs.reverse()
		for seqID in seqIDs:
			mutation_location, mutation_identity = seqID.split('-')
			norm = norm_dict[seqID]
			tally = tally_dict[mutation_location][mutation_identity]
			sequence = properties_dict[mutation_location][mutation_identity][0]
			
			if not filter_chars in sequence and tally >= min_counts: #filter based on counts
				number_of_mutations = len(mutation_location.split(','))
				
				#update stats (only for DNA)
				if "_PRO_" not in input_file:
					if mutation_location != 'NA':
						mutation_number[number_of_mutations] += 1
						mutation_tally[number_of_mutations] += tally
						
						mut_locs = mutation_location.split(',')
						mut_ids = mutation_identity.split(',')
						
						for l in range(len(mut_locs)):
							loc = int(mut_locs[l])
							mutation_position[loc-dna_offset] += 1 
							base_identity[loc-dna_offset][mut_ids[l]] += 1	
					else:
						mutation_number[0] += 1
						mutation_tally[0] += tally
					 
				print >>f, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])
			 
				if number_of_mutations == 1 and mutation_location != 'NA':
					print >>f_1, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])	  
				
				if number_of_mutations == 2:
					print >>f_2, '\t'.join([seqID] + properties_dict[mutation_location][mutation_identity] + [str(norm), str(tally)])  
				
				total_variants += 1 
						 
		f.close()
		f_1.close()
		f_2.close()
		
		if "_PRO_" not in input_file:
			qc_data = {'mutation-positions_'+input_file:mutation_position, 'mutations-per-sequence_'+input_file:mutation_number, 'total-variants_'+input_file:total_variants, 'total-variants-per-mutation-number_'+input_file:mutation_tally, 'mutation-identity_'+input_file:base_identity}
			enrich_util.update_qc_file( project_directory + 'qc/' + qc_filename, qc_data)
		
	except IndexError:
		print 'Error: write to output file'
		return(1)

	return(0)
			
if __name__ == '__main__':
	parser = optparse.OptionParser()
	parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'exact location from / project directory')
	parser.add_option('--infile', action = 'store', type = 'string', dest = 'infile', help = 'input filename')
	parser.add_option('--qcfile', action = 'store', type = 'string', dest = 'qcfile', help = 'qc filename')
	parser.add_option('--dna_offset', action = 'store', type = 'string', dest = 'offset', help = 'dna offset')
	parser.add_option('--local', action = 'store', type = 'string', dest = 'local', help = 'Is this a local (L) run or should an SGE (SGE) job be scheduled?')
	(option, args) = parser.parse_args()	

	main(option.path, option.infile, option.qcfile, option.offset, option.local)
	

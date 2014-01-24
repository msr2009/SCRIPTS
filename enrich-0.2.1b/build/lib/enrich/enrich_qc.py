#!/usr/bin/python

__author__ = "Douglas M. Fowler"
__copyright__ = "Copyright 2011"
__credits__ = ["Douglas M Fowler", "Carlos L. Araya"]
__license__ = "FreeBSD"
__version__ = "0.2"
__maintainer__ = "Douglas M. Fowler"
__email__ = "dfowler@uw.edu"

import os, sys, matplotlib, json
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser

'''
enrich_qc: The enrich_qc module compiles quality control analysis from data that was compiled during a pipeline run.
'''

# to convert a qc file to a dictionary
def qc_json(qc_file_name, analysis_name):	# qc_file_name is 'config_file_name.qc' in qc directory
	qc_file = open(qc_file_name, 'r')	 
	qc_str = qc_file.read()
	qc_dic = json.loads(qc_str)

	analysis_dic = {}						
	for key in qc_dic.keys():				# keys are 'analysis_name+source_file'
		key_split = key.split('_')
		summary = key_split[0]
		source_file = '_'.join(key_split[1:])
		if analysis_name in summary:
			analysis_dic[key] = [summary, source_file, qc_dic[key]]	
	qc_file.close()
	return analysis_dic

##############################

# to plot average quality at each base
def average_base_quality(path, qc_file_name):
	try:
		analysis_name = 'average-base-quality'
		analysis_dic = qc_json(qc_file_name, analysis_name)
		
		plt.clf()				# clear and initialize a figure
		fig = plt.figure()
		fig.suptitle('Average Qscore at Each Base (Cycle)', size=10)
		
		keys =  analysis_dic.keys()
		keys.sort()
		
		subplot_counter = 0		# there will be 4 subplots (input_f, input_r, selection_f, selection_r)
		for key in keys:
			subplot_counter += 1
			ax = fig.add_subplot(2, 2, subplot_counter)
			ax.plot(analysis_dic[key][2])	# analysis_dic[key][2] is a list of average qscore at each base
			plt.title(analysis_dic[key][1], size=10)	# title is source_file_name
			if subplot_counter == 3 or subplot_counter == 4:	# show xlabel only in bottom subplots
				plt.xlabel('base(cycle)', size=10)
			if subplot_counter%2 != 0:							# show ylabel only in left subplots
				plt.ylabel('Average Qscore', size=10)
			plt.ylim(0,50)					# to set ylim to be same
			plt.axhline(y=20, color='r')	# draw a line to show qscore = 20
			plt.axhline(y=30, color='g')	# draw a line to show qscore = 30
			plt.tick_params(axis='both', labelsize=10)
		plt.subplots_adjust(hspace=0.3)		# to have enough space between top and bottom subplots
		plt.savefig(path + 'qc/average_base_quality_' + qc_file_name.split('/')[-1].rstrip('.qc') + '.pdf')
		#print 'average_base_quality: done'
		return 0
		
	except:
		print "Error: could not create plot of average base quality."
		return 1
		
##############################

# to plot mutation_identity at each base	
def mutation_identity(path, qc_file_name): 
	try:
		analysis_name = 'mutation-identity'
		analysis_dic = qc_json(qc_file_name, analysis_name)
	
		keys = analysis_dic.keys()
		for key in keys:	
			plt.clf()		# clear and initialize a figure
			fig = plt.figure()
			fig.suptitle(key, size=10)	# key is analysis_name+source_file
			fig.text(0.05, 0.65, 'Counts of Each Mutation', size=10, rotation=90)
			
			sub_dic_list = analysis_dic[key][2]		# analysis_dic[key][2] is a list of dictionaries
			total_base = len(sub_dic_list)			# calculate the length of total base
			subplot_row = total_base/10 + 1			# to determine the number of subplots
	
			counts_list=[]							# to find the maximum counts from all data
			for base in sub_dic_list:				
				counts_list += base.values()
			counts_list.sort()
			max_counts = counts_list[-1]
			ylim = max_counts*1.2					# ylim in subplots is around 120% of max_counts
			
			subplot_counter = 0				# subplot_counter is also base number
			for base in sub_dic_list:
				base_keys = base.keys()		# sort the base to be order, 'ACGT'
				base_keys.sort()
	
				y_value = []				# list of counts of each mutant in the list, 'y_value',
				empty_tick_label=[]			# to remove tick labels
				
				for DNA in base_keys:		# in the order, 'ACGT'
					y_value.append(base[DNA])
				x_value = range(len(y_value))
				
				subplot_counter += 1
				ax = fig.add_subplot(subplot_row, 10, subplot_counter)	# now generate subplots for each base
				ax.bar(x_value, y_value, align='center')				
				plt.xlim(-1,4)				# have space before and after the bars
				plt.ylim(0,ylim)
				plt.text(-0.8, max_counts*0.95, subplot_counter, size=6)	# to indicate base number
				
				ax.set_xticks(x_value)			# set xtick labels
				ax.set_xticklabels(['A', 'C', 'G', 'T'])
	
				if subplot_counter%20-1!= 0 :	# set ytick labels (show only in every second most-left subplots)
					plt.yticks(empty_tick_label)	
	
				plt.tick_params(axis='x', labelsize=6, direction = 'out', top = 'off')
				plt.tick_params(axis='y', labelsize=6, direction = 'out', right = 'off')
		
			plt.subplots_adjust(wspace=0.0001, hspace=0.0001)	# arrange subplots to look similar to a lattice plot 		
			plt.savefig(path + 'qc/' + key + '_' + qc_file_name.split('/')[-1].rstrip('.qc') + '.pdf')
		#print 'mutation_identity: done'
		return 0
	
	except:
		print "Error: could not create plot of mutation identities at each base."
		return 1
	
##############################

# to plot number of mutations at each base
def mutation_positions(path, qc_file_name):
	try:
		analysis_name = 'mutation-positions'
		analysis_dic = qc_json(qc_file_name, analysis_name)
		
		plt.clf()				# clear and initialize a figure
		fig = plt.figure()
		fig.suptitle('Number of Mutations Assayed at Each Position', size=10)
		subplot_counter = 0
		
		keys =  analysis_dic.keys()
		keys.sort()
		
		for key in keys:
			subplot_counter += 1
			ax = fig.add_subplot(2, 1, subplot_counter)
			ax.plot(analysis_dic[key][2])
			plt.title(analysis_dic[key][1], size=10)
			if subplot_counter == 2 :
				plt.xlabel('number of mutations', size=10)
			plt.ylabel('counts')
			plt.xlim(0, len(analysis_dic[key][2])*1.01)	
			plt.tick_params(axis='both', labelsize=10)
		plt.subplots_adjust(hspace=0.3)
		plt.savefig(path + 'qc/mutation_positions_' + qc_file_name.split('/')[-1].rstrip('.qc') + '.pdf')		
		#print 'mutation_positions: done'
		return 0
	
	except:
		print "Error: could not create plot of number of mutations per position."
		return 1
		
##############################

# to plot number of mutations per variant
def mutations_per_sequence(path, qc_file_name):
	try:
		analysis_name = 'mutations-per-sequence'
		analysis_dic = qc_json(qc_file_name, analysis_name)
		
		plt.clf()									# clear and initialize a figure
		fig = plt.figure()
		fig.suptitle('Number of Mutations per Variant', size=10)
		subplot_counter = 0
	
		xrange_list = []							# to determine the length of list to show in the subplots
		for key in analysis_dic.keys():
			while analysis_dic[key][2][-1] == 0:	# remove 0s from the end of the list until it gets to non-0
				analysis_dic[key][2].pop(-1)
			xrange_list.append(len(analysis_dic[key][2]))
		xrange_list.sort()
		xrange_max = xrange_list[-1]				# take a bigger number  
	
		keys =  analysis_dic.keys()
		keys.sort()
		
		for key in keys:
			subplot_counter += 1
			ax = fig.add_subplot(2, 1, subplot_counter)
			ax.plot(analysis_dic[key][2])
			plt.title(analysis_dic[key][1], size=10)
			if subplot_counter == 2:
				plt.xlabel('number of mutations', size=10)
			plt.ylabel('counts', size=10)
			plt.xlim(-1, xrange_max)	
			plt.tick_params(labelsize=10)
			
			x = range(len(analysis_dic[key][2]))	# to annotate y values (counts)
			y = analysis_dic[key][2]
			plt.plot(x,y)
			for i,j in zip(x,y):
				ax.annotate(str(j), xy=(i-0.2, j+0.5), size=8, color='b')
			
		plt.subplots_adjust(hspace=0.3)
		plt.savefig(path + 'qc/mutations_per_sequence_' + qc_file_name.split('/')[-1].rstrip('.qc') + '.pdf')		
		#print 'mutations_per_sequence: done'
		return 0
		
	except:
		print "Error: could not create histogram of mutations per variant."
		return 1

##############################

# to plot absolute number of variants assayed with N mutations
def total_variants_per_mutation_number(path, qc_file_name):
	try:	
		analysis_name = 'total-variants-per-mutation-number'
		analysis_dic = qc_json(qc_file_name, analysis_name)
		
		plt.clf()									# clear and initialize a figure
		fig = plt.figure()
		fig.suptitle('Absolute Number of Variants Assayed with N Mutations', size=10)
		subplot_counter = 0
	
		xrange_list = []							# to determine the length of list to show in the subplots
		for key in analysis_dic.keys():
			while analysis_dic[key][2][-1] == 0:	# remove 0s from the end of the list until it gets to non-0
				analysis_dic[key][2].pop(-1)
			xrange_list.append(len(analysis_dic[key][2]))
		xrange_list.sort()
		xrange_max = xrange_list[-1]				# take a bigger number
	
		keys =  analysis_dic.keys()
		keys.sort()
		
		for key in keys:
			subplot_counter += 1
			ax = fig.add_subplot(2, 1, subplot_counter)
			ax.plot(analysis_dic[key][2])
			plt.title(analysis_dic[key][1], size=10)
			if subplot_counter == 2:
				plt.xlabel('number of mutations', size=10)
			plt.ylabel('counts', size=10)
			plt.xlim(-1, xrange_max)	
			plt.tick_params(labelsize=10)
	
			x = range(len(analysis_dic[key][2]))	# to annotate y values (counts)
			y = analysis_dic[key][2]
			plt.plot(x,y)
			for i,j in zip(x,y):
				ax.annotate(str(j), xy=(i-0.2, j+0.5), size=8, color='b')
	
		plt.subplots_adjust(hspace=0.3)
		plt.savefig(path + 'qc/total_variants_per_mutation_number_' + qc_file_name.split('/')[-1].rstrip('.qc') + '.pdf')		
		#print 'total_variants_per_mutation_number: done'
		return 0
		
	except:
		print "Error: could not create total variants histogram."
		return 1
	
##############################

# to save separate files containing counts file information of readIDs which do appear in input library, but not in selection library	
def variants_not_in_sel(path, qc_file_name):
	analysis_name = 'variants-not-in-sel'
	analysis_dic = qc_json(qc_file_name, analysis_name)
	
	try:
		for key in analysis_dic.keys():					# two keys (DNA and PRO)
		  
			counts_file = open(path + '/data/output/' + analysis_dic[key][1], 'r') # find matching counts file
			header = counts_file.readline().strip()
			f_output = open(path + 'qc/not_in_sel_' + key.split('_')[-2] + '_' + qc_file_name.split('/')[-1].rstrip('.qc'),  'w')
			print >> f_output, header								# above is DNA or PRO
				
			counts_file_list = []
			counts_file_dic = {}
				
			counts_file_line = counts_file.readline()	# from the 2nd line
			while counts_file_line:						# generate a dictionary of counts file
				counts_file_line_list = counts_file_line.strip().split('\t')	
				counts_file_dic[counts_file_line_list[0]] = counts_file_line.strip()	# key is the seqID
				counts_file_line = counts_file.readline()				
				
			not_in_sel_list = analysis_dic[key][2]		# analysis_dic[key][2] is a list of seqID
			not_in_sel_list.sort()						# need a better sorting!
		
			for seqID in not_in_sel_list:				# generate a file of seqIDs w/ their 'counts file info'	
				print >> f_output, counts_file_dic[seqID]
										
			f_output.close()
			return 0	
			#print 'variants_not_in_sel: done'

	except:
			print "Error: could not compile list of variants found in input library, but not selected library."
			return 1

##############################	

# to save a file w/ basic sequencing summary information such as number of total reads, etc. 
def sequencing_summary(path, qc_file_name):
	try:
		f_output = open(path + 'qc/sequencing_summary', 'w')
		print >> f_output, 'sequencing summary', '\n', '------------------'
	
		analysis_list = ['total-reads', 'reads-passed', 'reads-filtered', 'total-ratios-all', 'total-ratios-m1', 'total-ratios-m2']
	
		for analysis_name in analysis_list:
			analysis_dic = qc_json(qc_file_name, analysis_name)
			
			keys = analysis_dic.keys()
			keys.sort()
			
			for key in keys:
				print >> f_output, key, ' : ', analysis_dic[key][2]
	
		f_output.close()
		return 0
		#print 'sequencing_summary: done'
	
	except:
		print "Error: could not compile summary information."
		return 1
		
##############################

# to analyze qc_removed file.
def qc_removed(path, tile_info_choice):		
	
	try:
		f_output = open(path + 'qc/' + 'qc_removed_summary', 'w')	# to report the summary of qc analysis
		
		tmp_directory = path+'data/tmp/'
		file_list = os.listdir(tmp_directory)	
		for file_name in file_list:
			file_name_split = file_name.split('_')
		
			if file_name_split[-1] == 'removed':				
				qc_removed_file = open(tmp_directory + file_name, 'r')
				header = qc_removed_file.readline()			# skip the header.
				
				total_qc_removed = 0
				readID_line = qc_removed_file.readline()
			#	readID_line_dic = {}
				error_code_sum = [0,0,0,0,0,0,0,0,0]	# initialize a list to calculate the error_code 
				error_code_readID_list = [[], [], [], [], [], [], [], [], []]
			
				while readID_line:	# to generate a dic of readID and error_code pairs
					total_qc_removed += 1
					readID_line_list = readID_line.strip().split('\t')
					error_code = readID_line_list[1] 	# readID is the key, error_code is the value.
					readID_line = qc_removed_file.readline() 		
		
					for i in range(0,9):
						error_code_sum[i] += int(error_code[i])	
						if int(error_code[i]) == 1:
							error_code_readID_list[i].append(readID_line_list[0])
												
				print >> f_output, file_name
				summary_list = ['read1 quality', 'read2 quality', 'min_quality', 'gap_count', 'unresolvable', 'maxmutrun', 'read1_Ncount', 'read2_Ncount', 'chastity']	
				print >> f_output, 'total_qc_removed :', '\t', total_qc_removed 
				i = -1
				for summary_name in summary_list:
					i += 1
					print >> f_output, summary_name, ': \t', error_code_sum[i]
				print >> f_output, '-----------------------------------------'
	
				if tile_info_choice == 'y':
					tile_info(path, file_name, error_code_readID_list)
	
		f_output.close()
		return 0
		
	except:
		print "Error: Could not perform qc_removed file parsing."
		return 1
		
#	print 'qc_removed: done'
	
##############################

# this is optional to generate tile info files of each qc in error_code_readID_list from qc_removed().
def tile_info(path, file_name, error_code_readID_list):
	try:
		tpath = path + 'qc/tile_info/'
		try:
			os.mkdir(tpath)
		except:
			os.chdir(tpath)
	
		summary_list = ['read1 quality', 'read2 quality', 'min_quality', 'gap_count', 'unresolvable', 'maxmutrun', 'read1_Ncount', 'read2_Ncount', 'chastity']

		i = 0	
		for summary_name in summary_list:
			f_output = open(path + 'qc/tile_info/' + file_name + '_' + summary_name + '_tile_info', 'w')
			print >> f_output, '\t'.join(['tile_name', 'X_coordinate', 'Y_coordinate'])

			for read_ID in error_code_readID_list[i]:
				read_ID_list = read_ID.split(':')
				print >> f_output, read_ID_list[2], '\t', read_ID_list[3], '\t', read_ID_list[4]
			i += 1

			f_output.close()
		return 0
		
	except:
		print "Error: could not compile tile information for qc_removed reads."
		return 1	
		
##############################

def main(path, qcfilename, tile_info_choice):
		
		retcodes = []
	
		qcfile = path + 'qc/' + qcfilename
			
		retcodes.append( average_base_quality(path, qcfile) )
		retcodes.append( mutation_identity(path, qcfile) )
		retcodes.append( mutation_positions(path, qcfile) )
		retcodes.append( mutations_per_sequence(path, qcfile) )
		retcodes.append( total_variants_per_mutation_number(path, qcfile) )
		retcodes.append( variants_not_in_sel(path, qcfile) )
		retcodes.append( sequencing_summary(path, qcfile) )
		retcodes.append( qc_removed(path, tile_info_choice) )
		
		if 1 in retcodes:
			return 1
		else:
			return 0

##############################

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option('--path', action = 'store', type = 'string', dest = 'path', help = 'path to project directory')
	parser.add_option('--qc_filename', action = 'store', type = 'string', dest = 'qc_filename', help = 'QC file')
	parser.add_option('--tile_info', action = 'store', type = 'string', dest = 'tile_info', help = 'Create files containing tile information for removed reads?')
	(option, args) = parser.parse_args()
	
	main(option.path, option.qc_filename, option.tile_info)
	
#qc_plot.py -p path_to_project_directory -c config_file_name
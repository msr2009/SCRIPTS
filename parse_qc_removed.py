"""
parse_qc_removed.py

parses *_index_filtered_B_qc_removed file (from Enrich read_aligner output) and returns number of each qc removal reason, 
lists of scores for creating histograms

creates two output files: ..._qc_removed-reasons and ..._qc_removed-qualityscores

python parse_qc_removed.py -f <FILE> -o <whichfilestooutput>
"""

def main():
	#parse list of output files
	outputs = set(option.output.split(','))
	#initialize dictionary to keep track of reasons for qc removal
	reasons = {}
	#open output files
	if 'quality' in outputs:
		q_out = open(option.infile+'-quality', 'w')
	
	if 'maxrun' in outputs:
		m_out = open(option.infile+'-maxmutrun', 'w')
	
	if 'unresolve' in outputs: 
		un_out = open(option.infile+'-unresolvable', 'w')
	
	if 'Ncount' in outputs:
		nc_out = open(option.infile+'-ncount', 'w')
		
	if 'gap' in outputs:
		g_out = open(option.infile+'-gap', 'w')
	
	#iterate through lines in file
	for line in open(option.infile, 'r'):
		try:
			qc_reason, score = line.strip().split('\t')[1].split('=')	
			#add to/increment reasons dictionary
			if qc_reason in reasons:
				reasons[qc_reason] += 1
			else:	
				reasons[qc_reason] = 1
			
			#output score to appropriate file
			if 'quality' in outputs:
				if qc_reason == 'read1_quality' or qc_reason == 'read2_quality':
					print >> q_out, score
	
			if 'maxrun' in outputs:		
				if qc_reason == 'maxmutrun_exceeded':
					print >> m_out, score
	
			if 'unresolve' in outputs:
				if qc_reason == 'unresolvable_max_exceeded':
					print >> un_out, score
			
			if 'Ncount' in outputs:
				if qc_reason == 'read1_Ncount_max_exceeded' or qc_reason == 'read2_Ncount_max_exceeded':
					print >> nc_out, score
		
			if 'gap' in outputs:
				if qc_reason == 'gap_count_exceeded':
					print >> g_out, score
		except ValueError:
			continue
	
	#open outfile for reasons counts
	f_out = open(option.infile+'-reasons', 'w')
				
	#print counts of reasons for removal of read
	for r in sorted(reasons.keys()):
		print >> f_out, '\t'.join([r, str(reasons[r])]) 
	
	try:
		f_out.close()
	except: 
		pass
	try:
		q_out.close()
	except:
		pass
	try:
		m_out.close()
	except:
		pass
	try:
		un_out.close()
	except:
		pass	
	try:
		nc_out.close()
	except:
		pass
	try:
		g_out.close()
	except:
		pass	

if __name__ == '__main__':
	from optparse import OptionParser
	import sys, time

	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'path to qc_removed file')
	parser.add_option('-o', '--output', action = 'store', type = 'string', dest = 'output', help = 'what files to output? (comma-delimit multiple options) possible options: quality,maxrun,unresolve,gap,Ncount,all')
	(option, args) = parser.parse_args()
	
	main()

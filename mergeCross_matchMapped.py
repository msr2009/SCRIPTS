"""
mergeCross_matchMapped.py

script takes input of two crossmatch .mapped files (output for parseCross_matchOutput.py), 
and merges reads with their paired-end counterparts (forward reads in one file, reverse in the other)

Matt Rich 03/2013
"""

def main(forward, reverse):

	#open each file
	f = open(forward, 'r')
	r = open(reverse, 'r')
	
	unmatched_reads = {} #a dictionary to store unmatched reads
	
	while True:
		lf = f.readline()
		lr = r.readline()
		
		if lf == '' and lr == '':
			break
	
		f_dat = lf.split('\t')
		r_dat = lr.split('\t')
		
		f_name = f_dat[0].split('#')[0]
		r_name = f_dat[0].split('#')[0]
		
		if f_name == r_name and f_dat[1] != r_dat[1]:
			print f_name, r_name, f_dat[1], r_dat[1]
		
		"""
		if f_name == r_name:
			#these are matching reads, so merge them
			print mergeDiscrepancies(f_dat, r_dat)
			
		else:
			#check if there's match to any of the unmatched reads
			if f_name in unmatched_reads:
				#merge them...
				print mergeDiscrepancies(f_dat, unmatched_reads[f_dat])	
				#...delete f_dat and unmatched_reads[f_dat]
				f_dat = []
				del unmatched_reads[f_dat]
			else:
				#otherwise, add f_dat to unmatched reads
				unmatched_reads[f_name] = f_dat[1:]
				
			if r_name in unmatched_reads:
				#merge them...
				print mergeDiscrepancies(r_dat, unmatched_reads[r_dat])
				#...delete r_dat and unmatched_reads[r_dat]
				r_dat = []
				del unmatched_reads[r_dat]
			else:
				#otherwise, add r_dat to unmatched reads
				unmatched_reads[r_name] = r_dat[1:]	
			"""
				
def mergeDiscrepancies(fread, rread):
	#@MISEQ:1:1101:14784:1937:1#0/2	643,532-A,T	S,S	243	2	ASR1
	#@MISEQ:1:1101:14784:1937:1#0/1	532-T	S	211	1	ASR1

	
	merged_read = []
	return

		
if __name__ == "__main__":
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f','--forward', action='store', type='str', dest='forward', help='forward read .mapped file')
	parser.add_option('-r','--reverse', action='store', type='str', dest='reverse', help='reverse read .mapped file')
	(option, args) = parser.parse_args()

	main(option.forward, option.reverse)



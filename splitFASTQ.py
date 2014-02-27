"""
splitFastq.py

Takes an input of a FASTQ file name, and splits reads based on a list of index sequences, as well as a mismatch threshold

Matt Rich, 2/2013
"""

def main(forward, reverse, index, mismatch, index_list):
	#create outfiles for all indices (and one for headers that don't match any indices close enough)
	outfiles = {}
	for x in index_list:
		outfiles[x] = [ open(forward.strip('.fq')+'_'+x+'.fq','w'), open(index.strip('.fq')+'_'+x+'.fq','w') ]
		if reverse != None:
			outfiles[x].append(open(reverse.strip('.fq')+'_'+x+'.fq','w'))
	
	outfiles['NNN'] = [ open(forward.strip('.fq')+'_NNN.fq','w'), open(index.strip('.fq')+'_NNN.fq','w') ]
	if reverse != None:
		outfiles['NNN'].append(open(reverse.strip('.fq')+'_NNN.fq','w'))

	#open fastq files		
	forw = open(forward, 'r')
	ind = open(index, 'r')
	if reverse != None:
		rev = open(reverse, 'r')
			
	#use read_fastq_multi to read teh fastq files
	for record in read_fastq_multi([forward, index, reverse]):
		f_read = record[0][1]
		i_read = record[1][1]
		r_read = record[2][1]
		
		found_index = False
		
		#check index sequence
                if i_read in index_list:
			print_fastq(record[0], outfiles[i_read][0])
			print_fastq(record[1], outfiles[i_read][1])
			print_fastq(record[2], outfiles[i_read][2])

                        found_index = True

                if found_index == False and mismatch != 0:
                        #check for close matches (within threshold)
                        diffs = [ hammingDistance(i_read, seq) for seq in index_list ]
                        for d in range(len(diffs)):
                                if diffs[d] <= mismatch:
					print_fastq(record[0], outfiles[index_list[d]][0])
		                        print_fastq(record[1], outfiles[index_list[d]][1])
                		        print_fastq(record[2], outfiles[index_list[d]][2])

                                        found_index = True
                                        break

                #no close match
                if found_index == False:
                        print_fastq(record[0], outfiles["NNN"][0])
                        print_fastq(record[1], outfiles["NNN"][1])
                        print_fastq(record[2], outfiles["NNN"][2])

    	#close all the files once we're done splitting            
        for f in outfiles:      
                outfiles[f][0].close()
                outfiles[f][1].close()
                try:
                        outfiles[f][2].close()
                except:
                        pass

def hammingDistance(seq1, seq2):
	return len(seq1) - sum([ seq1[x] == seq2[x] for x in range(len(seq1))])
	
if __name__ == '__main__':
	from optparse import OptionParser
	from fastq_tools import read_fastq_multi, print_fastq
	parser = OptionParser()
	parser.add_option('--forward', action = 'store', type = 'string', dest = 'forward', help = 'path to forward reads file')
	parser.add_option('--reverse', action = 'store', type = 'string', dest = 'reverse', help = 'path to reverse reads file', default=None)
	parser.add_option('--index', action = 'store', type = 'string', dest = 'index', help = 'path to index reads file')
	parser.add_option('--mismatch', action = 'store', type = 'int', dest = 'mismatch', default = 1, help = 'number of permitted mismatches (Default=1)')
	parser.add_option('--index_list', action = 'store', type = 'string', dest = 'index_list', help = 'comma-delimited list of index seqeunces')
	(option, args) = parser.parse_args()
	
	main(option.forward, option.reverse, option.index, option.mismatch, option.index_list.split(','))	
	

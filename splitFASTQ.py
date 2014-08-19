"""
splitFastq.py

Takes an input of a FASTQ file name, and splits reads based on a list of index sequences, as well as a mismatch threshold

Matt Rich, 2/2013
"""

def main(forward, reverse, index, mismatch, index_list, name_list):
	#create outfiles for all indices (and one for headers that don't match any indices close enough)
	outfiles = {}
	for x in range(len(index_list)):
		outfiles[index_list[x]] = [ open(forward.strip('.fq')+'_'+name_list[x]+'.fq','w'), open(index.strip('.fq')+'_'+name_list[x]+'.fq','w') ]
		if reverse != None:
			outfiles[index_list[x]].append(open(reverse.strip('.fq')+'_'+name_list[x]+'.fq','w'))
	
	outfiles['NNN'] = [ open(forward.strip('.fq')+'_NNN.fq','w'), open(index.strip('.fq')+'_NNN.fq','w') ]
	if reverse != None:
		outfiles['NNN'].append(open(reverse.strip('.fq')+'_NNN.fq','w'))
	
	#open fastq files		
	forw = open(forward, 'r')
	ind = open(index, 'r')
	if reverse != None:
		rev = open(reverse, 'r')
	
	counter = 0		
	#use read_fastq_multi to read teh fastq files
	for record in read_fastq_multi([forward, index, reverse]):
		f_read = record[0][1]
		i_read = record[1][1]
		if reverse != None:
			r_read = record[2][1]
		
		found_index = False
		
		#check index sequence
                if i_read in index_list:
			print_fastq(record[0], outfiles[i_read][0])
			print_fastq(record[1], outfiles[i_read][1])
			if reverse != None:
				print_fastq(record[2], outfiles[i_read][2])

                        found_index = True

                if found_index == False and mismatch != 0:
			#check for close matches (within threshold)
                        diffs = [ hammingDistance(i_read, seq) for seq in index_list ]
                        for d in range(len(diffs)):
                                if diffs[d] <= mismatch:
					print_fastq(record[0], outfiles[index_list[d]][0])
		                        print_fastq(record[1], outfiles[index_list[d]][1])
                		        if reverse != None:
						print_fastq(record[2], outfiles[index_list[d]][2])

                                        found_index = True
                                        break

                #no close match
                if found_index == False:
                        print_fastq(record[0], outfiles["NNN"][0])
                        print_fastq(record[1], outfiles["NNN"][1])
			if reverse != None:
				 print_fastq(record[2], outfiles["NNN"][2])
		
		#increment counter
		counter += 1.0
		if counter % 100000 == 0:
			divcount = counter/1000000.0
			print 'split %0.1fM reads' % divcount

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

def readIndexFile(f): 
	d = []
	n = []

	for line in open(f, 'r'):
		l = line.strip().split('\t')
		d.append(l[1])
		n.append(l[0])
	return d, n

if __name__ == '__main__':
	from optparse import OptionParser
	from fastq_tools import read_fastq_multi, print_fastq

	parser = OptionParser()
	parser.add_option('--forward', action = 'store', type = 'string', dest = 'forward', help = 'path to forward reads file')
	parser.add_option('--reverse', action = 'store', type = 'string', dest = 'reverse', help = 'path to reverse reads file', default=None)
	parser.add_option('--index', action = 'store', type = 'string', dest = 'index', help = 'path to index reads file')
	parser.add_option('--mismatch', action = 'store', type = 'int', dest = 'mismatch', default = 1, help = 'number of permitted mismatches (Default=1)')
	parser.add_option('--index_list', action = 'store', type = 'string', dest = 'index_list', help = 'comma-delimited list of index sequences, or tab-delimited file of name,index pairs')
	parser.add_option('--name_list', action = 'store', type = 'string', dest = 'name_list', help = 'comma-delimited list of names for fastq files', default=None)
	(option, args) = parser.parse_args()
	
	id = []
	names = None	

	if option.reverse == "None":
		option.reverse = None
	
	try:
		id, names = readIndexFile(option.index_list)
	except IOError:
		print "IndexList isn't a file, parsing as comma-delimited list"
		try:
			names = option.name_list.split(',')
		except AttributeError:
			pass
		id = option.index_list.split(',')
	
#	print ",".join(map( lambda x: x[0] + ":" + x[1], zip(names, id) ))
	
	if names:	
		main(option.forward, option.reverse, option.index, option.mismatch, id, names)
	else:
		main(option.forward, option.reverse, option.index, option.mismatch, id, id)	

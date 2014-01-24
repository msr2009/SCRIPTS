#/usr/bin/python
"""
splitFASTQ_bySAM.py

Script splits FASTQ file (forward and reverse read files) depending on where reads map, as shown in SAM file

Matt Rich 7/2013

"""

def main(forward, reverse, sam, map_list):
	#create outfiles for all indices (and one for headers that don't match any indices close enough)
	outfiles = {}
	for x in map_list:
		outfiles[x] = [ open(forward.strip('.fq')+'_'+x+'.fq','w'), open(reverse.strip('.fq')+'_'+x+'.fq','w')]
	
	#open fastq files		
	forw = open(forward, 'r')
	rev = open(reverse, 'r')
	sam = open(sam, 'r')

	#start looping through files
	#sam file = 2 lines/iteration, FASTQs = 4 lines/iteration	
		
	while True:
		
		s_read = [ sam.readline().strip() for i in range(2) ]
		if s_read[0] == '':
			break
		
		f_read = [ forw.readline().strip() for i in range(4) ]
		r_read = [ rev.readline().strip() for i in range(4) ]  

		#check that SAM file has matching paired reads (they will ahve the same ID)
		s_split = [s_read[i].split('\t') for i in range(2)]
		if s_split[0][0] != s_split[1][0]:
			print "Paired IDs do not match! " + s_split[0][0] + ", " + s_split[1][0]
			break
		
		#get mapping locations 
		s_map = [ s_split[i][2] for i in range(2) ]
		#check that mapping locations are the same
		if s_map[0] != s_map[1]:
			print "Mapped locations do not match! " + s_split[0][0] + ", " + s_split[1][0]
			break
		
		#otherwise, print out FASTQ lines to appropriate files
		if s_map[0] in map_list:
			print >> outfiles[s_map[0]][0], '\n'.join(f_read)
			print >> outfiles[s_map[0]][1], '\n'.join(r_read)
		else:
			print "Found mapping location that wasn't in --map_list: " + s_map[0]
		
	#close files	
	for f in outfiles:	
		outfiles[f][0].close()
		outfiles[f][1].close()
	
if __name__ == '__main__':
	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('--forward', action = 'store', type = 'string', dest = 'forward', help = 'path to forward reads file')
	parser.add_option('--reverse', action = 'store', type = 'string', dest = 'reverse', help = 'path to reverse reads file')
	parser.add_option('--sam', action = 'store', type = 'string', dest = 'sam', help = 'path to SAM file')	
	parser.add_option('--map_list', action = 'store', type = 'string', dest = 'map_list', help = 'comma-delimited list of mapping locations')
	(option, args) = parser.parse_args()
	
	main(option.forward, option.reverse, option.sam, option.map_list.split(','))	

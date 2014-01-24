"""
illumina_qc.py

scavenging from Enrich to grab a couple qc functions to be able to qc-test your illumina FASTQs.

illumina_qc.py --fastq FASTQ --output OUTPUTPATH


also output reasonable summary stats (#reads, fraction passing read/base qscore, etc..)
variance at each position (keep track of histogram of quality scores for each base)

"""

import optparse, json, array
		
def main(filename):	

	prefix = filename.split('.')[0]

	lenread = 100000 #lenread specifies the number of bytes read in from each input file 
	lenbit = 1 #lenbit is a flag used to indicate when the last chunk of input from the file has been read
	tile_info = {}
	score_adjustment = 33
	read_counter = 1.0
	read1_means = []
	remainder_read1 = ''
	firstrun = 1
	read_range = []
	base_quality = []
	
	try:
		f_read1 = open(filename, 'U')
	except:
		print "Error: could not open FASTQ file."
		return
	
	while lenbit:
		#read in chunks of each input file and prepend the remainder, which is the incomplete set of lines left over from the previous loop
		read1_bytes = remainder_read1 + f_read1.read(lenread)

		#check to see if the end of file has arrived (i.e. if the length of the read from the input file is less than the specified length
		if len(read1_bytes) < lenread:
			lenbit = 0
	
		#split the read into lines 
		read1_lines = read1_bytes.split('\n')
		
		#iterate over the split reads in four line chunks (since each illumina read consists of four lines)
		for i in range(0, ((len(read1_lines)/4)-lenbit)):
			
			#take the first four lines from each read and then remove them from the list
			read1_entry = read1_lines[0:4]
			read1_lines[0:4] = []
			
			#check read lengths if first time, create read1/read2_means lists
			if firstrun == 1:
				len_read1 = len(read1_entry[1])
				read_range = range(len_read1)
				for baselength in read_range:
#					read1_means.append(0)
					base_quality.append({})
				firstrun = 0
		
			#convert ascii quality scores for read to qscore	
			read1_quality = [score - score_adjustment for score in array.array('b', read1_entry[3]).tolist()]
			#calculate average quality of read
			avg_quality = sum(read1_quality)/float(len_read1)
			
			#increment tile counters based on average quality of read
			read_ids = read1_entry[0].split(':')
			tile = read_ids[2]
			chaste = read_ids[-1][0]
			if tile not in tile_info:
				tile_info[tile] = [0,0,0,0]
			else:
				if avg_quality >= 30:
					tile_info[tile][2] += 1
				elif avg_quality < 20:
					tile_info[tile][0] += 1
				else:
					tile_info[tile][1] += 1
			tile_info[tile][3] += int(chaste)
			
			#increment average score for each base
#			foo = [ ((read_counter-1)/read_counter)*x for x in read1_means ]
#			bar = [ x/read_counter for x in read1_quality ]
#			read1_means = [ foo[i] + bar[i] for i in read_range ]
			
			for s in range(len(read1_quality)):
				if read1_quality[s] in base_quality:
					base_quality[read1_quality[s]] += 1
				else:
					base_quality[read1_quality[s]] = 1 

			read_counter += 1
			
		remainder_read1 = '\n'.join(read1_lines)
		
	#print QC information to JSON-formatted file
	qc_dict = {'average_quality': read1_means, 'tile_quality': tile_info}
	print >> open(prefix+'_qc-data', 'w'), json.dumps(qc_dict)
	plot_qc(json.dumps(qc_dict), prefix)

def plot_qc(json_str, prefix):
	from numpy import arange
	
	try:
		import matplotlib 
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
	except:
		print 'Error: could not load matplotlib.  Plots will not be produced.'
		return
	
	qc_dic = json.loads(json_str)
	
	#make average quality plot
	av_len = len(qc_dic['average_quality'])
	
	plt.clf()
	plt.plot(range(1,av_len+1), qc_dic['average_quality'], 'b')
	plt.axis((0,102,2,40))
	plt.ylabel('Average base quality')
	plt.xlabel('Position')
	plt.title('Average quality along read: '+prefix.split('/')[-1])
	plt.savefig(prefix+'_average-quality.pdf')
	
	#make tile info plot
	tile_dic = qc_dic['tile_quality']
	tiles = tile_dic.keys()
	tile_0 = []
	tile_20 = []
	tile_30 = []
	tile_ch = []
	sum_tile = []
	pos = arange(len(tiles))
	width = .9	
		
	for i in sorted(tiles):
		tile_0.append(tile_dic[i][0])
		tile_20.append(tile_dic[i][1])
		tile_30.append(tile_dic[i][2])
		tile_ch.append(tile_dic[i][3])	
		sum_tile.append(tile_dic[i][0]+tile_dic[i][1])
		
	plt.clf()
	pl0 = plt.bar(pos, tile_0, width, color='r', bottom=0 )
	pl20 = plt.bar(pos, tile_20, width, color='b', bottom=tile_0)
	pl30 = plt.bar(pos, tile_30, width, color='g', bottom=sum_tile)
	x1,x2,y1,y2 = plt.axis()
	plt.axis((pos[0],pos[-1]+width,y1,y2))
	plt.xlabel('Tile')
	plt.xticks(pos+width/2.0, sorted(tiles), rotation=90)
	plt.ylabel('Read frequency')
	plt.title('Read quality by tile: '+prefix.split('/')[-1])
	plt.legend( (pl0[0], pl20[0], pl30[0]), ('Q>20','20<Q<30','Q>30'), loc=9, ncol=3, prop={'size':8} )
	
	plt.savefig(prefix+'_tile-quality.pdf', bbox_inches='tight')

	#make chaste plot (by tile)
	#make chaste data %'s of total reads in lane
	tile_chaste = [ tile_ch[i]/(sum([float(tile_0[i]), float(tile_20[i]), float(tile_30[i])])) for i in range(len(tiles)) ]
	plt.clf()
	plch = plt.bar(pos, tile_chaste, width, color='grey')
	x1,x2,y1,y2 = plt.axis()
	plt.axis((pos[0],pos[-1]+width,y1,y2))
	plt.xlabel('Tile')
	plt.xticks(pos+width/2.0, sorted(tiles), rotation=90)
	plt.ylabel('% of reads passing chastity filter')
	plt.title('Read chastity by tile: ' + prefix.split('/')[-1])
	plt.savefig(prefix+'_tile-chastity.pdf',  bbox_inches='tight')
	
if __name__ == '__main__':
						
	parser = optparse.OptionParser()
	parser.add_option('--fastq', action = 'store', type = 'string', dest = 'fastq', help = 'Path to fastq file')
	parser.add_option('--plot_only', action = 'store', type = 'string', dest = 'plot', default=False, help = 'Do not run FASTQ analysis, only plot from existing data. Argument takes path to qc_data file.')
	(option, args) = parser.parse_args()
	
	if option.plot != False:
		print "Only plotting.  Using " + option.plot.split('/')[-1] + "."
		plot_qc(open(option.plot, 'r').readline(), option.plot.split('.')[0].split('_qc-data')[0])
	else:
		main(option.fastq)

	
"""
q_mutation.py

identifies quality scores of mutations in reads

python q_mutation.py --fastq FASTQFILE --merged MERGEDFILE --output OUTPUTFILE

"""

def main():

	print time.asctime(time.localtime())
	#read merged reads file first, make dictionary of readID: mutation positions
	d_merged = {}
	merged_file = open(option.merged,'r')
	for line in merged_file:
		if line.startswith('@') == True:
			l = line.strip().split('\t')
			if l[4] != 'NA':
				d_merged[l[0]] = l[4].split(',')
	merged_file.close()
	
	#read fastq file, store scores in list of lists to create histograms
	#create scores list
	scores = []
	for i in range(101):
		tmp = []
		for i in range(42):
			tmp.append(0)
		scores.append(tmp)

	#open file, read ID line, see if ID in d_merged, add position score data	
	count = 0
	found_ID = False

	fqfile = open(option.fastq, 'r')
	
	for line in fqfile:
		#find ID that's in d_merged
		if found_ID == False:
			if count % 4 == 0:
				readID = line.strip().split('\t')[0]
				if readID in d_merged:
					found_ID = d_merged[readID]
		#read next line
		elif found_ID != False:
			if count % 4 == 3:
				qs = list(line.strip())
				for pos in found_ID:
					pos = int(pos)
					
					try:
						sc = ord(qs[pos])-33
					except IndexError:
						continue
										
					if sc < 41:
						try:
							scores[pos][sc] += 1	
						except IndexError:
							continue
					
					else:
						try:
							scores[pos][42] += 1
						except IndexError:
							continue
				found_ID = False
		count += 1
	
	fqfile.close()

	#open out_file
	out_file = open(option.f_output_folder, 'w')
	for i in range(101):
		s = map(lambda x: str(x), scores[i])
		s.insert(0, i)
		print >> out_file, '\t'.join(map(lambda x: str(x), s))
		
	print time.asctime(time.localtime())
	
if __name__ == '__main__':
	from optparse import OptionParser
	import sys, time

	parser = OptionParser()
	parser.add_option('--fastq', action = 'store', type = 'string', dest = 'fastq', help = 'path to fastq reads file')
	parser.add_option('--merged', action = 'store', type = 'string', dest = 'merged', help = 'path to merged reads file')
	parser.add_option('--output', action = 'store', type = 'string', dest = 'f_output_folder', help = 'path to output file')
	(option, args) = parser.parse_args()
	
	main()
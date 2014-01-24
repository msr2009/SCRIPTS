"""
fastq_index_histogram.py

plots histogram of all sequences in a FASTQ file.

Matt Rich 11/13

"""

def main(fastq):
	#store indices
	indices = {}
	
	#read each fastq record
	for record in read_fastq(fastq, None, 250000):
		if record[1] in indices:
			indices[record[1]] += 1
		else:	
			indices[record[1]] = 1

	#plotting
	plt.figure()
	plt.axis([1, 100, 0.5, 1.0])
	plt.xlabel("Reads")
	plt.ylabel("Frequency")
	plt.title(fastq.split('/')[-1])
	plt.hist(indices.values(), bins= max(indices.values()), cumulative=True, normed=True)
	plt.savefig(fastq+".counts_histogram.pdf")
	
	
if __name__ == "__main__":
	from fastq_tools import read_fastq
	from optparse import OptionParser
	import matplotlib.pyplot as plt
	import operator

	parser = OptionParser()
	parser.add_option("-f", "--fastq", action = "store", type = "string", dest = "fastq", help = "fastq file")
	
	(option, args) = parser.parse_args()

	main(option.fastq)

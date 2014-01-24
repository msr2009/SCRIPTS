#/usr/bin/python

"""
simulate_barcodes.py

script creates library of barcodes of specified length, samples specified number of barcodes, 
and checks all-by-all comparisons for longest common substrings.

python simulate_barcodes.py -l LENGTH -n NUMBER_SAMPLED -i ITERATIONS -t TOTAL_LIB_SIZE(def=1e8)

"""

def main(libsize, bclength, samplesize, itera):
	
	#create bc library
	#allows for duplicated BCs, since that's how they'd be made
	barcodes = []
	for i in range(libsize):
		barcodes.append(createSequence(bclength))
	
	#list to store homology lengths
	microhomologies = []
	
	count = 0
	#do sampling, testing
	for countIterations in range(itera):
		#sample barcodes
		bcs = sample(barcodes, samplesize)
		#try all-by-all comparisons
		for i in range(len(bcs)):
			if count % 1 == 0:
				print count
			for j in range(len(bcs)):
				if i != j:
					LCS = LongestCommonSubstring(bcs[i], bcs[j])
					if LCS != 1:
						microhomologies.append(len(LCS))
			count += 1
	
	#make a histogram of the results
	homohist = makeHistogram(microhomologies)
	for i in sorted(homohist.keys()):
		print '\t'.join([str(i), str(homohist[i])])	
	
def makeHistogram(dat):
	h = {}
	for i in range(min(dat), max(dat)+1):
		h[i] = 0
	for i in dat:
		h[i] += 1
	return h

def createSequence(l):
	bases = ['A','C','T','G']
	seq = ''
	for i in range(l):
		seq += sample(bases,1)[0]
	return seq

def LongestCommonSubstring(S1, S2, offset=0):
	M = [[0]*(1+len(S2)) for i in xrange(1+len(S1))]
	longest, x_longest = 0, 0
	for x in xrange(1,1+len(S1)):
		for y in xrange(1,1+len(S2)):
			if S1[x-1] == S2[y-1]:
				M[x][y] = M[x-1][y-1] + 1
				if M[x][y]>longest:
					longest = M[x][y]
					x_longest  = x
			else:
				M[x][y] = 0
	return S1[x_longest-longest: x_longest]

if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	from random import sample
	
	parser = OptionParser()
	parser.add_option('-l', '--length', action = 'store', type = 'int', dest = 'barcode_length', help = 'barcode length')
	parser.add_option('-n', '--number', action = 'store', type = 'int', dest = 'number_sampled', help = 'number of barcodes to sample')
	parser.add_option('-i', '--iterations', action = 'store', type = 'int', dest = 'iterations', help = 'number of iterations to run')
	parser.add_option('-t', '--total', action = 'store', type = 'int', dest = 'total_size', help = 'size of initial barcode library (default = 1e8)')
	(option, args) = parser.parse_args()

	length = 100000000
	iterations = 1
	
	if option.barcode_length != None:
		length = option.barcode_length
	if option.iterations != None:
		iterations = option.iterations

	main(option.total_size, length, option.number_sampled, iterations)
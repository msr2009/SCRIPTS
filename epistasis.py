#/usr/bin/python

"""
epistasis.py

calculates epistatis scores (AB - A*B) for all mutation pairs. Using input of ratios .m1 and .m2 files, finds all possible epistatis scores

python epistatis.py -f PATHTOALLMUTSFILE

"""

def main():
	
	#lists to store fitness scores 
	singles = {} #can use a dictionary for speed, since there's only on entry for each variant
	doubles = []
	epistasis = []
	
	#open, store single mutants file
	for line in open(option.infile+'.m1', 'r').readlines()[1:]:
		l = line.strip().split('\t')
		sID = l[0].split('-')
		singles[sID[0]+sID[1]] = float(l[9])
		
	#open, store double mutants file
	for line in open(option.infile+'.m2', 'r').readlines()[1:]:
		l = line.strip().split('\t')
		sID = map(lambda x: x.split(','), l[0].split('-'))
		doubles.append([ ','.join([ sID[0][i]+sID[1][i] for i in range(len(sID[0])) ]), float(l[9]) ])
	
	for k in doubles:
		[a,b] = k[0].split(',')
		if a in singles and b in singles:
			epistasis.append( [ k[0], k[1], singles[a] * singles[b],  k[1] - (singles[a] * singles[b]) ] ) #[ k1_seqID, k1_fitness, pred_k1_fitness, epistasis ]				
			
	#open output file
	f_out = open(option.infile+'-epistasis', 'w')
	
	#print output	
	print >> f_out, '\t'.join(['seqID', 'double_fitness', 'predicted_fitness', 'epistasis'])
	for i in epistasis:
		print >> f_out, '\t'.join(map(str, i))

	f_out.close()
		
if __name__ == '__main__':

	from optparse import OptionParser
	
	parser = OptionParser()
	parser.add_option('-f', '--file', action = 'store', type = 'string', dest = 'infile', help = 'path to file containing all fitness values')
	(option, args) = parser.parse_args()

	main()
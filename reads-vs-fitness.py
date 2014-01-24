
def main():
	dat = {}
	
	bins = [0,1,5,10,50,100,500,1000,10000,100000,1000000,10000000]
	
	for line in open(option.fit,'r'):
		if line.startswith('seqID') == True:
			continue
		l = line.strip().split('\t')
		dat[l[0]] = [l[7]]
	
	for line in open(option.counts, 'r'):
		if line.startswith('seqID') == True:
			continue
		l = line.strip().split('\t')
		
		if l[0] in dat:
			dat[l[0]].append(l[8])
			for i in range(1, len(bins)):	
				if int(l[8]) <= bins[i] and int(l[8]) > bins[i-1]:	
					dat[l[0]].append(i)
			
	for i in dat:
		print '\t'.join([i, str(dat[i][0]), str(dat[i][1]), str(dat[i][2])])
	
if __name__ == '__main__':

	from optparse import OptionParser
		
	parser = OptionParser()
	parser.add_option('-f', '--fitness', action = 'store', type = 'string', dest = 'fit', help = 'ratios file')
	parser.add_option('-c', '--counts', action = 'store', type = 'string', dest = 'counts', help = 'counts file')
	(option, args) = parser.parse_args()

	main()
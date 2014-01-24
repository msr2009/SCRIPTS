def checkstarts(genes, o=0, oo=0):
	no_start_end = 0
	correct = 0
	o = 0
	oo = 0
	bad = {}
	stops = set(['TAA','TAG','TGA'])
	for i in genes:
		try:
			if genes[i][0].split('-')[0][-1] == 'W':
				start = int(genes[i][2])
				end = int(genes[i][3])	
				if c[genes[i][1]-1][start-1+o:start+2+o] == 'ATG' and c[genes[i][1]-1][end-3+oo:end+oo] in stops:
					correct += 1
					print genes[i][0], c[genes[i][1]-1][start-1+o], c[genes[i][1]-1][end-3+oo], len(genes[i][4])
				else:
					no_start_end += 1
					bad[i] = genes[i]
			else:
				start = int(genes[i][3])
				end = int(genes[i][2])	
				if complement(c[genes[i][1]-1][start-3+o:start+0+o][::-1]) == 'ATG' and complement(c[genes[i][1]-1][end-1+o:end+2+oo][::-1]) in stops:
					correct += 1
					print genes[i][0], complement(c[genes[i][1]-1][start-1][::-1]), complement(c[genes[i][1]-1][end+1+o][::-1]), len(genes[i][4])
				else:
					no_start_end += 1
					bad[i] = genes[i]
		except TypeError:
			pass
	print correct, no_start_end

"""
W: A = start-1, T = end-3
C: A = start-1, T = end+1
"""
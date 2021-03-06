#/usr/bin/python

"""

"""




def matrix(res1, res2):
	mat = { 'XXX': ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*'],
			'A': [5,-2,-2,-2,-1,-1,-1,0,-2,-2,-2,-1,-1,-3,-1,1,0,-3,-2,0,-6],
			'R': [-2,6,-1,-2,-4,1,-1,-3,0,-3,-3,2,-2,-4,-2,-1,-1,-4,-3,-3,-6],
			'N': [-2,-1,6,1,-3,0,-1,-1,0,-4,-4,0,-3,-4,-3,0,0,-4,-3,-4,-6],
			'D': [-2,-2,1,6,-4,-1,1,-2,-2,-4,-5,-1,-4,-4,-2,-1,-1,-6,-4,-4,-6],
			'C': [-1,-4,-3,-4,9,-4,-5,-4,-4,-2,-2,-4,-2,-3,-4,-2,-1,-3,-3,-1,-6],
			'Q': [-1,1,0,-1,-4,6,2,-2,1,-3,-3,1,0,-4,-2,0,-1,-3,-2,-3,-6],
			'E': [-1,-1,-1,1,-5,2,6,-3,0,-4,-4,1,-2,-4,-2,0,-1,-4,-3,-3,-6],
			'G': [0,-3,-1,-2,-4,-2,-3,6,-3,-5,-4,-2,-4,-4,-3,-1,-2,-4,-4,-4,-6],
			'H': [-2,0,0,-2,-4,1,0,-3,8,-4,-3,-1,-2,-2,-3,-1,-2,-3,2,-4,-6],
			'I': [-2,-3,-4,-4,-2,-3,-4,-5,-4,5,1,-3,1,-1,-4,-3,-1,-3,-2,3,-6],
			'L': [-2,-3,-4,-5,-2,-3,-4,-4,-3,1,4,-3,2,0,-3,-3,-2,-2,-2,1,-6],
			'K': [-1,2,0,-1,-4,1,1,-2,-1,-3,-3,5,-2,-4,-1,-1,-1,-4,-3,-3,-6],
			'M': [-1,-2,-3,-4,-2,0,-2,-4,-2,1,2,-2,6,0,-3,-2,-1,-2,-2,1,-6],
			'F': [-3,-4,-4,-4,-3,-4,-4,-4,-2,-1,0,-4,0,6,-4,-3,-2,0,3,-1,-6],
			'P': [-1,-2,-3,-2,-4,-2,-2,-3,-3,-4,-3,-1,-3,-4,8,-1,-2,-5,-4,-3,-6],
			'S': [1,-1,0,-1,-2,0,0,-1,-1,-3,-3,-1,-2,-3,-1,5,1,-4,-2,-2,-6],
			'T': [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-2,-1,-1,-2,-2,1,5,-4,-2,0,-6],
			'W': [-3,-4,-4,-6,-3,-3,-4,-4,-3,-3,-2,-4,-2,0,-5,-4,-4,11,2,-3,-6],
			'Y': [-2,-3,-3,-4,-3,-2,-3,-4,2,-2,-2,-3,-2,3,-4,-2,-2,2,7,-2,-6],
			'V': [0,-3,-4,-4,-1,-3,-3,-4,-4,3,1,-3,1,-1,-3,-2,0,-3,-2,4,-6],
			'*': [-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,1],
	 		}
	return mat[res1.upper()][mat['XXX'].index(res2.upper())]
	
def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
             'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
             'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
             'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
             'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
             'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
             'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
             'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
             'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
             'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
             'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
             'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
             'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
             'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
             'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
             'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' }
	return lookup[codon.lower()]

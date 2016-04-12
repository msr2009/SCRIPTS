"""
ligand_selection_simulation.py

simulates a selection experiment 

(the first case is for a two-sided growth selection, with positive effects (for
stable proteins), and then with negative effects for stable proteins)

Matt Rich, 3/2016
"""

def main():

	#set up population
	n_strains = 4
	population = [ 1.0/n_strains for n in range(n_strains) ]
	cost1 = .4
	cost2 = .2


	#we have these strains:
	#0) stable w/ ligand, stable w/o
	#1) stable w/ ligand, unstable w/o
	#2) unstable w/ ligand, stable w/o
	#3) unstable w/ ligand, unstable w/o
	strain_phenotypes = [	
		[1.1, 1, 1-cost1, 1-cost1],
		[1-cost2, 1, 1-cost2, 1]
	]
	
	sel_condition = 0
	generation = 0
	#run selection
	for i in range(10): 	#do ten cycles
		for j in range(5): 	#five generations per condition
			#calculate fitness change for generation
			fitness = [a*b for a,b in 
				zip(population, strain_phenotypes[sel_condition])]
			population = ( [ f/sum(fitness) for f in fitness] )
			generation += 1
			print "\t".join([str(generation)] + [str(x) for x in population])
	 
		sel_condition = 1 - sel_condition

if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('--OPT', action = 'store', type = str, dest = 'DEST', 
		help = "HELP")
	args = parser.parse_args()
	
	main()	


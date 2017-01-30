from itertools import chain, combinations, product

def get_powerset(vars):
	'''Gets all combinations of variables in a peptide.'''
	# Gotta get a little fancier
	
	# I want to keep consistent positions, so taking things out of the dictionary for now
	save_varkeys = vars.keys()
	save_varvals = []
	for key in save_varkeys:
		save_varvals.append(vars[key])
	
	# Add a blank to each
	for v in save_varvals:
		v.append([])
	
	# Get the full powerset, including blanks
	full_powerset = list(product(*save_varvals))
	
	variant_sets = []
	for variant in full_powerset:
		new_set = []
		for i in range(len(variant)):
			if variant[i] != []:
				new_set.append([save_varkeys[i], variant[i][1]])
		if new_set != []:
			variant_sets.append(new_set)
	
	return variant_sets
	
vars = {1: [['a','A']], 2: [['b','B'],['c','C']], 3: [['d','D'],['e','E']]}

print get_powerset(vars)

#[(['a', 'A'], ['b', 'B'], ['d', 'D']), (['a', 'A'], ['b', 'B'], ['e', 'E']), (['a', 'A'], ['b', 'B'], []), (['a', 'A'], ['c', 'C'], ['d', 'D']), (['a', 'A'], ['c', 'C'], ['e', 'E']), (['a', 'A'], ['c', 'C'], []), (['a', 'A'], [], ['d', 'D']), (['a', 'A'], [], ['e', 'E']), (['a', 'A'], [], []), ([], ['b', 'B'], ['d', 'D']), ([], ['b', 'B'], ['e', 'E']), ([], ['b', 'B'], []), ([], ['c', 'C'], ['d', 'D']), ([], ['c', 'C'], ['e', 'E']), ([], ['c', 'C'], []), ([], [], ['d', 'D']), ([], [], ['e', 'E']), ([], [], [])]
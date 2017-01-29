'''
I will probably never need this again, but I made it to check whether my fasta files (both the full protein one and the tryptic protein one) were keeping different amino acid substitutions at the same position (the tryptic one wasn't at the time, the full protein one was). Also found through this that the full protein one wasn't keeping two amino acid substitutions at the same position in consecutive genes, mostly because I didn't really think that would ever happen (it did). Fixed that.
'''

import re
saved = []
savedset = set([])
f = open('proteome.aa.var.bed.dna.fasta','r')
line = f.readline()
while line:
	if line[0] == '>':
		name = line.split('-')[0]
		var = line.split(':')[-2]
		varpos = int(re.findall(r'\d+', var)[0])
		tosave = name+'-'+str(varpos)
		saved.append(tosave)
		savedset.add(tosave)
	line = f.readline()
print len(saved), len(savedset)
f.close()

savedset_orig = savedset

saved = []
savedset = set([])
f = open('tryptic_proteome.fasta','r')
line = f.readline()
while line:
        if line[0] == '>':
                spline = line.split()
                pos = int(re.findall(r'\d+', spline[1])[0])
                tosave = spline[0]+'-'+str(pos)
                saved.append(tosave)
                savedset.add(tosave)
        line = f.readline()
print len(saved), len(savedset)
f.close()

print savedset.difference(savedset_orig)

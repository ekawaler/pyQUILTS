f = open('orig/proteome.aa.fasta','r')

origmap = {}
for i in range(1,24):
	origmap[str(i)] = {}
origmap['MT'] = {}
origmap['X'] = {}
origmap['Y'] = {}

line = f.readline()
while line:
	chr = line.split('chr=')[1].split('|')[0]
	nextline = f.readline()
	origmap[chr][line] = nextline
	line = f.readline()

f.close()

f = open('new/proteome.aa.fasta','r')

newmap = {}
for i in range(1,24):
	newmap[str(i)] = {}
newmap['MT'] = {}
newmap['X'] = {}
newmap['Y'] = {}

line = f.readline()
while line:
	chr = line.split('chr=')[1].split('|')[0]
	nextline = f.readline()
	newmap[chr][line] = nextline
	line = f.readline()

f.close()

only_orig = open('only_orig.fasta','w')
only_new = open('only_new.fasta','w')
diff = open('diffs.fasta','w')

for i in origmap.keys():
	for j in origmap[i].keys():
		if j in newmap[i]:
			if origmap[i][j] != newmap[i][j]:
				diff.write(j.rstrip()+'-orig\n')
				diff.write(origmap[i][j])
				diff.write(j.rstrip()+'-new\n')
				diff.write(newmap[i][j])
		else:
			only_orig.write(j)
			only_orig.write(origmap[i][j])

for i in newmap.keys():
	for j in newmap[i].keys():
		if j in origmap[i]:
			pass
		else:
			only_new.write(j)
			only_new.write(newmap[i][j])

only_orig.close()
only_new.close()
diff.close()
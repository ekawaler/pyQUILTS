'''
Compares two different proteome.bed.var files (just at the gene level - does a gene show up in one but not the other?)
Not very exciting.
'''

import sys

'''f = open(sys.argv[1],'r')
keep1 = {}
for line in f.readlines():
	if line.find('-') != -1:
		keep1[line.split('\t')[0]] = line
f.close()

f = open(sys.argv[2],'r')
keep2 = {}
for line in f.readlines():
	if line.find('-') != -1:
		name = line.split('\t')[0]
		if name in keep1.keys():
			del keep1[name]
		else:
			keep2[name] = line
f.close()


w = open('compare_vars.txt','w')
w.write("Only in %s\n" % sys.argv[1])
for k in keep1:
	w.write(keep1[k]+'\n')
w.write("\nOnly in %s\n" % sys.argv[2])
for k in keep2:
	w.write(keep2[k]+'\n')
w.close()'''

# Now I want to check if the entire line is the same.
f = open(sys.argv[1],'r')
keep1 = set([])
for line in f.readlines():
	if line.find('-') != -1:
		keep1.add(line)
f.close()

f = open(sys.argv[2],'r')
keep2 = set([])
for line in f.readlines():
	if line.find('-') != -1:
		if line in keep1:
			keep1.discard(line)
		else:
			keep2.add(line)
f.close()


w = open('compare_vars.txt','w')
w.write("Only in %s\n" % sys.argv[1])
for k in keep1:
	w.write(k)
w.write("\nOnly in %s\n" % sys.argv[2])
for k in keep2:
	w.write(k)
w.close()
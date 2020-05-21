f = open('old_proteome.bed','r')
old_prots = f.readlines()
f.close()

count_same = 0
f = open('proteome.bed','r')
line = f.readline()
while line:
	spline = line.split('\t')
	spline[4] = '1000'
	line = '\t'.join(spline)
	if line in old_prots:
		count_same += 1
		old_prots.remove(line)
	line = f.readline()
print count_same
print len(old_prots)
f.close()

w = open('missing_old.txt','w')
for prot in old_prots:
	w.write(prot)
w.close()

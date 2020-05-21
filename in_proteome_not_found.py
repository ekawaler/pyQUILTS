only_cpdb = []
f = open('only_cpdb.txt','r')
for line in f.readlines():
	only_cpdb.append(line.rstrip())
f.close()

proteome = set([])
f = open('proteome.bed','r')
for line in f.readlines():
	proteome.add(line.split()[3])
f.close()

refseqali = []
f = open('Karsten/rsa.txt','r')
for line in f.readlines():
	refseqali.append(line.rstrip())
f.close()

in_both = []
for thing in only_cpdb:
	if thing in proteome:
		in_both.append(thing)

print in_both, len(in_both)

in_count = 0
not_in_refseqali = []
for thing in only_cpdb:
	if thing not in refseqali:
		not_in_refseqali.append(thing)
	else:
		in_count += 1
	
print not_in_refseqali, len(not_in_refseqali), in_count
print refseqali[:5]

f = open('newfr.txt','r')
from sets import Set
nms = Set([])
for line in f.readlines():
	#nm = line.split()[0].rsplit('-',1)[0]
	nm = line.split()[0]
	nms.add(nm)
f.close()

f = open('untrans.newfr.txt','r')
for line in f.readlines():
	if line.split()[0] not in nms:
		print line.split()[0]
f.close()

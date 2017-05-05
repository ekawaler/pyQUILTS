'''Various snippets for finding duplicates in the tryptic peptide output file.'''

import sys

thingmap = {}
f = open(sys.argv[1],'r')

headline = f.readline().split()[0]
seqline = f.readline()

while headline:
	if seqline in thingmap:
		if headline in thingmap[seqline]:
			print seqline, headline
		thingmap[seqline].append(headline)
	else:
		thingmap[seqline] = [headline]
	headline = f.readline().split()[0]
	seqline = f.readline()

if seqline in thingmap:
	if headline in thingmap[seqline]:
		print seqline, headline
	thingmap[seqline].append(headline)
else:
	thingmap[seqline] = [headline]

f.close()

#for k in thingmap.keys():
#	if len(thingmap[k]) > 1:
#		print k, thingmap[k]


'''
Found duplicates:
VTARTNK
['>NP_001138777 START:237 END:243 (missed cleavage after 240) VAR:S239T\n', '>NP_079525 START:283 END:289 (missed cleavage after 286) VAR:S285T\n']
IKTGPQK
['>NP_001138777 START:559 END:565 (missed cleavage after 560) VAR:T564P\n', '>NP_079525 START:605 END:611 (missed cleavage after 606) VAR:T610P\n']
VRGVAFLPHQTVTIR
['>NP_001138777 START:185 END:199 (missed cleavage after 186) VAR:L189V\n', '>NP_079525 START:231 END:245 (missed cleavage after 232) VAR:L235V\n']
AVVHPR
['>NP_001138777 START:858 END:863 VAR:P862H\n', '>NP_079525 START:904 END:909 VAR:P908H\n']
NPIIQEEENIFK
['>NP_001138550 START:59 END:70 VAR:V62I\n', '>NP_001018857 START:192 END:203 VAR:V195I\n']
GVAFLPHQTVTIRFPCPVSLDAK
['>NP_001138777 START:187 END:209 (missed cleavage after 199) VAR:L189V\n', '>NP_079525 START:233 END:255 (missed cleavage after 245) VAR:L235V\n']
NKAVVHPR
['>NP_001138777 START:856 END:863 (missed cleavage after 857) VAR:P862H\n', '>NP_079525 START:902 END:909 (missed cleavage after 903) VAR:P908H\n']
VAEEGDKPPHVFVPVDMAVTLPR
['>NP_001138777 START:590 END:612 VAR:Y602F\n', '>NP_079525 START:636 END:658 VAR:Y648F\n']
TGPQKQAK
['>NP_001138777 START:561 END:568 (missed cleavage after 565) VAR:T564P\n', '>NP_079525 START:607 END:614 (missed cleavage after 611) VAR:T610P\n']
AVVHPRR
['>NP_001138777 START:858 END:864 (missed cleavage after 863) VAR:P862H\n', '>NP_079525 START:904 END:910 (missed cleavage after 909) VAR:P908H\n']
GVAFLPHQTVTIR
['>NP_001138777 START:187 END:199 VAR:L189V\n', '>NP_079525 START:233 END:245 VAR:L235V\n']
NPIIQEEENIFKCNECEK
['>NP_001138550 START:59 END:76 (missed cleavage after 70) VAR:V62I\n', '>NP_001018857 START:192 END:209 (missed cleavage after 203) VAR:V195I\n']
'''
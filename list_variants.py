import sys

# The line to parse to get the variants looks like this:
# >NP_003234-S-G2512A (MAP:chr1:92149295- 119,108,42,121,300,159,141,153,338,190,148,169,184,138,185,61 0,11933,14350,24924,28504,32497,32829,35573,36154,38216,43920,46066,51037,74874,113548,177732 S-G2512A:158.000000) (VAR:S-S15F:158.000000)

# The file in the pyQUILTS results folder you'll want to look at is called proteome.aa.var.bed.dna.fasta
# I know it's an ugly name. It'll eventually get combined with some other fasta files to make the final output
# and that'll be called something nicer.

try:
	f = open(sys.argv[1],'r')
except IndexError:
	raise SystemExit("ERROR: Must specify an input file on command line:\n\tpython list_variants.py <input.fasta>")
except IOError:
	raise SystemExit("ERROR: Cannot open input file %s" % sys.argv[1])

w = open('variant_list.txt','w')

line = f.readline()
while line:
	if line[0] == '>':
		id = line[1:].split('-')[0]
		var = line.split('VAR:')[1].split(':')[0]
		w.write('%s\t%s\n' % (id, var))
	line = f.readline()

f.close()
w.close()
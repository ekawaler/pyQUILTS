f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/tryptic_proteome.fasta",'r')
# Quilts v.Python
# >NP_005808 START:50 END:56 VAR:S52G
python_quilts = set([])
for line in f.readlines():
	if line[0] == '>':
		python_quilts.add(line[1:].split()[0])
f.close()

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/proteome_new.fasta",'r')
# Quilts v.Python with new proteome.bed
# >NP_005808 START:50 END:56 VAR:S52G
new_quilts = set([])
for line in f.readlines():
        if line[0] == '>':
                new_quilts.add(line[1:].split()[0])
f.close()

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/customProDB results/var.flt.vcf.hg19_COLO829_chr.vcf_nsSNV.fasta",'r')
# customProDB
# >NP_000005_N639D,I1000V |NM_000014|A2M|chr12|-|9248233,9232268|alpha-2-macroglobulin precursor
custom_pro = set([])
for line in f.readlines():
	if line[0] == '>':
		custom_pro.add('_'.join(line[1:].split('_')[:2]))
f.close()

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/proteome_v2.fasta",'r')
# Quilts v.2.0
# >NP_005808-variant-stop87 perilipin-3 isoform 1 [Homo sapiens]. GN=PLIN3 (MAP:chr19:4859829- 195,66 0,1499 S-T1151C:222.0,) (VAR:S-I52V:222.0,)
quilts2 = []
for line in f.readlines():
	if line[0] == '>':
		quilts2.append(line[1:].split('-')[0])
f.close()

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/proteome.aa.var.bed.dna.fasta",'r')
# Quilts v.2.0
# >NP_005808-variant-stop87 perilipin-3 isoform 1 [Homo sapiens]. GN=PLIN3 (MAP:chr19:4859829- 195,66 0,1499 S-T1151C:222.0,) (VAR:S-I52V:222.0,)
quilts_py_full = set([])
for line in f.readlines():
	if line[0] == '>':
		quilts_py_full.add(line[1:].split('-')[0])
f.close()

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/proteome_refgene.fasta",'r')
# Quilts v.2.0
# >NP_005808-variant-stop87 perilipin-3 isoform 1 [Homo sapiens]. GN=PLIN3 (MAP:chr19:4859829- 195,66 0,1499 S-T1151C:222.0,) (VAR:S-I52V:222.0,)
quilts_py_full_new = set([])
for line in f.readlines():
        if line[0] == '>':
                quilts_py_full_new.add(line[1:].split('-')[0])
f.close()

print "In Python QUILTS: %d" % len(python_quilts)
print "In QUILTS 2.0: %d" % len(quilts2)
print "In customProDB: %d" % len(custom_pro)
print "In full Python QUILTS proteins: %d" % len(quilts_py_full)
print "In Python QUILTS, tryptic, with new proteome: %d" % len(new_quilts)
print "In Python QUILTS, full, with new proteome: %d" % len(quilts_py_full_new)

in_both = []
in_python = []
in_custom = []

for thing in quilts_py_full_new:
	if thing in custom_pro:
		in_both.append(thing)
	else:
		in_python.append(thing)
for thing in custom_pro:
	if thing not in quilts_py_full_new:
		in_custom.append(thing)

w = open('only_cpdb.txt','w')
for thing in in_custom:
	w.write(thing+'\n')
w.close()

print "Examples of the %d genes in both: %s" % (len(in_both),', '.join(in_both[:5]))
print "Examples of the %d genes in QUILTS only: %s" % (len(in_python),', '.join(in_python[:5]))
print "Examples of the %d genes in customProDB only: %s" % (len(in_custom),', '.join(in_custom[:5]))

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/proteome.aa.var.bed.dna.fasta",'r')
# PyQUILTS full proteins
# >NP_001171800-S-G391A (MAP:chr1:161677003+ 130,296 0,5871 S-G391A:210.000000) (VAR:S-R131Q:210.000000)
quilts_py_full = {}
for line in f.readlines():
	if line[0] == '>':
		quilts_py_full.setdefault(line[1:].split('-')[0],[]).append(line.split('-')[-1].split(':')[0])
f.close()

f = open("/Users/Emily/Desktop/Fenyo Lab/Quilts/Karsten/customProDB results/var.flt.vcf.hg19_COLO829_chr.vcf_nsSNV.fasta",'r')
# customProDB
# >NP_000005_N639D,I1000V |NM_000014|A2M|chr12|-|9248233,9232268|alpha-2-macroglobulin precursor
custom_pro = {}
for line in f.readlines():
	if line[0] == '>':
		key = '_'.join(line[1:].split('_')[:2])
		variants = line.split('_')[2].split()[0].split(',')
		custom_pro[key] = variants
		#custom_pro.setdefault(key,[]).append(line.split()[-1].split(':')[0])
		#custom_pro.add('_'.join(line[1:].split('_')[:2]))
f.close()


same_variants = 0
sv = []
different_variants = 0
dv = []
more_in_cpdb = 0
mc = []
more_in_quilts = 0
mq = []
just_different = 0
jd = []

for i in quilts_py_full.keys():
	if i in custom_pro.keys():
		custom_pro[i].sort()
		quilts_py_full[i].sort()
		if custom_pro[i] == quilts_py_full[i]:
			same_variants += 1
			sv.append(i)
		else:
			different_variants += 1
			dv.append(i)
			#print i, custom_pro[i], quilts_py_full[i]
			if len(custom_pro[i]) > len(quilts_py_full[i]):
				more_in_cpdb += 1
				mc.append(i)
			elif len(custom_pro[i]) < len(quilts_py_full[i]):
				more_in_quilts += 1
				mq.append(i)
			else:
				just_different += 1
				jd.append(i)

for i in ['NP_056052', 'NP_061187', 'NP_060686', 'NP_112214', 'NP_775922', 'NP_945352']:
	print custom_pro[i], quilts_py_full[i]

print "Exactly the same: %d, %s" % (same_variants, ', '.join(sv[:5]))
print "Different: %d" % different_variants
print "More in cpDB: %d, %s" % (more_in_cpdb, ', '.join(mc[:5]))
print "More in QUILTS: %d, %s" % (more_in_quilts, ', '.join(mq[:5]))
print "Different variants: %d, %s" % (just_different, ', '.join(jd[:6]))

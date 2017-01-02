'''Trying to find a WHIM that has overlapping somatic/germline variants so I can test my merging script.'''

whims = [2,4,5,6,8,9,11,12,13,14,16,17,18,20,21,24,25,26,27,30,32,35,37,40,43,44,46]

for i in whims:
	somatic_variants = set([])
	folder = "/ifs/data/proteomics/tcga/samples/breast-xenografts/whim"+str(i)+"/rna/vcf"
	f = open(folder+'/somatic.vcf','r')
	line = f.readline()
	line = f.readline()
	while line:
		spline = line.split('\t')
		somatic_variants.add(spline[0]+'#'+spline[1])
		line = f.readline()
	f.close()
	
	f = open(folder+'/germline.vcf')
	line = f.readline()
	line = f.readline()
	line = f.readline()
	while line:
		spline = line.split('\t')
		if spline[0]+'#'+spline[1] in somatic_variants:
			print("WHIM%d: %s, %s" % (i, spline[0], spline[1]))
		line = f.readline()
	f.close()
	
'''Cool, my output: WHIM27: chr9, 33798017
WHIM27: chr12, 11461596
WHIM27: chr13, 111109289
WHIM30: chr11, 1017110
WHIM30: chr12, 11420454
WHIM32: chr1, 16890331
WHIM32: chr1, 16946437
WHIM35: chr4, 88536883
WHIM37: chr1, 16958024
WHIM37: chr6, 30954934
WHIM37: chr7, 151970859
WHIM37: chr20, 17639807
WHIM37: chr21, 11039055
WHIM40: chr15, 23685126
WHIM43: chr1, 16959834
WHIM43: chr16, 90095620
WHIM43: chr17, 21319007
WHIM44: chr1, 144994861
WHIM44: chr6, 29910717
WHIM44: chr9, 33798017'''
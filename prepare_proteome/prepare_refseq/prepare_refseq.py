import argparse
import shutil

parser = argparse.ArgumentParser(description="Prepare Proteome Reference")


f = open('refLink.txt','r')
genes = open('proteome-genes.txt','w')
descriptions = open('proteome-descriptions.txt','w')
trans_map = {}
for line in f.readlines():
	spline = line.split('\t')
	if spline[3] != '':
		trans_map[spline[2]] = spline[3]
		genes.write('%s\t%s\n' % (spline[3], spline[0].upper()))
		descriptions.write('%s\t%s\n' % (spline[3],spline[1]))
f.close()
genes.close()
descriptions.close()

f = open('refSeq.bed','r')
w = open('proteome.bed','w')

for line in f.readlines():
	spline = line.rstrip().split('\t')
	if spline[3] in trans_map.keys() and len(spline[0].split('_')) == 1:
		spline[3] = trans_map[spline[3]]
		# We'll only keep the exons that are in the CDS (between thickStart and thickEnd, I guess?
		# yeah, they're the start and stop codons), which is to say between spline[6] and spline[7]
		# So the fields that will be changed by this are: exonCount (9), exonLengths (10), and exonStarts (11)
		# If thickStart = thickEnd, it's a ncRNA, and we can dump it.
		if spline[6] == spline[7]:
			continue
		cdsStart, cdsEnd = int(spline[6]), int(spline[7])
		exonLengths = spline[10].rstrip(',').split(',')
		exonStarts = spline[11].rstrip(',').split(',')
		startDiff = cdsStart-int(spline[1])
		cdsLength = cdsEnd-cdsStart
		newExonStarts = [int(x)-startDiff for x in exonStarts]
		newExonLengths = [int(x) for x in exonLengths]
		for i in range(len(newExonStarts)-1,-1,-1):
			if newExonStarts[i] > cdsLength:
				del newExonStarts[i]
				del newExonLengths[i]
			elif newExonStarts[i] < 0 and newExonStarts[i]+newExonLengths[i] > cdsLength:
				newExonStarts[i] = 0
				newExonLengths[i] = cdsLength
			elif newExonStarts[i] < 0:
				if newExonStarts[i]+newExonLengths[i] <= 0:
					del newExonStarts[i]
					del newExonLengths[i]
				else:
					newExonLengths[i] = newExonLengths[i] + newExonStarts[i]
					newExonStarts[i] = 0
			elif newExonStarts[i]+newExonLengths[i] > cdsLength:
				newExonLengths[i] = cdsLength-newExonStarts[i]

		# Let's now change the line and write it out
		spline[1], spline[2] = spline[6], spline[7]
		newExonLengths = [str(x) for x in newExonLengths]
		newExonStarts = [str(x) for x in newExonStarts]
		spline[9] = str(len(newExonLengths))
		spline[10] = ','.join(newExonLengths)
		spline[11] = ','.join(newExonStarts)
		w.write('\t'.join(spline)+'\n')

f.close()
w.close()

# I'm using "copy" here instead of "move" in case someone wants to run this again and 
# doesn't remember where their refSeq.bed went.
shutil.copy('refSeq.bed','transcriptome.bed')
import argparse
import shutil
import sys
import os

parser = argparse.ArgumentParser(description="Prepare Proteome Reference")

parser.add_argument('--uniprot_bed', type=str, help="Name of your UniProt .bed file (step 1A in the README)", required=True)
parser.add_argument('--uniprot_fasta', type=str, help="Name of your .fasta file (step 1B in the README)", required=True)
parser.add_argument('--ucsc_bed', type=str, help="Name of your UCSC .bed file (step 1C in the README)", required=True)
parser.add_argument('--version', type=str, default="swiss", help="Version of UniProt to use (either 'swiss', 'trembl', or 'both')", choices=['swiss','trembl','both'], required=True)

args = parser.parse_args()

# Make proteome.bed (this is just the ucsc bed file but with the identifiers fixed to match the other files)
f = open(args.ucsc_bed,'r')
w = open('tmp_proteome.bed','w') # just in case they named the original proteome.bed
for line in f.readlines():
	spline = line.split('\t')
	spline[3] = spline[3].split('-')[0]
	w.write('\t'.join(spline))
f.close()
w.close()
os.system('mv %s ucsc_proteome.bed' % args.ucsc_bed)
os.system('mv tmp_proteome.bed proteome.bed')

# Make proteome.fasta (this is just the proteome fasta file)
if args.version == 'both':
	os.system('cp %s proteome.fasta' % args.uniprot_fasta)
else:
	kept = []
	if args.version == 'swiss':
		starter = '>sp'
	else:
		starter = '>tr'
	w = open('tmp_proteome.fasta','w')
	f = open(args.uniprot_fasta,'r')
	line = f.readline()
	while line:
		if line.split('|')[0] == starter:
			kept.append(line.split('|')[1])
			w.write(line)
			line = f.readline()
			while line and line[0] != '>':
				w.write(line)
				line = f.readline()
		else:
			line = f.readline()
			while line and line[0] != '>':
				line = f.readline()
	w.close()
	f.close()
	os.system('mv %s uniprot_proteome.fasta' % args.uniprot_fasta)
	os.system('mv tmp_proteome.fasta proteome.fasta')

# Make transcriptome.bed (this is the uniprot bed file with all the swissprot/trembl entries removed depending on which version you're using)
if args.version == 'both':
	os.system('cp %s transcriptome.bed' % args.uniprot_bed) # copying in case they want to run it again and wonder where their uniprot bed file went
else:
	w = open('transcriptome.bed','w')
	f = open(args.uniprot_bed, 'r')
	for line in f.readlines():
		if line.split()[3] in kept:
			w.write(line)
	w.close()
	f.close()

# Make proteome-genes.txt and proteome-descriptions.txt
f = open('proteome.fasta','r')
genes = open('proteome-genes.txt','w')
desc = open('proteome-descriptions.txt','w')

line = f.readline()
while line:
	if line[0] == '>':
		g_id = line.split('|')[1]
		if 'GN=' in line:
			g_name = line.split('GN=')[1].split('PE=')[0].rstrip()
		else:
			g_name = ''
		g_desc = line.split(' ',1)[1].split('OS=')[0].rstrip()
		genes.write('%s\t%s\n' % (g_id, g_name))
		desc.write('%s\t%s\n' % (g_id, g_desc))
		line = f.readline()
	else:
		line = f.readline()
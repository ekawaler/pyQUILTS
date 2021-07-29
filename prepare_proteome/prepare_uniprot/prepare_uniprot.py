import argparse
import os

parser = argparse.ArgumentParser(description="Prepare Proteome Reference")

parser.add_argument('--bed_file', type=str, help="Name of your UniProt .bed file (step 1A in the README)", required=True)
parser.add_argument('--input_fasta', type=str, help="Name of your .fasta file (step 1B in the README)", required=True)
parser.add_argument('--info_file', type=str, help="Name of your .tab info file (step 1C in the README)", required=True)
parser.add_argument('--version', type=str, default="swiss", help="Version of UniProt to use (either 'swiss', 'trembl', or 'both'; defaults to 'swiss')", choices=['swiss','trembl','both'])

args = parser.parse_args()

# Make proteome.fasta (this is just the proteome fasta file)
if args.version == 'both':
	os.system('cp %s proteome.fasta' % args.input_fasta)
	kept = []
	for line in open("proteome.fasta"):
		if line.startswith('>'):
			kept.append(line.split('|')[1])
else:
	kept = []
	if args.version == 'swiss':
		starter = '>sp'
	else:
		starter = '>tr'
	w = open('tmp_proteome.fasta','w')
	f = open(args.input_fasta,'r')
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
	os.system('mv tmp_proteome.fasta proteome.fasta')

kept = set(kept) # for speeding up the lookups below
# Make transcriptome.bed and proteome.bed (transcriptome.bed is the uniprot bed file with all the swissprot/trembl entries removed depending on which version you're using, proteome.bed is that but untranslated sections are removed)
if True: # only so that the diff for pull request is easier
	w = open('transcriptome.bed','w')
	p = open('proteome.bed','w')
	f = open(args.bed_file, 'r')
	for line in f.readlines():
		if line.split()[3].split('-')[0] in kept: # Keeping only SwissProt/TrEMBL
			# Basic stuff for transcriptome file
			spline = line.split()
			spline[3] = spline[3].replace('-','@')
			w.write('\t'.join(spline[:12])+'\n')
			
			# Fancier stuff for proteome file
			pstart = int(spline[6])
			pend = int(spline[7])
			tstart = int(spline[1])
			tend = int(spline[2])
			lens = []
			starts = []
			ex_lens = list(map(int, spline[10].split(',')))
			ex_starts = list(map(int, spline[11].split(',')))
			pflag = False
			for i in range(len(ex_lens)):
				if not pflag:
					# check if it's a start and then set pflag
					if pstart <= tstart+ex_starts[i]+ex_lens[i]:
						pflag = True
						lens.append(tstart+ex_starts[i]+ex_lens[i]-pstart)
						starts.append(0)
				elif pflag:
					# Check for end, break
					if pend <= tstart+ex_starts[i]+ex_lens[i]:
						lens.append(pend-tstart-ex_starts[i])
						starts.append(tstart+ex_starts[i]-pstart)
						break
					# No end? Just add the thing and change the start
					else:
						lens.append(ex_lens[i])
						starts.append(tstart+ex_starts[i]-pstart)
			spline[1] = str(pstart)
			spline[2] = str(pend)
			spline[9] = str(len(lens))
			spline[10] = ','.join(map(str,lens))
			spline[11] = ','.join(map(str,starts))			
			p.write('\t'.join(spline[:12])+'\n')
			
	w.close()
	p.close()
	f.close()

# Make proteome-genes.txt and proteome-descriptions.txt
f = open(args.info_file,'r')
g = open('proteome-genes.txt','w')
desc = open('proteome-descriptions.txt','w')

header = f.readline().rstrip().split('\t')
entry_pos = header.index("Entry")
prot_pos = header.index("Protein names")
gene_pos = header.index("Gene names")
#print header
genes = {}

line = f.readline()
while line:
	spline = line.split('\t')
	uid = spline[entry_pos].replace('-','@')
	pn = spline[prot_pos]
	gn = spline[gene_pos]
	if gn != '\n':
		gn = gn.split()[0]
	else:
		gn = "NOSYMBOL"
	g.write("%s\t%s\t%s\n" % (uid, uid, gn))
	desc.write("%s\t%s\n" % (uid, pn))
	line = f.readline()

f.close()
g.close()
desc.close()

# proteome-genes.txt: ENSP00000000233 ENST00000000233 ENSG00000004059 ARF5
# proteome-descriptions.txt: ENSP00000000233-ARF5    ADP ribosylation factor 5 [Source:HGNC Symbol;Acc:HGNC:658]
# proteome.bed: chr7	127588498	127591299	ENSP00000000233	1000	+	127588498	127591299	0	6	67,81,110,72,126,87	0,584,986,1567,2464,2714
# proteome.bed: chr1	924431	939291	ENSP00000411579	1000	+	924431	939291	0	7	517,92,182,51,125,90,17	0,1490,5723,6607,11340,14608,14843

import sys

class Gene:
	def __init__(self, chromosome, gene_id, protein_id, transcript_id, gene_name, gene_desc, strand):
		self.chromosome = 'chr'+chromosome
		self.gene_id = gene_id
		self.protein_id = protein_id
		self.transcript_id = transcript_id
		self.gene_name = gene_name
		self.gene_desc = gene_desc
		self.strand = '+' if strand=='1' else '-'
		self.exons = {}
	
	def add_exon(self, exon_no, exon_start, exon_end):
		self.exons[int(exon_no)] = [int(exon_start), int(exon_end)]
	
	def print_bed_line(self):
		exon_nos = self.exons.keys()
		if self.strand == '+':
			exon_nos.sort()
			start = self.exons[min(exon_nos)][0]
			end = self.exons[max(exon_nos)][1]
		else:
			exon_nos.sort()
			exon_nos.reverse()
			end = self.exons[min(exon_nos)][1]
			start = self.exons[max(exon_nos)][0]
		lengths = []
		starts = []
		for i in exon_nos:
			ex_start, ex_end = self.exons[i][0], self.exons[i][1]
			lengths.append(str(ex_end-ex_start+1))
			starts.append(str(ex_start-start))
		print_str = "%s\t%d\t%d\t%s\t1000\t%s\t%d\t%d\t0\t%d\t%s\t%s\n" % (self.chromosome, start, end, self.protein_id, self.strand, start, end, len(self.exons), ','.join(lengths), ','.join(starts))
		# Negative strand!?
		return print_str
			
f = open(sys.argv[1], 'r')
header = f.readline().rstrip().split('\t')
chrm_pos = header.index("Chromosome/scaffold name")
gene_id_pos = header.index("Gene stable ID")
transcript_id_pos = header.index("Transcript stable ID")
protein_id_pos = header.index("Protein stable ID")
gene_name_pos = header.index("Gene name")
gene_desc_pos = header.index("Gene description")
strand_pos = header.index("Strand")
exon_no_pos = header.index("Exon rank in transcript")
exon_start_pos = header.index("Exon region start (bp)")
exon_end_pos = header.index("Exon region end (bp)")
cds_start_pos = header.index("Genomic coding start")
cds_end_pos = header.index("Genomic coding end")
print header
genes = {}

line = f.readline()
while line:
	spline = line.rstrip().split('\t')
	if spline[protein_id_pos] != '' and spline[cds_start_pos] != '':
		gene_name = spline[protein_id_pos]
		if gene_name not in genes:
			new_gene = Gene(spline[chrm_pos], spline[gene_id_pos], spline[protein_id_pos], spline[transcript_id_pos], spline[gene_name_pos], spline[gene_desc_pos], spline[strand_pos])
			genes[gene_name] = new_gene
		genes[gene_name].add_exon(spline[exon_no_pos], spline[cds_start_pos], spline[cds_end_pos])
	line = f.readline()	
f.close()

for i in genes.keys()[:5]:
	print i,genes[i].exons

wgenes = open('proteome-genes.txt','w')
wdesc = open('proteome-descriptions.txt','w')
wbed = open('proteome.bed','w')
for g in genes:
	gene = genes[g]
	if gene.exons != {}:
		wgenes.write('%s\t%s\t%s\t%s\n' % (gene.protein_id,gene.transcript_id,gene.gene_id,gene.gene_name))
		wdesc.write('%s\t%s\n' % (gene.protein_id, gene.gene_desc))
		wbed.write(gene.print_bed_line())
wgenes.close()
wdesc.close()
wbed.close()
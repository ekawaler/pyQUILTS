You need three files:

1. refLink.txt from http://hgdownload.cse.ucsc.edu/goldenpath/hgFixed/database/refLink.txt.gz

2. refSeq.bed from https://genome.ucsc.edu/cgi-bin/hgTables
-Genome: Human, group: Genes & Gene Predictions, Track: RefSeq Genes, table: refGene, output format: BED, then on the next page create one BED record per whole gene

3. proteome.fasta from https://genome.ucsc.edu/cgi-bin/hgTables
-Genome: Human, group: Genes & Gene Predictions, Track: RefSeq Genes, table: refGene, output format: sequence, then on the next page: proteins
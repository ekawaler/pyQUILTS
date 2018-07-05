1. You need three files:

	A. The UniProt proteome .bed file. You can find this at ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/
		-Enter UPXXXXXXXXX_XXXX_beds folder, and download UPXXXXXXXXX_XXXX_proteome.bed (human only)

	B. A proteome .fasta file from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/
		-Download the file UPXXXXXXXXX_XXXX.fasta.gz, where the UP number matches the one from the .bed file you downloaded in step 1A. Unzip this.

	C. A .bed file from the UCSC genome browser at https://genome.ucsc.edu/cgi-bin/hgTables
		-Genome: Human, group: Genes & Gene Predictions, Track: UniProt, table: SwissProt Aln. or TrEMBL Aln., depending on your needs, output format: BED, then on the next page create one BED record per whole gene
		-If you want to use the combined SwissProt and TrEMBL, I think you can probably download both and concatenate them. I haven’t tested this, though.

2. Next, you’ll run the python script:

	python prepare_uniprot.py --uniprot_bed [uniprot bed file name] --uniprot_fasta [uniprot fasta file name] --ucsc_bed [ucsc bed file name] --version [‘swiss’,’trembl’, or ‘both’]

This will leave you with the five reference files you need to run QUILTS: proteome.bed, transcriptome.bed, proteome.fasta, proteome-genes.txt, proteome-description.txt
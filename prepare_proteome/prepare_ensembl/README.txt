
1. Get your info files from Biomart.

	Go to: http://useast.ensembl.org/biomart/martview/

	Choose the correct database and organism from the dropdown menus.

	In the sidebar, under “Attributes”, click the “Structures” radio button and then select the following boxes:
	-Under GENE: Gene stable ID, Transcript stable ID, Protein stable ID, Gene name, Gene description, Chromosome/Scaffold name, Strand
	-Under EXON: Exon region start (bp), Exon region end (bp), Exon rank in transcript, Genomic coding start, Genomic coding end

	Hit “Results” on the top left, then export your results to a TSV.


2. Get your reference proteome FASTA file.

	Go to: http://useast.ensembl.org/info/data/ftp/index.html/

	In the “Single species data” table, find your species and click the FASTA link in the column “Protein sequence (FASTA)”

	Download and unzip the file with the extension .pep.all.fa.gz

	At least for the human and mouse proteomes, this only works with release 84 onward. (They added some new data to the FASTA header that we use.)

3. Run the prepare_ensembl.py script:

	python prepare_ensembl.py <the file you downloaded in step 1> <the .pep.all.fa file>

This should output four files: transcriptome.bed, proteome.bed, proteome-genes.txt, proteome-descriptions.txt. 

4. Once this is done, rename your .pep.all.fa file to “proteome.fasta” (the file QUILTS will be looking for). Congratulations, you’ve got your reference proteome ready for QUILTS!

1. Get your chromosome data files.

	Go to: http://useast.ensembl.org/info/data/ftp/index.html/

	In the “Single species data” table, find your species and click the GenBank link in the column “Annotated sequence (GenBank)”

	Download all of the .dat.gz files here (you can use wget from the command line to do this more efficiently, for example: “wget ftp://ftp.ensembl.org/pub/release-92/genbank/mus_musculus/*”)

	Unzip the files (“gunzip *.dat.gz”)

2. Get your reference proteome FASTA file.

	Go to: http://useast.ensembl.org/info/data/ftp/index.html/

	In the “Single species data” table, find your species and click the FASTA link in the column “Protein sequence (FASTA)”

	Download and unzip the file with the extension .pep.all.fa.gz

	At least for the human and mouse proteomes, this only works with release 84 onward. (They added some new data to the FASTA header that we use.)

3. Run the prepare_ensembl.py script:

	python prepare_ensembl.py <directory with your .dat files> <the .pep.all.fa file>

This should output four files: transcriptome.bed, proteome.bed, proteome-genes.txt, proteome-descriptions.txt. 

4. Once this is done, rename your .pep.all.fa file to “proteome.fasta” (the file QUILTS will be looking for). Congratulations, you’ve got your reference proteome ready for QUILTS!
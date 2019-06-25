Welcome to the alpha version of QUILTS 3.0! This version of QUILTS is able to find variants (both single amino acid variants and indels) and most forms of splice junctions, though the junction functionality has not been extensively tested, so use at your own risk.

TO RUN:
The required arguments are ‘output_dir’, ’genome’, and ‘proteome’. Also must have one or more of ‘somatic’, ‘germline’, ‘fusion’, and ‘junction’. An example genome and proteome (HG19, RefSeq from September 22, 2017) can be downloaded from https://nyumc.box.com/v/quilts-sample-references (please let me know if this link doesn’t work). To create your own refSeq proteome, go into the prepare_proteome folder and follow the instructions in the README there. To create your own genome, go into the prepare_genome folder and follow the instructions in the README there. You can also look at the sample running script included as runtest.sh (though it won’t run until you change the paths to reflect where the files are located on your system).

python quilts.py
	—-output_dir (path to the desired location of the results)
	—-genome (path to a folder containing the reference genome)
	—-proteome (path to a folder containing the reference proteome)
	—-somatic (path to a folder containing .vcf files of somatic mutations)
	—-germline (path to a folder containing .vcf files of germline mutations)
	—-fusion (path to a folder containing .txt files of fusion genes, format specified below)
	—-junction (path to a folder containing junction files, .bed if tophat, .txt if MapSplice, or .tab if STAR)
	—-junction_file_type (can be either ‘mapsplice’ (default), ‘tophat’, or ‘star’)
	—-threshA (integer - minimum number of reads to support a splice junction with conserved exon boundaries)
	—-threshAN (integer - minimum number of reads to support a splice junction with one conserved exon boundary)
	—-threshN (integer - minimum number of reads to support a completely novel splice junction)
	—-variant-quality-threshold (quality threshold for variants)
	—-no-missed-cleavage (defaults to allow for one missed cleavage; used for generating tryptic peptides which is not yet ready)

The output will be in a folder named results_(date).(time). The ‘log’ folder contains a lot of files that are generated in the process of running QUILTS, but what you’re looking for will be in the ‘fasta’ folder. ‘fasta/proteome.fasta’ combines all of the variants with the reference proteome and is what you’ll most likely want to use for your search. ‘fasta/parts/‘ contains a fasta file for each of the variant types (SAAV, indel, each splice junction type) - these files are what got combined to create your proteome.fasta file.

Fusion file input format: For the moment, since there are so many possibilities, I’m using the format we’re using in our projects. This will require your first five columns to be:
FusionName	LeftBreakpoint	RightBreakpoint	Sample	JunctionReadCount
GENE1—GENE2	chr19:1000000:-	chr16:500003:+	Sample1	5

You can threshold on JunctionReadCount, though that’s not an option - you have to go into the code to do it. In addition, you’re welcome to have whatever other columns you’d like after the first five.

An explainer for the headers in the output fasta can be found in OUTPUT FASTA HEADER KEY.txt
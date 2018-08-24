'''
quilts.py

Remaking QUILTS in Python.

It's in its early stages, so half the lines are just me asking myself questions
that I will forget to answer in the future.

Also, trying to be super meticulous about checking input for correctness/existence.
Hold me to this in the future when I start to get lazy.

A possible Important Thing for the future: write something to status/log whenever SystemExit happens?

A definite thing for the future: making a script that will call this one, since these filepaths
are real unwieldy. Should be able to move various "does this file actually exist" checks to that script.

Another definite thing: change all merged_pytest to merged

Complications to watch out for when defining tryptic peptides: right now we assume the missed cleavages will account for substitutions that add/remove an R/K. Might be worth checking if this actually happens.

Variant Start Codons: Right now, if there's a SNP in the start codon, we just throw that out. At some point, perhaps we should go through and look for the next possible start codon, assuming the transcription will start from there, but honestly, that is probably just going to make a garbage protein, so we're ignoring it for the time being.

Right now, there are no proteins/peptides that have both an indel and some other variant, unlike the SAAVs where we have tryptic peptides with all possible combinations of variants. This will probably change in the future, but for now if you've got an indel that's all you get.

Am I ever going to get rid of duplicate tryptic peptides (peptides that show up in multiple proteins)? Maybe that can be a postprocessing step. Or maybe that's a search problem. There aren't very many.

Emily Kawaler
'''

import argparse
import os
import sys
import shutil
from datetime import datetime
from subprocess import call, check_call, CalledProcessError
import warnings
from exonSearchTree import ExonSearchTree
from string import maketrans
import re
from sets import Set
from itertools import product

# ahhhh look at these hideous global variables
global logfile
global statusfile
global results_folder
global codon_map

### These functions are used in the setup phase.

def parse_input_arguments():
	''' Sets up an argument parser and checks a couple of the inputs for correctness.'''
	
	# Set up the argument parser
	# Do I want to have some sort of root directory so they don't have to enter full paths? Probably not - 
	# full paths are a pain but they're more flexible.
	# Also, at some point, make some of these arguments not optional.
	parser = argparse.ArgumentParser(description="QUILTS") # What even is this description for, anyway? Maybe I'll remove it eventually
	#parser.add_argument('--output_dir', type=str, default="/ifs/data/proteomics/tcga/scripts/quilts/pyquilts", help="full path to output folder")
	#parser.add_argument('--proteome', type=str, default="/ifs/data/proteomics/tcga/databases/refseq_human_20160914", help="full path to folder containing reference proteome")
	#parser.add_argument('--genome', type=str, default="/ifs/data/proteomics/tcga/databases/genome_human", help="full path to folder containing reference genome")
	parser.add_argument('--output_dir', type=str, default=".", help="full path to output folder")
	parser.add_argument('--proteome', type=str, default=".", help="full path to folder containing reference proteome")
	parser.add_argument('--genome', type=str, default=".", help="full path to folder containing reference genome")
	#parser.add_argument('--somatic', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna/vcf/TCGA-20130502-S", help="VCF file of somatic variants")
	#parser.add_argument('--germline', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna-germline/vcf/GATK26-G-Nature", help="VCF file of germline variants")
	#parser.add_argument('--junction', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/rna/tophat/junction_tophat208bowtie2_gMx", help="BED file of splice junctions [in progress]")
	parser.add_argument('--germline',type=str)
	parser.add_argument('--somatic',type=str)
	parser.add_argument('--junction', type=str)
	parser.add_argument('--mapsplice', action='store_true', default=False, help="Use this if your junctions were created by MapSplice rather than Tophat. Automatically converts the junctions.txt file at the junction location to a file called junctions.bed which can be found by this script.")
	parser.add_argument('--fusion', type=str, help="full path to folder containing fusion file")
	parser.add_argument('--threshA', type=int, default=2)
	parser.add_argument('--threshAN', type=int, default=3)
	parser.add_argument('--threshN', type=int, default=3)
	parser.add_argument('--variant_quality_threshold', type=float, default=0.0, help="Quality threshold for variants")
	parser.add_argument('--no_missed_cleavage', action='store_true', default=False, help="Tryptic peptide fasta by default allows for a single missed cleavage; adding this argument will tell the virtual trypsinizer to assume perfect cleavage")
	
	# Pull out the arguments
	args = parser.parse_args()
	# Check that we have a somatic and/or germline and/or junction file. Abort if not.
	# Will add fusion to this check later.
	if not args.somatic and not args.germline and not args.junction and not args.fusion:
		raise SystemExit("ERROR: Must have at least one variant file (somatic, germline, fusion, and/or junctions).\nAborting program.")
	
	return args

def set_up_output_dir(output_dir, args):
	''' Sets up the results folder (creates it, sets up a log file and status file).
		Because this will probably confuse people (me) later, "output directory" refers
		to the directory the user specifies in which the "results folder" will reside. The
		results folder contains the actual results.'''
		
	# I hate this global variable stuff but is there another way to do it
	# without having to pass the logfile/statusfile to all the functions so they can be used?
	global logfile
	global statusfile
	global results_folder
	
	# Will the directory they give us be where we put the results, or where we put a folder
	# containing the results?
	
	# Checks to make sure the output directory exists.
	# Currently aborts if this is the case, but I could also just create the folder 
	# in the working directory? Maybe allow user input for this? Probably not though.
	# Wouldn't be good if they were running it through qsub or whatever.
	if not os.path.isdir(output_dir):
		#raise SystemExit("ERROR: Output directory does not exist.\nAborting program.")
		os.makedirs(output_dir)
	
	# Makes the results folder. Calls it results_(date and time). Looks ugly, but I'm cool with that.
	# Lets us avoid the problem of dealing with multiple results folders in the same output directory.
	today = datetime.today()
	day_time_string = str(today.year)+str(today.month).zfill(2)+str(today.day).zfill(2)+'.'+str(today.hour).zfill(2)+str(today.minute).zfill(2)+str(today.second).zfill(2)
	results_folder = output_dir+'/results_'+day_time_string
	os.makedirs(results_folder)
	
	# Gives us addresses for the logfile and statusfile,
	# and writes the starting date and time to them.
	logfile = results_folder+'/log.txt'	
	statusfile = results_folder+'/status.txt'
	write_to_log("Logfile created: "+str(today), logfile)
	write_to_status("Status file created: "+str(today))
	
	# Creates some folders within the results folder
	# Currently just copying what's in the Perl version.
	# Likely to change as I start building out other functionality.
	os.makedirs(results_folder+"/log")
	os.makedirs(results_folder+"/fasta")
	os.makedirs(results_folder+"/fasta/parts")
	
	# Moves reference proteome to the working area.
	try:
		#shutil.copy(args.proteome+"/proteome.fasta",results_folder+"/fasta/"+args.proteome.split("/")[-1]+".fasta")
		shutil.copy(args.proteome+"/proteome.fasta",results_folder+"/fasta/reference_proteome.fasta")
	except IOError:
		raise SystemExit("ERROR: Reference proteome .fasta file not found at "+args.proteome+"/proteome.fasta.\nAborting program.")	
	try:
		shutil.copy(args.proteome+"/proteome.bed",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome .bed file not found at "+args.proteome+"/proteome.bed.\nAborting program.")
	# Not sure if these four should be quite so important.
	try:
		shutil.copy(args.proteome+"/proteome-descriptions.txt",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome gene descriptions file not found at "+args.proteome+"/proteome-descriptions.txt.\nAborting program.")
	try:
		shutil.copy(args.proteome+"/proteome-genes.txt",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome gene names file not found at "+args.proteome+"/proteome-genes.txt.\nAborting program.")
	if args.junction:
		try:
			shutil.copy(args.proteome+"/proteome-descriptions.txt",results_folder+"/log/alternative-descriptions.txt")
		except IOError:
			raise SystemExit("ERROR: Reference proteome gene descriptions file not found at "+args.proteome+"/proteome-descriptions.txt.\nAborting program.")
		try:
			shutil.copy(args.proteome+"/proteome-genes.txt",results_folder+"/log/alternative-genes.txt")
		except IOError:
			raise SystemExit("ERROR: Reference proteome gene names file not found at "+args.proteome+"/proteome-genes.txt.\nAborting program.")
	
	return results_folder

### These functions are for merging and quality checking the variant files.

def set_up_vcf_output_dir(vcf_dir):
	'''Make the directory we'll put our merged VCF files in. Returns with an error if directory doesn't exist.'''
	if not os.path.isdir(vcf_dir):
		# No VCF files are going to be found. Gotta leave the function.
		return 1
	try:
		# Change to merged later, but for now...
		os.makedirs(vcf_dir+"/merged_pytest")
	except OSError:
		# Right now I'm just going to overwrite things in there.
		warnings.warn("VCF directory "+vcf_dir+"/merged already exists!\nOverwriting contents...")
	return None

def pull_vcf_files(vcf_dir):
	'''Finds all .vcf files in the directory.'''
	files = os.listdir(vcf_dir)
	vcf_files = []
	for f in files:
		if f.endswith('.vcf'):
			vcf_files.append(f)
	return vcf_files	

def merge_and_qual_filter(vcf_dir, quality_threshold):
	'''Returns None if successful, or 1 if unable to find variant files. If duplicates appear, pick the one with the highest quality score.'''
	
	# Set up the output directory, if possible. If not found, return with a warning.
	output_dir_flag = set_up_vcf_output_dir(vcf_dir)
	if output_dir_flag:
		warnings.warn("VCF directory not found at %s" % vcf_dir)
		return 1

	# Pull the list of .vcf files in the directory. If none found, return with a warning.
	vcf_files = pull_vcf_files(vcf_dir)
	if len(vcf_files) == 0:
		warnings.warn("Unable to find any .vcf files in %s" % vcf_dir)
		return 1

	# Cool, we have a valid file, let's open our log and output files.
	vcf_log_location = vcf_dir+'/merged_pytest/merged.log'
	vcf_log = open(vcf_log_location,'w')
	w = open(vcf_dir+'/merged_pytest/merged.vcf','w')
	w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")	
	

	# Set up some tracking stats etc.
	qual_removed_all = 0
	total_variants_all = 0
	kept_variants_all = 0
	duplicates_all = 0
	#header_written = False
	existing_variants = {}

	for vf in vcf_files:
		f = open(vcf_dir+'/'+vf, 'r')
		write_to_log(vcf_dir+'/'+vf, vcf_log_location)
		
		# I read one line at a time instead of going with f.readlines() because
		# if the file is long, it's bad to hold all of the lines in memory.
		# In case you were wondering, which you almost certainly weren't.
		line = f.readline()
		qual_removed = 0
		total_variants = 0
		kept_variants = 0
		duplicates_removed = 0 # Duplicates in this file that were lower-quality than older duplicates
		old_duplicates_removed = 0 # Duplicates in this file that were higher-quality than older duplicates
		while line:
			spline = line.split('\t')
			if line[0] == '#':
				# It's a header line. Ignore it.
				#if not header_written:
				#	w.write('\t'.join(spline[0:6])+'\n')
				#	header_written = True
				line = f.readline()
				continue
			if len(spline) < 6:
				# Bummer, we don't have six tab-separated fields, this isn't a valid VCF line
				write_to_log("Error parsing: ", vcf_log_location)
				line = f.readline()
				continue
			if spline[4] == '.':
				# It's a "no variant". Ignore it.
				# This hasn't been tested too rigorously because it didn't appear in any of my test data,
				# but someone else had some of these so I had him insert this and it seemed to work for him.
				#write_to_log("No variant: "+line, vcf_log_location) # Uncomment this if you want to write these lines to the log for testing
				line = f.readline()
				continue

			# Not sure why we subtract one from pos, but it happened in the original.
			# Also, if there's no QUAL, will it always be a dot? If not, we need to figure out how to deal
			# If there's no QUAL and there's a dot, make it just over the threshold.
			try:
				if spline[5] == '.':
					spline[5] = quality_threshold+.01
				chr, pos, id, old, new_array, qual = spline[0].lstrip('chr'), int(spline[1])-1, spline[2], spline[3], spline[4].split(','), float(spline[5])
			except ValueError:
				write_to_log("Error parsing: "+line.rstrip(), vcf_log_location)
				line = f.readline()
				continue				
					
			# If we passed those checks, let's say it's good...
			total_variants += len(new_array) # some lines have multiple variants, annoyingly
					
			# Quality check!
			if qual < quality_threshold:
				qual_removed += len(new_array)
				line = f.readline()
				continue

			for new in new_array:
				# Is there a higher-quality version of this variant already found?
				line_map_key = "%s#%d#%s#%s" % (chr, pos, old, new)
				if line_map_key in existing_variants:
					if existing_variants[line_map_key][5] > qual:
						duplicates_removed += 1
						continue
					else:
						old_duplicates_removed += 1

				# Passed the other checks? Cool, store the first six fields in our map and continue.
				# Overwrite the duplicate if it exists.
				kept_variants += 1
				existing_variants[line_map_key] = [chr, pos, id, old, new, qual]

			line = f.readline()

		# Add everything up and write a status line to the log
		#header_written = True # Don't want to write a header unless it's the first file
		qual_removed_all += qual_removed
		total_variants_all += total_variants
		kept_variants_all += kept_variants
		duplicates_all += (duplicates_removed+old_duplicates_removed)
		write_to_log("\nFrom variant file %s: %d total variants, %d failed quality check at threshold %f, %d duplicates found \(%d overwritten\), %d variants kept in final version" % (vf, total_variants, qual_removed, quality_threshold, (duplicates_removed+old_duplicates_removed), old_duplicates_removed, kept_variants), vcf_log_location)
		write_to_log("------------------------", vcf_log_location)
		f.close()
		
	# Write all variants to merged.vcf
	for key in sorted(existing_variants.keys()):
		w.write('chr'+'\t'.join(str(i) for i in existing_variants[key])+'\n')

	# Write final status line to the log and close all files, we're done.
	write_to_log("\nTotal for all variant files: %d total variants, %d failed quality check at threshold %f, %d duplicates removed, %d variants kept in final version" % (total_variants_all, qual_removed_all, quality_threshold, duplicates_all, len(existing_variants.keys())), vcf_log_location)
	w.close()
	vcf_log.close()
	return None

### This function removes any duplicate germline/somatic variants from the somatic file.

def remove_somatic_duplicates(germ_dir, som_dir):
	'''Removes duplicate germline/somatic variants from the somatic file.'''
	# Shouldn't have to check and make sure the files are there - that should have happened earlier
	# Open the germ file, save all of the chr/pos/old/new in a hash table for easy lookup
	# Then open the somatic file, go through each line, if one matches just don't write it back out
	# Afterwards, overwrite the old somatic merged.vcf with the new one
	
	global logfile
	
	# Open the germline file, save all of the information in a hash table for easy lookup, close it
	f = open(germ_dir+'/merged_pytest/merged.vcf','r')
	germ_variants = set([])
	line = f.readline()
	while line:
		if line[0] == '#':
			# found a header
			line = f.readline()
			continue
		spline = line.split('\t')
		variant = "%s#%s#%s#%s" % (spline[0], spline[1], spline[3], spline[4]) # a.k.a. chr, pos, old, new
		germ_variants.add(variant)
		line = f.readline()
	f.close()
	
	w = open(som_dir+'/merged_pytest/temp.vcf','w')
	f = open(som_dir+'/merged_pytest/merged.vcf','r')
	duplicates_found = 0
	kept_variants = 0
	line = f.readline()
	while line:
		if line[0] == '#':
			# found a header
			w.write(line)
			line = f.readline()
			continue
		spline = line.split('\t')
		variant = "%s#%s#%s#%s" % (spline[0], spline[1], spline[3], spline[4]) # a.k.a. chr, pos, old, new
		if variant in germ_variants:
			duplicates_found += 1
		else:
			kept_variants += 1
			w.write(line)
		line = f.readline()
		
	w.close()
	f.close()
	write_to_log("Found %d duplicates between somatic and germline variant files. Removed from somatic file. %d somatic variants remain." % (duplicates_found, kept_variants),logfile)
	shutil.move(som_dir+'/merged_pytest/temp.vcf', som_dir+'/merged_pytest/merged.vcf')

### This function is used to create .bed files of only variants that exist in exons.

def get_variants(vcf_file, proteome_file, type):
	f = open(proteome_file, 'r')
	
	# Basically goes through proteome.bed and finds all positions in exons
	# Then goes through the .vcf file and checks if they're in any exons
	# Output line: NP_XX (name) \t type-old(position in genome)new:(quality),(repeat)\t type-old(position in exome)new:(quality),(repeat) Looks like instead of quality it might have taken the wrong field at some point? Sometimes it looks like instead of quality we just get a dot, even when it shouldn't be missing. I'll leave quality in for now and see if we use it later. Output all genes, even if they have no variants.
	
	# QUESTION: In the original, it appears that in proteome.bed.var the pos-- doesn't happen? Why?
	# Also, somatic hasn't been filtered yet at this step in the original. I think it should be okay that that's changed here - 
	# the earlier you filter it, the faster the other steps are.
	
	est = ExonSearchTree()
	line = f.readline() # Going to assume no headers for now
	all_names = set([])
	while line:
		spline = line.rstrip().split('\t')
		try:
			chr, start, name, lengths, offsets = spline[0].lstrip('chr'), int(spline[1]), spline[3], spline[-2].rstrip(','), spline[-1].rstrip(',')
		except ValueError:
			# Maybe write this to a log somewhere?
			warnings.warn("Failed to parse %s" % line)
			line = f.readline()
			continue

		all_names.add(name)
		
		splengths = lengths.split(',')
		spoffsets = offsets.split(',')
		total_exon_length = 0
		for i in range(len(spoffsets)):
			# The original saved a map with every position that ever appears in an exon
			# and the names of the genes it goes with. This feels slow and inefficient, so instead
			# I made some wacked-out tree thing for storage and search. See if you like it!
			# Whoa okay, tested it, it's much faster.
			est.add_exon(chr, int(spoffsets[i])+start, int(spoffsets[i])+start+int(splengths[i])-1, total_exon_length, name)	
			total_exon_length += int(splengths[i])		
		line = f.readline()
	#print est.total_exons
	f.close()

	f = open(vcf_file, 'r')
	line = f.readline() # Still assuming no headers
	all_genes = {}
	while line:
		spline = line.rstrip().split('\t')
		try:
			chr, pos, id, old, new, qual = spline[0].lstrip('chr'), int(spline[1]), spline[2], spline[3], spline[4], float(spline[5])
		except ValueError:
			# Maybe write this to a log somewhere?
			warnings.warn("Failed to parse %s" % line)
			line = f.readline()
			continue
		exon = est.find_exon(chr,pos)
		if exon != []:
			# Save pos in chr and pos in gene
			# Each exon returned is a [name, position in gene] pair
			for ex in exon:
				in_chr = "%s-%s%d%s:%f" % (type, old, pos, new, qual) 
				in_gene = "%s-%s%d%s:%f" % (type, old, ex[1], new, qual) 
				try:
					all_genes[ex[0]].append([in_chr, in_gene])
				except KeyError:
					all_genes[ex[0]] = []
					all_genes[ex[0]].append([in_chr, in_gene])
		line = f.readline()
	f.close()

	# Write variants out to file.
	# Need to let it write genes with no variants, also.
	w = open(proteome_file+"."+type+".var", 'w')	
	for key in all_names:
		in_chr = []
		in_gene = []
		try:
			vars = all_genes[key]
		except KeyError:
			vars = []
		for var in vars:
			in_chr.append(var[0])
			in_gene.append(var[1])
		w.write("%s\t%s\t%s\n" % (key, ','.join(in_chr), ','.join(in_gene)))
	w.close()

### These functions are used to sort variants by type.

def process_gene(header_line, second_header, exon_headers, exon_seqs, variants):
	'''Checks all variants in a gene for whether or not they cause a single-AA non-stop substitution. Returns a list of the ones that do, and the AA substitutions they cause.'''
	global codon_map
	global logfile
	
	# Ready the translation table for the reverse strand.
	translate_table = maketrans("ACGTacgt","TGCAtgca")
	
	# For each variant in our variants, check if it will become an AA substitution.
	# Do I check only a single reading frame? That's what they did. 
	# Also: if the length of old != length of new, it's an indel, so ignore for now.	
	full_seq = ''.join(exon_seqs[1:-1])
	
	# Set a flag for a reverse gene
	reverse_flag = False
	if header_line.split('\t')[5] == '-':
		reverse_flag = True
	
	changes = []
	aa_substs = []
	indels = []
	indel_substs = []

	# With hg38, we get a weird thing where some of the genes in proteome.bed.dna have no sequence associated with them.
	# So...I guess I'll just skip those for now until I figure out why
	# I think this is fixed, and it was a problem in a different script. Keeping this in anyway, just in case.
	if '0' in full_seq  or 'N' in full_seq:
		write_to_log(header_line+full_seq, logfile)
		return changes, aa_substs, indels
	
	prev_start = -1 # Previous codon start position
	prev_subst = [-1, ''] # Previous substitution position and nucleotide
	for var in variants:
		# Remove the ones with multiple nt in the variant
		# I'm going to remove the ones where both sides have the same number but that number is greater than 1
		# Just for now.
		pos = int(re.findall(r'\d+', var)[0]) # using the -1 so we're counting from 0, not 1
		orig_nt = var.split('-')[1].split(str(pos))[0]
		new_nt = var.split(':')[0].split(str(pos))[1]
		
		# Here we removed the ones where it wasn't a 1 to 1 substitution.
		# Now we automatically add it to the changes.
		if len(orig_nt) != 1 or len(new_nt) != 1:
			if reverse_flag:
				#pos = len(full_seq) - pos
				#orig_nt = orig_nt[::-1].translate(translate_table)
				#new_nt = new_nt[::-1].translate(translate_table)
				qual = var.split(':')[1]
				prefix = var.split('-')[0]
				indels.append('%s-%s%d%s:%s' % (prefix, orig_nt, pos, new_nt, qual))
			else:
				indels.append(var)
			continue
		
		triplet_start = ((pos/3)*3) # start of the triplet containing pos. counting from 0, not 1
		triplet_orig = full_seq[triplet_start:(triplet_start+3)].upper()
		triplet_subst = var.split(':')[0].split(str(pos))[-1]
		if triplet_orig == '':
			print '\n'+'Empty codon: '+var, full_seq, len(full_seq), triplet_start, header_line
			continue
		if len(triplet_orig) != 3:
			print '\n'+'Triplet with length <3 (length of sequence likely not divisible by 3, for whatever reason): '+var, full_seq, len(full_seq), triplet_start, header_line
			continue
		subst_pos = pos%3
		triplet_new = triplet_orig[:subst_pos] + triplet_subst + triplet_orig[subst_pos+1:] # This did! change it in both.
		# If it's the reverse strand, check the reversed and translated triplet instead
		# I guess the .vcf files are based on proteome.bed.dna? So the position and nucleotides refer
		# to the negative strand if proteome.bed.dna only has the negative strand.
		if reverse_flag:
			triplet_orig = triplet_orig[::-1].translate(translate_table)
			triplet_new = triplet_new[::-1].translate(translate_table)
		AA_old = codon_map[triplet_orig]
		AA_new = codon_map[triplet_new]
		# If this was a start codon, ignore it for the time being.
		if AA_old == 'M' and triplet_start == 0:
			#print "Variant in start codon"
			#print header_line.split('\t')[3]+'\t'+var+'\t'+triplet_orig+'\t'+triplet_new+'\t'+full_seq[subst_pos-3:subst_pos+6]
			write_to_log("Start codon variant found: %s" % (header_line.split('\t')[3]+'\t'+var+'\t'+triplet_orig+'\t'+triplet_new+'\t'+full_seq[subst_pos-3:subst_pos+6]),logfile)
			continue
		# If neither the old nor the new codon is a stop codon and the AA changes,
		# count it!
		#if AA_old != AA_new and AA_new != '*' and AA_old != '*':
		# Now we're just checking whether the old and new AAs are different. Leaving stop codon things in...
		if AA_old != AA_new:
			changes.append(var)
			if reverse_flag:
				total_AA = len(full_seq)/3
				aa_substs.append((AA_old, total_AA-(pos/3), AA_new))
			else:
				aa_substs.append((AA_old, pos/3+1, AA_new))
		# If the new codon is a stop codon, 
		# uh...I found this comment here and it looks like I decided not to do whatever I was planning to do.
		
		# Check to see if this is the second substitution in a codon - if so, does the amino acid change with both?
		if prev_start == triplet_start:
			triplet_orig = full_seq[triplet_start:(triplet_start+3)].upper()
			triplet_new = list(triplet_orig)
			triplet_new[subst_pos] = triplet_subst
			triplet_new[prev_subst[0]] = prev_subst[1]
			triplet_new = ''.join(triplet_new)
			if reverse_flag:
				triplet_orig = triplet_orig[::-1].translate(translate_table)
				triplet_new = triplet_new[::-1].translate(translate_table)
			AA_old = codon_map[triplet_orig]
			AA_new = codon_map[triplet_new]
			#print header_line.split('\t')[3]+'\t'+var+'\t'+str(prev_subst[0])+'\t'+prev_subst[1]+'\t'+triplet_orig+'\t'+triplet_new+'\t'+full_seq[subst_pos-3:subst_pos+6]
			#print full_seq
			#if AA_old != AA_new and AA_new != '*' and AA_old != '*':
			# Now we're just checking whether the old and new AAs are different. Leaving stop codon things in...
			if AA_old != AA_new:
				# G-A2158G:2281.770000
				change = "X-%s%d%s:0.0" % (triplet_orig, triplet_start, triplet_new)
				changes.append(change)
				if reverse_flag:
					total_AA = len(full_seq)/3
					aa_substs.append((AA_old, total_AA-(pos/3), AA_new))
				else:
					aa_substs.append((AA_old, pos/3+1, AA_new))
		prev_start = triplet_start
		prev_subst = [subst_pos, triplet_subst]
	return changes, aa_substs, indels	
	
def write_out_aa(name, changed_vars, aa_substs, out_aa):
	'''Creates the proteome.bed.aa.var file, which is just a list of all variants in the style of proteome.bed.var.'''
	aa_sub_to_write = []
	for aa_sub in aa_substs:
		aa_sub_to_write.append(aa_sub[0]+str(aa_sub[1])+aa_sub[2])
	out_aa.write("%s\t%s\t%s\n" % (name, ','.join(changed_vars), ','.join(aa_sub_to_write)))

def write_out_aa_bed(header_line, changed_vars, out_aa_bed):
	'''Creates the proteome.aa.var.bed file, which is basically a collection of the header lines for each entry in proteome.aa.var.bed.dna.'''
	to_change = header_line[3]
	for var in changed_vars:
		v = var.split(':')[0]
		header_line[3] = to_change+'-'+v
		out_aa_bed.write('\t'.join(header_line))

def write_out_aa_bed_dna(header_line, second_header, exon_headers, changed_vars, aa_substs, out_aa_bed_dna):
	'''Creates the proteome.aa.bed.var.dna file, which is the same as proteome.bed.dna but instead of an entry for each gene, has an entry for each variant of each gene. It's a bit unwieldy, but it is what it is.'''
	# Maybe also add a line with all of the variants?
	header_to_change = header_line[3]
	second_header_name = second_header[0]
	second_header_map_end = second_header[-1]
	for i in range(len(changed_vars)):
		var = changed_vars[i]
		subst = aa_substs[i]
		# First line:
		v = var.split(':')[0]
		header_line[3] = header_to_change+'-'+v
		out_aa_bed_dna.write('\t'.join(header_line))
		# Second line:
		second_header[0] = second_header_name+'-'+v
		second_header[-1] = second_header_map_end[:-1]+' '+var+')'
		out_aa_bed_dna.write(' '.join(second_header)+(' (VAR:%s-%s%d%s:%s)\n' % (var[0], subst[0], subst[1], subst[2], var.split(':')[1]))) # Will write the AA position starting from 0, not 1. 
		# Sequence lines:
		out_aa_bed_dna.write('\t'.join(exon_headers[0]))
		var_pos = int(re.findall(r'\d+', var)[0])
		seq_length = 0
		edit_made = False
		for i in range(1,len(exon_headers)-1):
			chunk_length = int(exon_headers[i][0].split('\t')[4])
			if edit_made or seq_length+chunk_length <= var_pos:
				out_aa_bed_dna.write('\t'.join(exon_headers[i]))
			else:
				seq_chunk = exon_headers[i][-1]
				seq_var = var.split(':')[0].split(str(var_pos))[1]
				chunk_pos = var_pos-seq_length
				seq_chunk = seq_chunk[:chunk_pos]+seq_var+seq_chunk[chunk_pos+len(seq_var):]
				out_aa_bed_dna.write(exon_headers[i][0]+'\t'+seq_chunk)
				edit_made = True
			seq_length += chunk_length
		out_aa_bed_dna.write('\t'.join(exon_headers[-1]))

def test_indels(exon_headers, indels):
	'''Goes through all of the indels to make sure they are actually deleting parts of the exon.
	Some of them don't - their first character, which gets added back in (G-AC576A, S-ACGGT239A, etc) 
	is part of an exon, and the rest are not. Also removes any indel that's part of a start codon, because that causes a messy cascade of changes.'''
	true_indels = []
	for indel in indels:
		indel_pos = int(re.findall(r'\d+', indel)[0])
		if indel_pos < 3: # it's in a start codon
			continue
		seq_length = 0
		for i in range(1, len(exon_headers)-1):
			ex_head = exon_headers[i][0].split('\t')
			chunk_length = int(ex_head[4])
			if seq_length+chunk_length > indel_pos:
				chunk_pos = indel_pos - seq_length
				if chunk_pos != chunk_length-1:
					true_indels.append(indel)
				#else:
				#	write_to_status("Removing %s %s" % (ex_head[0], indel))
				break
			seq_length = seq_length+chunk_length
	return true_indels
		
def write_out_indel_bed_dna(header_line, second_header, exon_headers, indels, out_indel_bed_dna):
	'''Creates the proteome.indel.bed.var.dna file, which is the same as proteome.bed.dna but instead of an entry for each gene, has an entry for each indel.'''
	header_to_change = header_line[3]
	second_header_name = second_header[0]
	second_header_map_end = second_header[-1]
	for indel in indels:
		# First line:
		header_line[3] = header_to_change+'-indel'
		out_indel_bed_dna.write('\t'.join(header_line))
		# Second line:
		second_header[0] = second_header_name+'-indel'
		second_header[-1] = second_header_map_end[:-1]+' '+indel+')'
		out_indel_bed_dna.write(' '.join(second_header)+'\n')
		# Sequence lines:
		out_indel_bed_dna.write('\t'.join(exon_headers[0])) # Pre-sequence
		indel_pos = int(re.findall(r'\d+', indel)[0])
		seq_length = 0
		to_delete = indel.split(':')[0].split(str(indel_pos))[0].split('-')[1]
		to_insert = indel.split(':')[0].split(str(indel_pos))[1]
		edit_made = False
		for i in range(1,len(exon_headers)-1):
			ex_head = exon_headers[i][0].split('\t')
			chunk_length = int(ex_head[4])
			if edit_made or seq_length+chunk_length <= indel_pos:
				out_indel_bed_dna.write('\t'.join(exon_headers[i]))
			else:
				orig_seq_chunk = exon_headers[i][-1]
				chunk_pos = indel_pos-seq_length
				# What to do in this eventuality? 
				# Because the deletion extends into the intron and not the next exon, just ignore it and make sure you edit the length appropriately.
				seq_chunk = orig_seq_chunk[:chunk_pos]+to_insert+orig_seq_chunk[chunk_pos+len(to_delete):]
				if chunk_pos+len(to_delete) <= chunk_length:
					new_length = chunk_length-len(to_delete)+len(to_insert)
				else:
					new_length = chunk_pos+1 # adds one for the base we're adding back in
					try:
						write_to_status("Found a deletion that goes beyond the exon. %s\t%s" % (ex_head[0], indel))
					except OSError:
						print "Found a very strange deletion: " + '\t'.join(ex_head) + '\t' + indel
				tmp_ex_head = ex_head
				tmp_ex_head[4] = str(new_length)
				ex_head_to_write = '\t'.join(tmp_ex_head)
				out_indel_bed_dna.write(ex_head_to_write+'\t'+seq_chunk.rstrip()+'\n')
				edit_made = True
				# Error testing, remove later
				# Which ones are wrong and why? Which ones are right and why? Print everything
				# Also, print what the variant is.
				if orig_seq_chunk[chunk_pos:chunk_pos+len(to_delete)].upper() != to_delete.upper() and chunk_pos+len(to_delete) <= chunk_length:
					write_to_status("Wrong sequence being deleted! %s\t%s\t%s\t%s" % (ex_head[0], indel, seq_chunk[chunk_pos:chunk_pos+len(to_delete)].upper(), orig_seq_chunk))
			seq_length += chunk_length
		out_indel_bed_dna.write('\t'.join(exon_headers[-1])) # Post-sequence	

def sort_variants(proteome_file, variant_file, output_prefix):
	'''Goes through the variants and sorts them by type, writing out a bunch of intermediate files in the process. Right now the types are "single-AA non-stop variant" and "not that", but eventually will have more.'''
	global codon_map
	
	# Grab all the variants and their locations within their genes
	variants_map = {}
	f = open(variant_file, 'r')
	line = f.readline()
	while line:
		spline = line.rstrip().split('\t')
		if len(spline) > 1:
			split_variants = spline[2].split(',')
			variants_map[spline[0]] = variants_map.get(spline[0], []) + split_variants
		line = f.readline()
	f.close()
	
	# Write those to a file, just for kicks/error checking (might remove this later)
	file_base = proteome_file.rsplit('/',1)[0]
	w = open(variant_file+".SG.combined", 'w')
	for var in variants_map:
		w.write("%s\t%s\n" % (var, ','.join(variants_map[var])))
	w.close()
	
	# Open some files to write to. Right now, only looking at single AA changes.
	# I hate these filenames. Do they actually mean anything?
	out_aa = open(file_base+"/"+output_prefix+".bed.aa.var", 'w')
	# This one will basically just be that first line. Might want an extra line for each variant though, depends how it's used. Actually, almost certainly do.
	out_aa_bed = open(file_base+"/"+output_prefix+".aa.var.bed", 'w')
	# proteome.bed.dna, but with the variants included. Def want an extra line for each variant.
	out_aa_bed_dna = open(file_base+"/"+output_prefix+".aa.var.bed.dna", 'w')
	# Right now just proteome.bed.aa.var but with the non-AA variants. Doesn't do anything right now because I don't care about it very much.
	out_indel = open(file_base+"/"+output_prefix+".bed.indel.var", 'w')
	out_indel_bed = open(file_base+"/"+output_prefix+".indel.var.bed", 'w')
	out_indel_bed_dna = open(file_base+"/"+output_prefix+".indel.var.bed.dna", 'w')
	out_other = open(file_base+"/"+output_prefix+".bed.other.var", 'w')
	
	f = open(proteome_file, 'r')
	line = f.readline()
	while line:
		# Line type one: First gene header line
		if line[:3] == 'chr':
			# Finish processing the previous gene - focus first on proteome.bed.aa.var
			# Using a try/except block here to deal with the possibility of this being the first line
			# There has to be a more graceful way to do that, right?
			try:
				changed_vars, aa_substs, indels = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []))
				if changed_vars != []:
					write_out_aa(name, changed_vars, aa_substs, out_aa)
					write_out_aa_bed(header_line.split('\t'), changed_vars, out_aa_bed)
					write_out_aa_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, changed_vars, aa_substs, out_aa_bed_dna)
				if indels != []:
					# Make sure none of them are just removing intron stuff
					indels = test_indels(exon_headers, indels)
					# If there are any left...
					if indels != []:
						out_indel.write("%s\t%s\n" % (name, ','.join(indels)))
						out_indel_bed.write("chr%s\t%d\t%d\t%s-indel\t1000\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (chr, start, end, name, strand, start, end, spline[-4], exon_count, ','.join(exon_lengths), ','.join(exon_offsets)))
						write_out_indel_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, indels, out_indel_bed_dna)
			except UnboundLocalError:
				pass
			# Start processing the new gene
			exon_headers = []
			exon_seqs = []
			header_line = line # To be printed to proteome.aa.var.bed, in some form
			spline = line.rstrip().split('\t')
			chr, start, end, name, strand, exon_count, exon_lengths, exon_offsets = spline[0].lstrip('chr'), int(spline[1]), int(spline[2]), spline[3], spline[5], int(spline[-3]), spline[-2].split(','), spline[-1].split(',')
		# Line type two: Second gene header line
		elif line[0] == '>':
			second_header = line
		# Line type three: Exon line (can be pre-100, post-100 or internal)
		else:
			exon_headers.append(line.rsplit('\t',1))
			exon_seqs.append(line.rstrip().split('\t')[-1])
		line = f.readline()
		
	# Do the last one
	try:
		changed_vars, aa_substs, indels = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []))
		if changed_vars != []:
			write_out_aa(name, changed_vars, aa_substs, out_aa)
			write_out_aa_bed(header_line.split('\t'), changed_vars, out_aa_bed)
			write_out_aa_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, changed_vars, aa_substs, out_aa_bed_dna)
		if indels != []:
			# Make sure none of them are just removing intron stuff
			indels = test_indels(exon_headers, indels)
			# If there are any left...
			if indels != []:
				out_indel.write("%s\t%s\n" % (name, ','.join(indels)))
				out_indel_bed.write("chr%s\t%d\t%d\t%s-indel\t1000\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (chr, start, end, name, strand, start, end, spline[-4], exon_count, ','.join(exon_lengths), ','.join(exon_offsets)))
				write_out_indel_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, indels, out_indel_bed_dna)
	except UnboundLocalError:
		pass

	# Close everything
	f.close()
	out_aa.close()
	out_aa_bed.close()
	out_aa_bed_dna.close()
	out_indel.close()
	out_indel_bed.close()
	out_indel_bed_dna.close()
	out_other.close()
	
### These functions are used to translate DNA variants into variant proteins.

def translate_seq(sequence, strand, return_all = False):
	'''Translates a DNA sequence to an AA sequence'''
	global codon_map
	
	if strand == '-':
		translate_table = maketrans("ACGTacgt","TGCAtgca")
		sequence = sequence[::-1].translate(translate_table)
	
	# Didn't find any of these yet, just takes up time. Should probably have a check like this in there more formally.
	#if len(sequence)%3 != 0:
	#	warnings.warn("Found sequence with length indivisible by 3!")
	
	translated = ''
	
	'''
	for i in range(0,len(sequence)-5,3):
		# Translates the whole thing and removes the last character (a stop codon)
		# Obsolete, since now we can have stop codons ANYWHERE
		translated += codon_map[sequence[i:i+3]]
	return translated
	'''
	
	for i in range(0,len(sequence)-2,3):
		translated += codon_map.get(sequence[i:i+3],'X')
	
	# Now this should return the translated sequence up to the first stop codon.	
	# That could either be at the actual end of the sequence, or at a variant that adds a stop.
	if return_all:
		# In some cases, I want to return the whole sequence, including internal stops.
		return translated
	else:
		return translated.split('*')[0]

def calculate_chr_pos(map_section):
	'''Finding the chromosomal position of a variant'''
	gene_start, strand = int(map_section.split()[0][:-1]),map_section.split()[0][-1]
	lengths = map(int,map_section.split()[1].rstrip(',').split(','))
	starts = map(int,map_section.split()[2].rstrip(',').split(','))
	snp = map_section.split()[-1].split('-')[-1]
	snp_pos = int(re.findall(r'\d+', snp)[0])
	orig_nt, new_nt = snp.split(str(snp_pos))[0],snp.split(str(snp_pos))[-1]
	tot_len = 0
	for i in range(len(lengths)):
		if tot_len+lengths[i] > snp_pos:
			start_sec = gene_start+starts[i]
			start_sec += snp_pos-tot_len
			return "%s%d%s" % (orig_nt,start_sec,new_nt)
		else:
			tot_len += lengths[i]

def format_aa_header(orig_header,abbr,desc):
	#>NP_001909-G-T299C (MAP:chr1:100661810- 168,72,192,78,167,217,122,182,76,124,51 0,9975,10190,14439,18562,19728,22371,34478,39181,44506,53515 G-T299C:10378.200000) (VAR:G-S384G:10378.200000)
	#to: >NP_XXXX-V500G|proteinDescription|GN=GeneID|chr=N|type=G|SNP=G12345C|QUAL=8|SAAP=V500G
	id = orig_header.split('-')[0]
	chr = orig_header.split(':')[1].lstrip('chr')
	type = orig_header.split('-')[1]
	#SNP_orig = orig_header.split('-')[2].split()[0]
	SNP = calculate_chr_pos(orig_header.split(':')[2])
	qual = orig_header.split(':')[-1].split(')')[0]
	SAAP = orig_header.split('-')[-1].split(':')[0]
	return "%s-%s|%s|GN=%s|chr=%s|type=%s|SNP=%s|qual=%s|SAAP=%s\n" % (id,SAAP,desc,abbr,chr,type,SNP,qual,SAAP)

def format_indel_header(orig_header,abbr,desc):
	#>NP_005057-indel (MAP:chr1:35650056- 138,122,49,118,85,197,96,302,189,828 0,2545,2747,3517,4547,4730,6042,6238,6885,7766 G-ACCA2082CCCC:0.000000)
	#to: >NP_005057-ACCA2082CCCC|proteinDescription|GN=GeneID|chr=N|type=G|QUAL=?
	id = orig_header.split('-')[0]
	chr = orig_header.split(':')[1].lstrip('chr')
	type = orig_header.split('-')[-2][-1]
	#indel_orig = orig_header.split('-')[-1].split(':')[0]
	indel = calculate_chr_pos(orig_header.split(':')[2])
	qual = orig_header.split(':')[-1].split(')')[0]
	return "%s-%s|%s|GN=%s|chr=%s|type=%s|indel=%s|qual=%s\n" % (id,indel,desc,abbr,chr,type,indel,qual)
	
def format_juncA_header(orig_header,abbr,desc):
	#>NP_064587-A-5-7-3-chr3-100064522-100067646-75 (MAP:chr3:100053635+ 7,119,121,89,94,79,99,56,92 0,4295,5023,6281,10793,14011,17612,19980,20385)
	# to: >NP_005057-bothJunc-100464971-100476920|proteinDescription|GN=GeneID|chr=N|type=bothConserved|donorExon=2|donorSite=100464971|acceptorExon=4|acceptorSite=100476920|reads=8	
	# New header: >NP_003318-3-5-6-1147518-1148372-JUNC_369 (MAP:chr1:1146934- 71,129,197,102,123,145 0,149,387,1437,2108,2428)
	id = orig_header.split('-')[0]
	donor = orig_header.split('-')[1]
	acceptor = orig_header.split('-')[2]
	reads = orig_header.split('-')[3]
	start = orig_header.split('-')[4]
	end = orig_header.split('-')[5]
	chr = orig_header.split('chr')[1].split(':')[0]
	strand = orig_header.split()[1][-1]
	
	return "%s-bothJunc-chr%s-%s-%s|%s|GN=%s|chr=%s|strand=%s|type=bothConserved|donorExon=%s|donorSite=%s|acceptorExon=%s|acceptorSite=%s|reads=%s\n" % (id, chr, start, end, desc, abbr, chr, strand, donor, start, acceptor, end, reads)
		
def format_juncAN_header(orig_header,abbr,desc):
	#>NP_060779-AN1-10-2-chr3-100018175-100020919 (MAP:chr3:99979862+ 53,112,106,205,124,125,47,104,123,93,163,43,107,140,89,136,195,82 0,18630,20723,22588,29559,34068,34283,35153,36904,38220,41057,43827,45418,49384,55080,58050,59758,62568)
	#>NP_055635-AN2-9-1-chr3-100100578-100103322 (MAP:chr3:100084407- 154,123,98,117,108,135,208,120,110,127,174,324 0,2480,3474,7042,7974,9454,12141,16051,18915,20654,21241,35062)
	#to: >NP_005057-donorJunc-100464971-100476920|proteinDescription|GN=GeneID|chr=N|strand=-|type=donorConserved|donorExon=2|donorSite=100464971|acceptorSite=100476920|reads=8
	# New: >NP_001034300-4-4-1389880-1421162-JUNC_635 (MAP:chr1:1386063+ 75,77,70,156,46 0,1362,1681,3661,35099)
		
	id = orig_header.split('-')[0]
	conserved_exon = orig_header.split('-')[1]
	reads = orig_header.split('-')[2]
	start = orig_header.split('-')[3]
	end = orig_header.split('-')[4].split()[0]
	chr = orig_header.split('chr')[1].split(':')[0]
	strand = orig_header.split()[1][-1]
	
	if strand == '+':
		return "%s-donorJunc-chr%s-%s-%s|%s|GN=%s|chr=%s|strand=%s|type=donorConserved|donorExon=%s|donorSite=%s|acceptorSite=%s|reads=%s\n" % (id, chr, start, end, desc, abbr, chr, strand, conserved_exon, start, end, reads)
	else:
		return "%s-donorJunc-chr%s-%s-%s|%s|GN=%s|chr=%s|strand=%s|type=donorConserved|donorExon=%s|donorSite=%s|acceptorSite=%s|reads=%s\n" % (id, chr, start, end, desc, abbr, chr, strand, conserved_exon, end, start, reads)

def format_juncN_header(orig_header):
	#>NO_GENE-N-1-chr3-j268-101520340-101520702 (MAP:chr3:101520289+ 50,50 0,414)-2-
	#>NP_055230-AN2-9-1-chr3-101384962-101389973 (MAP:chr3:101384911+ 50,50 0,5062)-2-
	#>NP_055575-AN1-5-12-chr3-10320146-10320347 (MAP:chr3:10320096+ 50,50 0,251)-0-
	# new: >NO_GENE-4-chr1-18369-18913-JUNC_18 (MAP:chr1:18267+ 102,82 0,645)
	#to: >NovelJunc-100464971-100476920|chr=N|type=novel|site1=100464971|site2=100476920|reads=8|Frame=1+
	start = orig_header.split('-')[3]
	end = orig_header.split('-')[4]
	chr = orig_header.split('chr')[1].split('-')[0]
	reads = orig_header.split('-')[1]
	frame = orig_header.rstrip().split(')-')[-1]
	return ">NovelJunc-chr%s-%s-%s|chr=%s|type=novel|site1=%s|site2=%s|reads=%s|frame=%s\n"% (chr,start,end,chr,start,end,reads,frame)

def format_header(orig_header,seq_type,abbr,desc):
	if seq_type=='aa':
		return format_aa_header(orig_header,abbr,desc)
	elif seq_type=='indel':
		return format_indel_header(orig_header,abbr,desc)
	elif seq_type=='juncA':
		return format_juncA_header(orig_header,abbr,desc)
	elif seq_type=='juncAN':
		return format_juncAN_header(orig_header,abbr,desc)
	elif seq_type=='juncN':
		return format_juncN_header(orig_header)
	return orig_header

def translate(log_dir, bed_file, logfile, seq_type):
	'''Reads in a bed.dna file and translates it to a protein .fasta file.'''
	# Get the descriptions
	desc = {}
	f = open(log_dir+"proteome-descriptions.txt",'r')
	line = f.readline()
	while line:
		spline = line.rstrip().split('\t')
		if len(spline) > 1:
			desc[spline[0]] = spline[1]
		else:
			desc[spline[0]] = ""
		line = f.readline()
	f.close()
	# Get the gene abbreviations
	abbr = {}
	f = open(log_dir+"proteome-genes.txt",'r')
	line = f.readline()
	while line:
		spline = line.rstrip().split()
		abbr[spline[0]] = spline[-1]
		line = f.readline()
	f.close()
	
	# The rest of it
	f = open(log_dir+bed_file,'r')
	out_fasta = open(log_dir+bed_file+".fasta",'w') # The original QUILTS writes two other files but they're just duplicates of proteome.aa.var.bed and proteome.aa.var.bed.dna, it seems. For now I'm leaving them out.
	line = f.readline()
	sequence = ""
	prev_AA_subst = []
	prev_gene = ''
	while line:
		# Line type one: First gene header line
		if line[:3] == 'chr':
			# Finish processing the previous gene - focus first on proteome.bed.aa.var
			# Using a try/except block here to deal with the possibility of this being the first line
			# There has to be a more graceful way to do that, right?
			if "Error Reading Sequence" in sequence: # This error comes from read_chr_bed, which I did not write and am afraid to touch.
				write_to_log("Error in sequence for %s. Skipping..." % gene, logfile)
			else:
				try:
					if prev_gene != gene:
						prev_AA_subst = []
					if AA_subst not in prev_AA_subst:
						translated = translate_seq(sequence.upper(), strand)
						out_fasta.write(format_header(second_header, seq_type, abbr[gene], desc[gene]))
						#out_fasta.write(second_header)
						out_fasta.write(translated+'\n')
					prev_AA_subst.append(AA_subst)
					prev_gene = gene
				except UnboundLocalError:
					pass
			# Start processing the new gene
			sequence = ""
			strand = line.split('\t')[5]
			gene = line.split('\t')[3].split('-')[0]
		# Line type two: Second gene header line
		elif line[0] == '>':
			second_header = line # This will be the header in the fasta file
			AA_subst = line.rsplit('-',1)[-1].split(':')[0]
		# Line type three: Exon line (can be pre-100, post-100 or internal)
		else:
			spline = line.rstrip().split('\t')
			# Basically, if the strand is reversed, we want the -1 part because that'll give us the extra stuff
			# that gets added in if the stop codon is deleted. But if the strand is forward, we keep the +1 instead.
			if spline[1] != '+1' and spline[1] != '-1':
				sequence += spline[-1]
			elif spline[1] == '-1' and spline[5] == '-':
				sequence += spline[-1]
			elif spline[1] == '+1' and spline[5] =='+':
				sequence += spline[-1]
		line = f.readline()
		
	# And finish the last one
	if "Error Reading Sequence" in sequence:
		write_to_log("Error in sequence for %s. Skipping..." % logfile)
	else:
		try:
			if prev_gene != gene or AA_subst != prev_AA_subst:
				translated = translate_seq(sequence.upper(), strand)
				out_fasta.write(format_header(second_header, seq_type, abbr[gene], desc[gene]))
				out_fasta.write(translated+'\n')
		except UnboundLocalError:
			pass

	f.close()
	out_fasta.close()

### These functions are used to make the peptide fasta files.

def trypsinize(sequence, pos = 0, end_pos = False, allow_missed_cleavage = False):
	'''Virtually trypsinizes a protein and returns its tryptic peptides. 
	If a position is given, it only returns the peptides following that AA position (used for indels
	and splice junctions, where anything before the variant is unchanged). Also can specify an end_pos.
	Right now, just chops it after any K or R that isn't followed by a P.'''
	if not end_pos:
		end_pos = len(sequence)
	break_next_flag = False # Says we've found the end of our desired sequence and want to keep one more peptide (in case of missed cleavage)
	tryptic_peptides = []
	seq = ''
	prev_peptide = (0,'*')
	for i in range(len(sequence)):
		letter = sequence[i]
		seq += letter
		try:
			next_letter = sequence[i+1]
		except IndexError:
			next_letter = '*'
		if ((letter == 'K' or letter == 'R') and next_letter != 'P') or letter == '*': #should only find a stop codon if we're passing it a fully-translated sequence including stops, as in junctions
			if pos <= i and prev_peptide[0] < pos and allow_missed_cleavage: 
				# This is the first peptide of our desired sequence. Keep the previous one around in case of missed cleavage.
				if prev_peptide[1][-1] != '*':  # If it ends in a stop codon we don't want it
					tryptic_peptides.append(prev_peptide[1])
			if pos <= i:
				tryptic_peptides.append(seq)
			prev_peptide = (i,seq) # Reset the previous peptide and its end position
			seq = ''
			if break_next_flag:
				# We just added a peptide after the end of our desired sequence and can now finish
				break_next_flag = False
				break
			if i >= end_pos:
				if not allow_missed_cleavage or prev_peptide[1][-1] == '*':
					break
				else:
					break_next_flag = True
	if seq != '' and pos <= i:
		tryptic_peptides.append(seq)
	return tryptic_peptides

def get_powerset(vars):
	'''Returns all unique combinations of variants in a peptide.'''
	
	# I want to keep consistent positions, so taking things out of the dictionary for now
	save_varkeys = vars.keys()
	save_varvals = []
	for key in save_varkeys:
		save_varvals.append(vars[key])
	
	# Add a blank to each
	for v in save_varvals:
		v.append([])
	
	# Get the full powerset, including blanks
	full_powerset = list(product(*save_varvals))
	
	variant_sets = []
	for variant in full_powerset:
		new_set = []
		for i in range(len(variant)):
			if variant[i] != []:
				new_set.append([save_varkeys[i], variant[i][0], variant[i][1]])
		if new_set != []:
			new_set = sorted(new_set, key=lambda x: x[0])
			# May or may not work, haven't tested yet
			repeat_flag = False
			for i in range(len(new_set)):
				if new_set[i][2] == '*' and i != len(new_set)-1:
					repeat_flag = True
			if not repeat_flag:		
				variant_sets.append(new_set)
	
	return variant_sets

def write_peptides(gene, vars, peptide_start_pos, peptide, out_file):
	'''For each peptide, writes out a copy with all possible combinations of variants.'''
	var_combos = get_powerset(vars)
	new_peptides = [[peptide, []]]
	for var_set in var_combos:
		new_pep = peptide
		var_string = []
		for var in var_set:
			pos = var[0]-peptide_start_pos-1
			new_pep = new_pep[:pos]+var[2]+new_pep[pos+1:]
			var_string.append(var[1]+str(var[0])+var[2])
		pep_to_write = new_pep.split('*')[0]
		# Letting there be a missed cleavage just in case the variant changes the cleaving.
		if 25 >= len(pep_to_write) >= 6 and len(trypsinize(pep_to_write)) <= 2:
			out_file.write('>%s START:%d END:%d VAR:%s\n' % (gene, peptide_start_pos+1, peptide_start_pos+len(new_pep), ','.join(var_string)))
			out_file.write('%s\n' % pep_to_write)
		new_peptides.append([new_pep, var_string])
	return new_peptides

def write_missed_cleavage_peptides(gene, prev_peptides, cur_peptides, cur_start_pos, out_file):
	'''Allows for a single missed cleavage - basically does a Cartesian join between all variants of the previous peptide and all variants of the current peptide (including the case where one, but not both, has no variants)'''
	if prev_peptides != [] and cur_peptides != []:
		prev_peptide_start = cur_start_pos-len(prev_peptides[0][0])
		cur_peptide_end = cur_start_pos+len(cur_peptides[0][0])-1
		for i in prev_peptides:
			# Make sure we don't get two missed cleavages due to an added R/K variant by trypsinizing again.
			if '*' not in i[0]:
				for j in cur_peptides:
					if i[1] != [] or j[1] != []:
						vars = i[1]+j[1]
						pep_to_write = (i[0]+j[0]).split('*')[0]
						if 25 >= len(pep_to_write) >= 6 and len(trypsinize(pep_to_write)) <=2 :
							out_file.write('>%s START:%d END:%d (missed cleavage after %d) VAR:%s\n' % (gene, prev_peptide_start+1, cur_peptide_end+1, cur_start_pos, ','.join(vars)))
							out_file.write('%s\n' % pep_to_write)

def assign_variants(gene, tryptic_peptides, variants, out_file, no_missed_cleavage):
	'''Figures out which variants belong to each tryptic peptide then calls the function that writes the peptide out to the file'''
	# Get positions and original/changed AAs in a workable format
	var_set = set(variants) # remove duplicates
	vars = {}
	for var in var_set:
		pos = int(re.findall(r'\d+', var)[0])
		tmp = vars.get(pos, [])
		tmp.append(var.split(str(pos)))
		vars[pos] = tmp
	# For each tryptic peptide, decide which variants belong and then call a function to write them out
	total_aas = 0
	prev_peptides = []
	# The last one is an extra peptide in case the stop codon gets erased.
	for i in range(len(tryptic_peptides)-1):
		peptide = tryptic_peptides[i]
		vars_in_peptide = [var for var in vars.keys() if total_aas <= var-1 < total_aas+len(peptide)]
		# Make sure to get the stop codon if that's a variant, and add the next chunk onto the peptide to be written out if the stop codon gets erased
		for var in vars.keys():
			if vars[var][0][0] == '*' and var-1==total_aas+len(peptide):
				vars_in_peptide.append(var)
				peptide += '*'+tryptic_peptides[i+1]
		cur_peptides = write_peptides(gene, {v: vars[v] for v in vars_in_peptide}, total_aas, peptide, out_file)
		if not no_missed_cleavage:
			# If we are allowing missed cleavages...
			write_missed_cleavage_peptides(gene, prev_peptides, cur_peptides, total_aas, out_file)
		prev_peptides = cur_peptides
		total_aas += len(peptide)

def make_aa_peptide_fasta(log_dir, no_missed_cleavage):
	'''Main function for the tryptic peptide fasta file.'''
	# Grab all of the amino acid variants
	vars = {}
	f = open(log_dir+'proteome.bed.aa.var','r')
	line = f.readline()
	while line:
		vars[line.split('\t')[0]] = line.rstrip().split('\t')[-1].split(',')
		line = f.readline()
	f.close()
	
	f = open(log_dir+'proteome.bed.dna','r')
	tryp_fasta = open(log_dir+'tryptic_proteome.fasta','w')
	sequence = ''
	extra_seq = ''
	line = f.readline()
	while line:
		# Line type one: First gene header line
		if line[:3] == 'chr':
			# Finish processing the previous gene
			# Using a try/except block here to deal with the possibility of this being the first line
			# There has to be a more graceful way to do that, right?
			try:
				if gene in vars:
					translated = translate_seq(sequence.upper(), strand) # no variants here yet
					tryptic_peptides = trypsinize(translated)
					tryp_tran_extra = trypsinize(translate_seq(extra_seq.upper(), strand))
					if len(tryp_tran_extra) == 0:
						tryp_tran_extra = ['*']
					tryptic_peptides.append(tryp_tran_extra[0]) # add one extra peptide in case the last stop codon gets erased
					assign_variants(gene, tryptic_peptides, vars[gene], tryp_fasta, no_missed_cleavage)
			except UnboundLocalError:
				pass
			# Start processing the new gene
			sequence = ''
			extra_seq = ''
			strand = line.split('\t')[5]
			gene = line.split('\t')[3]
		# Line type two: Second gene header line
		elif line[0] == '>':
			second_header = line # This will be the header in the fasta file
		# Line type three: Exon line (can be pre-100, post-100 or internal)
		else:
			spline = line.rstrip().split('\t')
			if spline[1] != '-1' and spline[1] != '+1':
				sequence += spline[-1]
			elif spline[1] == '-1' and spline[5] == '-':
				extra_seq = spline[-1]
			elif spline[1] == '+1' and spline[5] =='+':
				extra_seq = spline[-1]
		line = f.readline()
	# And finish the last one
	try:
		if gene in vars:
			translated = translate_seq(sequence.upper(), strand) # no variants here yet
			tryptic_peptides = trypsinize(translated)
			tryp_tran_extra = trypsinize(translate_seq(extra_seq.upper(), strand))
			if len(tryp_tran_extra) == 0:
				tryp_tran_extra = ['*']
			tryptic_peptides.append(tryp_tran_extra[0]) # add one extra peptide in case the last stop codon gets erased
			assign_variants(gene, tryptic_peptides, vars[gene], tryp_fasta, no_missed_cleavage)
	except UnboundLocalError:
		pass
		
	f.close()
	tryp_fasta.close()

def make_indel_peptide_fasta(log_dir, logfile):
	'''Main function for adding indels to the tryptic peptide fasta.
	Basically, add the peptide that includes the indel and everything that comes after it.
	Unless the indel changes the length by a non-frameshift amount,
	in which case just save the portion that includes the indel'''
	# Grab the indel variants
	vars = {}
	f = open(log_dir+'proteome.indel.var.bed.dna.fasta','r')
	tryp_fasta = open(log_dir+'tryptic_proteome.fasta','a') # will append to the file
	line = f.readline()
	while line:
		indel = line.split()[-1].split(':')[0]
		gene = line.split()[0]
		nt_pos = int(re.findall(r'\d+', indel)[0]) # The DNA position, starting from zero
		strand = line.split()[1][-1]
		to_delete = indel.split(str(nt_pos))[0].split('-')[1]
		to_insert = indel.split(str(nt_pos))[1]
		seq = f.readline().rstrip() # Next line is the sequence
		if strand == '+':
			aa_pos = nt_pos/3 # Again, indexes from zero
		else:
			aa_pos = len(seq) - nt_pos/3
		tryptic_peptides = trypsinize(seq, aa_pos)
		if (len(to_insert) - len(to_delete))%3 == 0: # no change in frame, so only the indel is changed
			if len(tryptic_peptides) > 0:
				tryptic_peptides = [tryptic_peptides[0]]
			else:
				write_to_log("No tryptic peptides: %s, %s, %s, %s" % (gene[1:], indel, seq, strand), logfile)
		for i in range(len(tryptic_peptides)):
			tp = tryptic_peptides[i]
			if 25 >= len(tp) >= 6:
				tryp_fasta.write("%s\tPeptide %d\n" % (line.rstrip(), i))
				tryp_fasta.write("%s\n" % tp)
		line = f.readline() # Next line is the header
	
	f.close()
	tryp_fasta.close()

### These functions are used to merge junction files.

def pull_junc_files(junc_dir):
	'''Finds all *.txt files in the directory.'''
	files = os.listdir(junc_dir)
	junc_files = []
	for f in files:
		if f.endswith('.txt'):
			junc_files.append(f)
	return junc_files

def pull_junc_bed_files(junc_dir):
	'''Finds all *.bed files in the directory.'''
	files = os.listdir(junc_dir)
	junc_files = []
	for f in files:
		if f.endswith('.bed'):
			junc_files.append(f)
	return junc_files

def convert_tophat_to_mapsplice(junc_dir):
	# Tries to pull a list of junction files. Returns with a warning if unable to find any.
	if not os.path.isdir(junc_dir):
		# No junction files are going to be found. Gotta leave the function.
		warnings.warn("Unable to open junction folder %s. Giving up on finding splice junctions." % junc_dir)
		return 1
	junc_files = pull_junc_bed_files(junc_dir)
	if len(junc_files) == 0:
		warnings.warn("Unable to find any .bed files in %s." % junc_dir)
		return 1
	w = open(junc_dir+'/junctions.txt')
	for fil in junc_files:
		f = open(fil,'r')
		for line in f.readlines():
			spline = line.split()
			new_start = str(int(spline[1])+int(spline[10].split(',')[0]))
			new_end = str(int(spline[2])-int(spline[10].split(',')[1])+1)
			spline[1] = new_start
			spline[6] = new_start
			spline[2] = new_end
			spline[7] = new_end
			new_len = str(new_end-new_start+int(spline[11].split(',')[1]))
			spline[11] = spline[11].split(',')[0]+','+new_len
			w.write('\t'.join(spline)+'\n')
		f.close()
	w.close()

def merge_junction_files(junc_dir, log_dir):
	'''Merges junction files found in junc_dir. Returns 1 if unable to open any junction files or nothing if it runs to completion.'''
	
	# May or may not allow this to change. Not a thousand percent sure what it is.
	min_length = 8
	min_sum_length = 45
	
	num_observ_map = {}
	exon_len_map = {}
	
	# Tries to pull a list of junction files. Returns with a warning if unable to find any.
	if not os.path.isdir(junc_dir):
		# No junction files are going to be found. Gotta leave the function.
		warnings.warn("Unable to open junction folder %s. Giving up on finding splice junctions." % junc_dir)
		return 1
	junc_files = pull_junc_files(junc_dir)
	if len(junc_files) == 0:
		warnings.warn("Unable to find any .txt files in %s." % junc_dir)
		return 1
	
	for file in junc_files:
		label = ''
		if '-' in file:
			label = '-'+file.split('-')[0]
			
		f = open(junc_dir+'/'+file,'r')
		for line in f.readlines():
			spline = line.split()
			if len(spline) > 11:
				chrm, end_ex_1, begin_ex_2, junc_num, num_observ, len_exons, begin_exons = spline[0], int(spline[1]), int(spline[2]), spline[3], int(spline[4]), spline[10].rstrip(','), spline[11].rstrip().rstrip(',')
				splen_exons = map(int,len_exons.split(','))
				spbegin_exons = map(int,begin_exons.split(','))
				begin_intron = end_ex_1
				end_intron = begin_ex_2
				spbegin_exons[1] = spbegin_exons[1]-splen_exons[1]+splen_exons[0]-2 # -2 for some weird read_chr_bed nonsense
				if 'chr' not in chrm:
					chrm = 'chr'+chrm
				if splen_exons[0] >= min_length and splen_exons[1] >= min_length and sum(splen_exons) >= min_sum_length:
					#key = '%s#%d#%d' % (chrm,begin_intron,end_intron)
					num_observ_map.setdefault(chrm,{})
					num_observ_map[chrm].setdefault(begin_intron,{})
					num_observ_map[chrm][begin_intron][end_intron] = num_observ_map[chrm][begin_intron].get(end_intron,0) + int(num_observ)
					exon_len_map.setdefault(chrm,{})
					exon_len_map[chrm].setdefault(begin_intron,{})
					exon_len_map[chrm][begin_intron][end_intron] = [splen_exons, spbegin_exons, junc_num]

		f.close()
		
	# Time to write out our results
	w = open(log_dir+'/merged-junctions.bed','w')
	count = 0
	for chr in sorted(num_observ_map.keys()):
		for beg_int in sorted(num_observ_map[chr].keys()):
			for end_int in sorted(num_observ_map[chr][beg_int].keys()):
				[lens, begins, junc_num] = exon_len_map[chr][beg_int][end_int]
				begin_ex_1 = beg_int-lens[0]
				end_ex_2 = end_int+lens[1]
				start_ex_2 = end_int-begin_ex_1
				if start_ex_2 > 0:
					w.write("%s\t%d\t%d\t%s\t%d\t+\t%d\t%d\t0\t2\t%d,%d\t%d,%d\n" % (chr,begin_ex_1,end_ex_2,junc_num,num_observ_map[chr][beg_int][end_int],begin_ex_1,end_ex_2, lens[0], lens[1], begins[0], begins[1]))
					count += 1
	w.close()
							
def filter_known_transcripts(transcriptome_bed, results_folder, logfile):
 	'''Filters out spliceform transcripts that are known from those that aren't, and annotates the ones that can be annotated. Haven't done any serious testing on this yet, so...yeah, use at yr own risk'''

	junctions_model = {}
	junctions_alternative = {}
	junctions_pos_donors = {}
	junctions_neg_donors = {}

	try:
		f = open(transcriptome_bed,'r')
	except IOError:
		raise_warning("Could not find transcriptome file at %s!\n Treating all junctions as novel." % transcriptome_bed)
	
	if f:	
		line_number = 0 # Might never get used
		for line in f.readlines():
			spline = line.split('\t')
			chrm, begin_exon, trans_id, strand, block_count, block_sizes, block_starts = spline[0], int(spline[1]), spline[3], spline[5], int(spline[9]), spline[10].rstrip(','), spline[11].rstrip().rstrip(',')
			sizes = map(int, block_sizes.split(','))
			starts = map(int, block_starts.split(','))
			start = begin_exon

			#Inclusive of (start+block lengths), EXCLUSIVE of (start+block starts) even 0

			for i in range(0,block_count-1):
				end_ex_1 = begin_exon+starts[i]+sizes[i]
				beg_ex_2 = begin_exon+starts[i+1]+1
				
				# Set up the normal splice
				junctions_model.setdefault(chrm,{})
				junctions_model[chrm].setdefault(end_ex_1,{})
				junctions_model[chrm][end_ex_1].setdefault(beg_ex_2,[]).append(trans_id) # this is saving a lot of info :/
				
				# Set up all alternative splices
				for j in range(i+1, block_count-1):
					new_end = begin_exon+starts[j]+1
					junctions_alternative.setdefault(chrm,{})
					junctions_alternative[chrm].setdefault(end_ex_1,{})
					junctions_alternative[chrm][end_ex_1].setdefault(new_end,[]).append(trans_id)
				
				if strand == '-':
					junctions_neg_donors.setdefault(chrm,{})
					junctions_neg_donors[chrm].setdefault(beg_ex_2,[]).append(trans_id)
				else:
					junctions_pos_donors.setdefault(chrm,{})
					junctions_pos_donors[chrm].setdefault(end_ex_1,[]).append(trans_id)
				
				start = beg_ex_2
			line_number += 1
		f.close()
	
	f = open(results_folder+'/merged-junctions.bed','r')
	w = open(results_folder+'/merged-junctions.known.bed.log','w')
	w2 = open(results_folder+'/merged-junctions.novel.bed','w')
	w3 = open(results_folder+'/merged-junctions.alt.bed','w')
	w4 = open(results_folder+'/merged-junctions.donor.bed','w')
	matches = 0
	non_matches = 0
	alt_match = 0
	donor_match = 0
	
	line = f.readline()
	while line:
		spline = line.split()
		if len(spline) > 11:
			chrm, beg_ex_1, end_ex_2, junc_num, num_observ, len_exons, begin_exons = spline[0], int(spline[1]), int(spline[2]), spline[3], int(spline[4]), spline[10], spline[11]
			splen_exons = map(int,len_exons.split(','))
			#spbegin_exons = map(int,begin_exons.split(','))
			if num_observ >= 1:	
				# Figure out whether it already exists in the database
				end_ex_1 = beg_ex_1+splen_exons[0]
				beg_ex_2 = end_ex_2-splen_exons[1]
				if chrm in junctions_model.keys() and end_ex_1 in junctions_model[chrm].keys():
					if beg_ex_2 in junctions_model[chrm][end_ex_1]:
						# This is a known splice site. Writing out the bed file just to have a record - won't do anything with it
						matches += 1
						spline[3] = spline[3]+'-'+'-'.join(junctions_model[chrm][end_ex_1][beg_ex_2])
						w.write('\t'.join(spline)+'\n')
					elif end_ex_1 in junctions_alternative[chrm].keys() and beg_ex_2 in junctions_alternative[chrm][end_ex_1]:
						# This is an unannotated splice site, but both donor and acceptor sites 
						# are known (so it's just skipped a couple of exons in one gene)
						alt_match += 1
						spline[3] = spline[3]+'-'+'-'.join(junctions_alternative[chrm][end_ex_1][beg_ex_2])
						w3.write('\t'.join(spline)+'\n')
					elif chrm in junctions_pos_donors.keys() and end_ex_1 in junctions_pos_donors[chrm]:
						# The donor splice site is known but the acceptor is completely new (pos strand)
						donor_match += 1
						spline[3] = spline[3]+'-'+'-'.join(junctions_pos_donors[chrm][end_ex_1])
						w4.write('\t'.join(spline)+'\n')
				elif chrm in junctions_neg_donors.keys() and beg_ex_2 in junctions_neg_donors[chrm]:
					# The donor splice site is known but the acceptor is completely new (neg strand)
					donor_match += 1
					spline[3] = spline[3]+'-'+'-'.join(junctions_neg_donors[chrm][beg_ex_2])
					w4.write('\t'.join(spline)+'\n')
				else:
					# Completely novel. Change the name and go ahead.
					nov_name = "NO_GENE-%d-%s-%d-%d-%s" % (num_observ, chrm, end_ex_1, beg_ex_2, junc_num)
					spline[3] = nov_name
					spline[8] = "0,0,0"
					spline[5] = '+'
					non_matches += 1
					w2.write('\t'.join(spline)+'\n')
					
		line = f.readline()
				
	f.close()
	w.close()
	w2.close()
	w3.close()
	w4.close()
	
	write_to_log("\n---Transcript Filtering Results---",logfile)
	write_to_log("Total known junctions: %d" % matches, logfile)
	write_to_log("Candidate novel junctions: %d" % non_matches, logfile)
	write_to_log("Alternative splice junctions with known endpoints: %d" % alt_match, logfile)
	write_to_log("Alternative splice junctions with known donor site: %d" % donor_match, logfile)

def filter_alternative_splices(log_dir, threshA, threshAN, threshN, logfile):
	
	# Get translations from transcript ID to protein ID
	transmap = {}
	f = open(log_dir+'/proteome-genes.txt','r')
	for line in f.readlines():
		spline = line.split('\t')
		transmap[spline[1]] = spline[0]
	f.close()
	
	# Grab each line from proteome.bed
	prot_lines = {}
	f = open(log_dir+'/proteome.bed','r')	
	for line in f.readlines():
		spline = line.split('\t')
		prot_lines[spline[3]] = line
	f.close()
	
	untranslated = open(log_dir+'/merged-junctions.untranslated.bed.log','w')
	under_threshold = open(log_dir+'/merged-junctions.under.threshold.bed.log','w')
	
	# Start with the known endpoints
	w = open(log_dir+'/merged-junctions.alt.filtered.bed','w')
	f = open(log_dir+'/merged-junctions.alt.bed','r')
	line = f.readline()
	while line:
		spline = line.split()
		if int(spline[4]) < threshA:
			under_threshold.write(line)
		else:
			junc_ids = spline[3].split('-')[1:]
			for id in junc_ids:
				if id in transmap.keys() and transmap[id] in prot_lines.keys():
					to_write = write_known_endpoints(line, prot_lines[transmap[id]])
					if to_write: # everything works as planned
						w.write(to_write+'\n')
					else: # one end of the junction is in an untranslated region
						spline[3] = spline[3].split('-')[0]+'-'+id
						untranslated.write('\t'.join(spline)+'\n')
				else: # this transcript is untranslated mRNA
					spline[3] = spline[3].split('-')[0]+'-'+id
					untranslated.write('\t'.join(spline)+'\n')
		line = f.readline()
	f.close()
	w.close()
	
	# Next: the conserved donor sites (this'll be harder)
	w = open(log_dir+'/merged-junctions.donor.filtered.bed','w')
	f = open(log_dir+'/merged-junctions.donor.bed','r')
	line = f.readline()
	while line:
		spline = line.split()
		if int(spline[4]) < threshAN:
			under_threshold.write(line)
		else:
			junc_ids = spline[3].split('-')[1:]
			for id in junc_ids:
				if id in transmap.keys() and transmap[id] in prot_lines.keys():
					to_write = write_donor_endpoint(line, prot_lines[transmap[id]])
					if to_write: # everything works as planned
						w.write(to_write+'\n')
					else: # the donor end of the junction is in an untranslated region
						spline[3] = spline[3].split('-')[0]+'-'+id
						untranslated.write('\t'.join(spline)+'\n')
				else: # this transcript is untranslated mRNA
					spline[3] = spline[3].split('-')[0]+'-'+id
					untranslated.write('\t'.join(spline)+'\n')
		line = f.readline()
	f.close()
	w.close()
	
	# Finally: novels (this only removes those under the read count threshold)
	w = open(log_dir+'/merged-junctions.novel.filtered.bed','w')
	f = open(log_dir+'/merged-junctions.novel.bed','r')
	line = f.readline()
	while line:
		spline = line.split()
		if int(spline[4]) < threshN:
			under_threshold.write(line)
		else:
			w.write(line)
		line = f.readline()
	f.close()
	w.close()
	
	untranslated.close()
	under_threshold.close()

def write_known_endpoints(junc_line, prot_line):
	'''Writing out alternate splice junctions where both donor and acceptor sites are known sites
	(essentially, we're just skipping a couple exons)'''
	
	# Get information from the junction bed
	spjunc = junc_line.split()
	junc_name = spjunc[3].split('-')[0]
	num_observ = int(spjunc[4])
	len_exons = spjunc[10]
	splen_exons = map(int,len_exons.split(','))
	beg_ex_1 = int(spjunc[1])
	end_ex_2 = int(spjunc[2])
	end_ex_1 = beg_ex_1+splen_exons[0]
	beg_ex_2 = end_ex_2-splen_exons[1]

	# Get information from the proteome bed
	spline = prot_line.split()
	chr, start, end, prot_id, strand, block_count, block_sizes, block_starts = spline[0], int(spline[1]), int(spline[2]), spline[3], spline[5], int(spline[9]), spline[10].rstrip(','), spline[11].rstrip().rstrip(',')
	spblock_sizes = block_sizes.split(',')
	spblock_starts = block_starts.split(',')
	int_block_sizes= map(int,spblock_sizes)
	int_block_starts = map(int,spblock_starts)
	
	# Calculate: donor exon #, acceptor exon #, new block count/sizes/starts
	donor = None
	acceptor = None
	for i in range(0,block_count):
		if start+int_block_sizes[i]+int_block_starts[i] == end_ex_1:
			donor = i+1
		elif start+int_block_starts[i]+1 == beg_ex_2:
			acceptor = i+1
			
	if donor and acceptor:	
		new_block_count = block_count-acceptor+donor+1
		new_block_sizes = ','.join(spblock_sizes[:donor] + spblock_sizes[acceptor-1:])
		new_block_starts = ','.join(spblock_starts[:donor] + spblock_starts[acceptor-1:])
	
		if strand=='-':
			tmp_donor = block_count-acceptor+1
			acceptor = block_count-donor+1
			donor = tmp_donor
			
		new_junc_name = "-".join([prot_id, str(donor), str(acceptor), str(num_observ), str(end_ex_1), str(beg_ex_2), junc_name])
		to_write = '\t'.join([chr, str(start), str(end), new_junc_name, str(num_observ), strand, str(start), str(end), "0,0,255", str(new_block_count), new_block_sizes, new_block_starts])
		return to_write
	
	else:
		# This splice site goes into an untranslated region of the protein. We're gonna toss it.
		return None

def write_donor_endpoint(junc_line, prot_line):
	'''When we have the donor site conserved but a new acceptor site, things get more complicated.'''
	
	# Get information from the junction bed
	spjunc = junc_line.split()
	junc_name = spjunc[3].split('-')[0]
	num_observ = int(spjunc[4])
	len_exons = spjunc[10]
	splen_exons = map(int,len_exons.split(','))
	beg_ex_1 = int(spjunc[1])
	end_ex_2 = int(spjunc[2])
	end_ex_1 = beg_ex_1+splen_exons[0]
	beg_ex_2 = end_ex_2-splen_exons[1]
	
	# Get information from the proteome bed
	spline = prot_line.split()
	chr, start, end, prot_id, strand, block_count, block_sizes, block_starts = spline[0], int(spline[1]), int(spline[2]), spline[3], spline[5], int(spline[9]), spline[10].rstrip(','), spline[11].rstrip().rstrip(',')
	spblock_sizes = block_sizes.split(',')
	spblock_starts = block_starts.split(',')
	int_block_sizes= map(int,spblock_sizes)
	int_block_starts = map(int,spblock_starts)
	
	# Calculate: donor exon #, acceptor exon #, new block count/sizes/starts
	low_exon = None
	high_exon = None
	for i in range(0,block_count):
		if start+int_block_sizes[i]+int_block_starts[i] == end_ex_1:
			low_exon = i+1
		elif start+int_block_starts[i]+1 == beg_ex_2:
			high_exon = i+1
	
	if strand == '+' and end_ex_2 > end:
		# This is actually easier - we take all blocks before beg_ex_1 
		
		if not low_exon:
			# The junction starts in an untranslated region
			return None
		
		end = end_ex_2
		
		new_block_count = low_exon+1
		new_block_sizes = ','.join(spblock_sizes[:low_exon]) + ',' + str(splen_exons[1])
		new_block_starts = ','.join(spblock_starts[:low_exon]) + ',' + str(beg_ex_2-start-1)
	
		new_junc_name = "-".join([prot_id, str(low_exon), str(num_observ), str(end_ex_1), str(beg_ex_2), junc_name])
		to_write = '\t'.join([chr, str(start), str(end), new_junc_name, str(num_observ), strand, str(start), str(end), "0,100,0", str(new_block_count), new_block_sizes, new_block_starts])
		
	elif strand == '-' and beg_ex_1 < start:
		# We have a negative strand donor - we just take all blocks after "end_ex_2"
		
		if not high_exon:
			# The junction starts in an untranslated region
			return None
		
		add_to_start = start-beg_ex_1
		start = beg_ex_1
		
		new_block_count = block_count - high_exon + 2
		new_block_sizes = str(splen_exons[0])+','+','.join(spblock_sizes[high_exon-1:])
		int_block_starts = [x+add_to_start for x in int_block_starts]
		spblock_starts = map(str,int_block_starts)
		new_block_starts = "0,"+','.join(spblock_starts[high_exon-1:])
		new_junc_name = "-".join([prot_id, str(block_count-high_exon+1), str(num_observ), str(end_ex_1), str(beg_ex_2), junc_name])
		
		to_write = '\t'.join([chr, str(start), str(end), new_junc_name, str(num_observ), strand, str(start), str(end), "0,100,0", str(new_block_count), new_block_sizes, new_block_starts])
		
	else:
		# The toughest one - we need to resume in the middle of a block
		if (strand == '-' and not high_exon) or (strand == '+' and not low_exon):
			# The junction starts in an untranslated region
			return None
		if strand == "+":
			for i in range(low_exon, block_count+1):
				if start+int_block_starts[i] + int_block_sizes[i] >= beg_ex_2:
					# I think we can simplify it all the way down to this!?
					next_block_num = i+1
					new_block_count = block_count-next_block_num+low_exon+1
					new_block_size = start + int_block_starts[i] + int_block_sizes[i] - beg_ex_2
					new_block_start = beg_ex_2-start-1
					break
			
			new_block_sizes = spblock_sizes[:low_exon] + [str(new_block_size)] + spblock_sizes[next_block_num:]
			new_block_starts = spblock_starts[:low_exon] + [str(new_block_start)] + spblock_starts[next_block_num:]
			new_junc_name = "-".join([prot_id, str(low_exon), str(num_observ), str(end_ex_1), str(beg_ex_2), junc_name])
		else:
			for i in range(high_exon-2,-1,-1):
				if start+int_block_starts[i] < end_ex_1:
					next_block_num = i+1
					new_block_count = block_count-high_exon+next_block_num+1
					new_block_size = end_ex_1 - start - int_block_starts[i]
					break
		
			new_block_sizes = spblock_sizes[:next_block_num-1] + [str(new_block_size)] + spblock_sizes[high_exon-1:]
			new_block_starts = spblock_starts[:next_block_num] + spblock_starts[high_exon-1:]
			donor = block_count-high_exon+1
			new_junc_name = "-".join([prot_id, str(donor), str(num_observ), str(end_ex_1), str(beg_ex_2), junc_name])
			
		to_write = '\t'.join([chr, str(start), str(end), new_junc_name, str(num_observ), strand, str(start), str(end), "0,100,0", str(new_block_count), ','.join(new_block_sizes), ','.join(new_block_starts)])
	
	return to_write
		
### These functions are used to translate junction files.

def translate_novels_with_trypsin(log_dir, bed_file, logfile):
	'''NOT IN USE CURRENTLY
	Main function for adding non-frameshifted junctions to the tryptic peptide fasta.
	Basically, add the peptide that includes the junction and everything that comes after it.
	Later: maybe remove the parts that come after it if it isn't frameshifted? That might be complicated.'''
	# Get the descriptions
	desc = {}
	f = open(log_dir+"proteome-descriptions.txt",'r')
	line = f.readline()
	while line:
		spline = line.rstrip().split('\t')
		desc[spline[0]] = spline[1]
		line = f.readline()
	f.close()
	# Get the gene abbreviations
	abbr = {}
	f = open(log_dir+"proteome-genes.txt",'r')
	line = f.readline()
	while line:
		spline = line.rstrip().split()
		abbr[spline[0]] = spline[-1]
		line = f.readline()
	f.close()
	
	# The rest of it
	f = open(log_dir+bed_file,'r')
	out_fasta = open(log_dir+bed_file+".fasta",'w') # The original QUILTS writes two other files but they're just duplicates of proteome.aa.var.bed and proteome.aa.var.bed.dna, it seems. For now I'm leaving them out.
	line = f.readline()
	sequence = ""
	prev_AA_subst = []
	prev_gene = ''
	while line:
		# Line type one: First gene header line
		if line[:3] == 'chr':
			# Finish processing the previous gene - focus first on proteome.bed.aa.var
			# Using a try/except block here to deal with the possibility of this being the first line
			# There has to be a more graceful way to do that, right?
			if "Error Reading Sequence" in sequence: # This error comes from read_chr_bed, which I did not write and am afraid to touch.
				write_to_log("Error in sequence for %s. Skipping..." % gene, logfile)
			else:
				try:
					if prev_gene != gene:
						prev_AA_subst = []
					if AA_subst not in prev_AA_subst:
						for i in range(3): # Gotta get rid of the 600 preceding
							translated = translate_seq(sequence.upper()[i:], '+', True)
							if i==1:
								# The splice happens between codons
								tryp = trypsinize(translated, 215, 216, True) # this math may be bad
							else:
								# The splice happens in the middle of a codon
								tryp = trypsinize(translated, 216, 216, True) # this math may be bad
							to_write = ''.join(tryp).rstrip('*')
							if len(to_write) > 6:
								out_fasta.write(format_header('%s-%s%s\n' % (second_header.rstrip(), str(i), '+'), 'juncN', None, None))
								out_fasta.write(to_write+'\n')
							translated = translate_seq(sequence.upper()[:len(sequence)-i], '-', True)
							if i==2:
								# Splice happens between codons
								tryp = trypsinize(translated, 214, 215, True) # this math may be bad
							else:
								# Splice happens in the middle of a codon
								tryp = trypsinize(translated, 215, 215, True) # this math may be bad
							to_write = ''.join(tryp).rstrip('*')
							to_write = to_write.split('*')[0]
							#if '>NP_001299602-AN2-3-3' in second_header:
							#	print translated, tryp, to_write, translated[214:216]
							if len(to_write) > 6:
								#out_fasta.write('%s-%s%s\n' % (second_header.rstrip(), str(i), '-'))
								out_fasta.write(format_header('%s-%s%s\n' % (second_header.rstrip(), str(i), '-'), 'juncN', None, None))
								out_fasta.write(to_write+'\n')
					prev_AA_subst.append(AA_subst)
					prev_gene = gene
				except UnboundLocalError:
					pass
			# Start processing the new gene
			sequence = ""
			strand = line.split('\t')[5]
			gene = line.split('\t')[3].split('-')[0]
		# Line type two: Second gene header line
		elif line[0] == '>':
			second_header = line # This will be the header in the fasta file
			AA_subst = line.rsplit('-',1)[-1].split(':')[0]
		# Line type three: Exon line (can be pre-100, post-100 or internal)
		else:
			spline = line.rstrip().split('\t')
			# Here we want to keep both the preceding and ...anteceding? sequence chunks
			sequence += spline[-1]
		line = f.readline()
		
	# And finish the last one
	if "Error Reading Sequence" in sequence: # This error comes from read_chr_bed, which I did not write and am afraid to touch.
		write_to_log("Error in sequence for %s. Skipping..." % gene, logfile)
	else:
		try:
			if prev_gene != gene or AA_subst != prev_AA_subst:
				for i in range(3): # Gotta get rid of the 600 preceding
					translated = translate_seq(sequence.upper()[i:], '+', True)
					if i==1:
						# The splice happens between codons
						tryp = trypsinize(translated, 215, 216, True) # this math may be bad
					else:
						# The splice happens in the middle of a codon
						tryp = trypsinize(translated, 216, 216, True) # this math may be bad
					to_write = ''.join(tryp).rstrip('*')	
					if len(to_write) > 6:
						#out_fasta.write('%s-%s%s\n' % (second_header.rstrip(), str(i), '+'))
						out_fasta.write(format_header('%s-%s%s\n' % (second_header.rstrip(), str(i), '+'), 'juncN', None, None))
						out_fasta.write(to_write+'\n')
					translated = translate_seq(sequence.upper()[:len(sequence)-i], '-', True)
					if i==2:
						# Splice happens between codons
						tryp = trypsinize(translated, 214, 215, True) # this math may be bad
					else:
						# Splice happens in the middle of a codon
						tryp = trypsinize(translated, 215, 215, True) # this math may be bad
					to_write = ''.join(tryp).rstrip('*')
					to_write = to_write.split('*')[0] # Not quite sure why, but was still getting stops in the middle of the output sequences.
					if len(to_write) > 6:
						#out_fasta.write('%s-%s%s\n' % (second_header.rstrip(), str(i), '-'))
						out_fasta.write(format_header('%s-%s%s\n' % (second_header.rstrip(), str(i), '-'), 'juncN', None, None))
						out_fasta.write(to_write+'\n')
		except UnboundLocalError:
			pass

	f.close()
	out_fasta.close()

def translate_novels(log_dir, bed_file, logfile):
	'''In the bed file we read in, the sequences come in pairs - first sequence in the pair is the left gene, second is the right gene.'''
	
	translate_table = maketrans("ACGTacgt","TGCAtgca")
	
	f = open(log_dir+bed_file,'r')
	out_fasta = open(log_dir+bed_file+".fasta",'w') # The original QUILTS writes two other files but they're just duplicates of proteome.aa.var.bed and proteome.aa.var.bed.dna, it seems. For now I'm leaving them out.
	
	line = 'a'
	
	while line:
		line = f.readline() # header 1
		second_header = f.readline() # header 2
		restart_pointer = f.tell() # if it turns out those headers didn't find a sequence in the bed file, we can restart from here in the next line
		line = f.readline() # buffer
		upstream = f.readline() # upstream sequence
		downstream = f.readline() # downstream sequence
		line = f.readline() # buffer
		
		if not upstream: # we're done and are going to error out if we don't break now
			break
		
		if upstream[0] == '>': # this sequence wasn't found, go backwards. might break for the last entry...
			f.seek(restart_pointer)		
		else: # we're good to go
			seq = upstream.rstrip().split()[-1] + downstream.rstrip().split()[-1]
			seq = seq.upper()
			for i in range(3):
				# Positive strand
				translated = translate_seq(seq[i:],'+',True)
				to_write = translated.split('*')[0]
				if len(to_write) > 16: #must include the junction peptide(s)
					out_fasta.write(format_header('%s-%s%s\n' % (second_header.rstrip(), str(i), '+'), 'juncN', None, None))
					out_fasta.write(to_write+'\n')
				# Negative strand
				translated = translate_seq(seq[:len(seq)-i], '-', True)
				to_write = translated.split('*')[0]
				if len(to_write) > 16:
					out_fasta.write(format_header('%s-%s%s\n' % (second_header.rstrip(), str(i), '-'), 'juncN', None, None))
					out_fasta.write(to_write+'\n')
	f.close()
	out_fasta.close()

### These functions are used for fusions

def pull_fusion_files(fusion_dir):
	'''Finds all *.txt files in the directory.'''
	files = os.listdir(fusion_dir)
	fusion_files = []
	for f in files:
		if f.endswith('.txt'):
			fusion_files.append(f)
	return fusion_files

def merge_and_filter_fusion_files(fusion_dir, result_dir):
	'''Merges fusion files found in fusion_dir, filtering out anything that doesn't match a threshold for Junction Read Count. Returns 1 if unable to open any junction files or nothing if it runs to completion.'''
	
	num_observ_map = {}
	jrc_threshold = 1
	
	# Tries to pull a list of junction files. Returns with a warning if unable to find any.
	if not os.path.isdir(fusion_dir):
		# No junction files are going to be found. Gotta leave the function.
		warnings.warn("Unable to open fusion folder %s. Giving up on finding gene fusions." % junc_dir)
		return 1
	fusion_files = pull_fusion_files(fusion_dir)
	if len(fusion_files) == 0:
		warnings.warn("Unable to find any .txt files in %s." % fusion_dir)
		return 1
	
	# Write a header line
	w = open(result_dir+'/log/merged-fusions.txt','w')
	f = open(fusion_dir+'/'+fusion_files[0],'r')
	w.write(f.readline())
	f.close()
	
	count_all = 0
	count_thresholded = 0
	
	for file in fusion_files:
		f = open(fusion_dir+'/'+file,'r')
		line = f.readline() # skip headers
		line = f.readline()
		while line:
			count_all += 1
			spline = line.split()
			if int(spline[4]) >= jrc_threshold:
				w.write(line)
				count_thresholded += 1
			line = f.readline()
		f.close()
	
	write_to_log("\n---Fusion Filtering Results---",result_dir+'/log.txt')
	write_to_log("Found %d fusions for this sample" % count_all,result_dir+'/log.txt')
	write_to_log("%d fusions passed threshold of %d JunctionReadCount" % (count_thresholded, jrc_threshold), result_dir+'/log.txt')

def create_fusion_bed(result_dir):
	'''Translate the fusion file to a DNA bed file. One entry for each side, since we can't guarantee they're on the same chromosome'''

	f = open(result_dir+'/log/merged-fusions.txt','r')
	#chr1	67000041	67000091	NP_001337147	0	+	67000041	67000091	0	1	50,	0,
	w = open(result_dir+'/log/fusions.bed','w')
	
	header = f.readline()
	line = f.readline()
	while line:
		spline = line.split()
		# FASTA header: fusion name, left break, right break, JRC
		fus_name = '%s|leftBreak=%s|rightBreak=%s|junctionReadCount=%s' % (spline[0], spline[1], spline[2], spline[4])
		# Grab 50 bases before the left breakpoint
		first_chr = spline[1].split(':')[0]
		first_chr_pos = int(spline[1].split(':')[1])
		first_chr_strand = spline[1].split(':')[2]
		if first_chr_strand == '+':
			first_chr_start = first_chr_pos-50
		else:
			first_chr_pos -= 1
			first_chr_start = first_chr_pos+50
		w.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t1\t50\t0\n" % (first_chr,min(first_chr_start, first_chr_pos),max(first_chr_start, first_chr_pos),fus_name,first_chr_strand,min(first_chr_start, first_chr_pos),max(first_chr_start, first_chr_pos)))
		# ...and 50 after the right
		second_chr = spline[2].split(':')[0]
		second_chr_pos = int(spline[2].split(':')[1])
		second_chr_strand = spline[2].split(':')[2]
		if second_chr_strand == '+':
			second_chr_pos -= 1
			second_chr_end = second_chr_pos+50
		else:
			second_chr_end = second_chr_pos-50
		w.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t1\t50\t0\n" % (second_chr,min(second_chr_end, second_chr_pos),max(second_chr_end, second_chr_pos),fus_name,second_chr_strand,min(second_chr_end, second_chr_pos),max(second_chr_end, second_chr_pos)))
		
		line = f.readline()
		
	f.close()
	w.close()

def translate_fusions(result_dir):
	'''In the bed file we read in, the sequences come in pairs - first sequence in the pair is the left gene, second is the right gene.'''
	
	translate_table = maketrans("ACGTacgt","TGCAtgca")
	
	f = open(result_dir+'/log/fusions.bed.dna','r')
	w = open(result_dir+'/log/fusions.fasta','w')
	# We want to take the 4th line and every 5th line after that (that's the line with the 50 bases we want to keep, the rest is various buffers and metadata stuff)
	pos = 2
	second = False
	seq = ""
	line = f.readline()
	while line:
		if pos%5 == 0:
			if line.split()[5] == '-':
				seq += line.rstrip().split()[-1][::-1].translate(translate_table)
			else:
				seq += line.rstrip().split()[-1]
			if second:
				seq = seq.upper()
				for i in range(0,3):
					w.write('>%s|%d+\n' % (line.split()[0],i))
					w.write('%s\n' % translate_seq(seq[i:], '+', True))
					w.write('>%s|%d-\n' % (line.split()[0],i))
					w.write('%s\n' % translate_seq(seq[:len(seq)-i], '-', True))
				second = False
				seq = ""
			else:
				second = True
		line = f.readline()
		pos += 1

### This function is used at the end to combine all of the output fasta files into one.

def combine_output_fastas(out_dir):
	'''Combines all of the fastas in the output directory into one'''
	files = os.listdir(out_dir+'parts/')
	for f in files:
		if f.endswith('.fasta'):
			with open(out_dir+'parts/'+f,'r') as src:
				with open(out_dir+"variant_proteome.fasta", 'a') as dest:
					shutil.copyfileobj(src, dest)

### These functions are used everywhere.

def write_to_log(message, log_file):
	'''Writes a message to the output log.'''
	# I think this is more efficient than opening and writing to the end of the file with the 
	# Python I/O tools, but if not it would probably be more convenient to use those.
	# Fixed it so we can write messages with returns in them!
	msg = message.split('\n')
	for m in msg:
		call("echo "+m+" >> "+log_file, shell=True)

def write_to_status(message):
	'''Writes a message to the status log.'''
	# Same comment as in write_to_log
	msg = message.split('\n')
	for m in msg:
		call("echo "+m+" >> "+statusfile, shell=True)

def raise_warning(warn_message):
	'''The default warning is ugly! I'm making a better one. Okay I'm not, this is a waste of time right now.'''
	my_warning = warnings.warn(warn_message)
	print my_warning
	
def quit_if_no_variant_files(args):
	if not (args.germline or args.somatic or args.junction or args.fusion):
		raise SystemExit("ERROR: Couldn't find any variant files!\nAborting program.")

# Main function!
if __name__ == "__main__":
	# Parse input, make sure we have at least one variant file.
	args = parse_input_arguments()
	script_dir = os.path.dirname(os.path.realpath(__file__)) # can this really be the best way to do this!?

	# Set up log/status files
	output_dir = args.output_dir
	results_folder = set_up_output_dir(output_dir, args)
	write_to_status("Started")
	write_to_log("Version Python.0", logfile)
	write_to_log("Reference DB used: "+args.proteome.split("/")[-1], logfile)
	
	# Set up codon map
	codon_map = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
	
	# Time to merge and quality-threshold the variant files!
	if args.somatic:
		som_flag = merge_and_qual_filter(args.somatic, args.variant_quality_threshold)
		if som_flag:
			args.somatic = None
	if args.germline:
		germ_flag = merge_and_qual_filter(args.germline, args.variant_quality_threshold)
		if germ_flag:
			args.germline = None
	quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
	write_to_status("Merge and qual filter finished")
		
	# Now let's remove everything in the somatic file that is duplicated in the germline file
	# since if it shows up in both, it's a germline variant.
	if args.somatic and args.germline:
		remove_somatic_duplicates(args.germline, args.somatic)
	write_to_status("Somatic duplicates removed")

	# Call read_chr_bed.c, which takes the reference genome and proteome.bed file as input and produces
	# a fasta file of exomes.
	# Possible but unlikely future work: rewrite the C file (still in C though) so it's more efficient?
	# I dunno, it seems fine for now.
	try:
		check_call("%s/read_chr_bed %s/log/proteome.bed %s" % (script_dir, results_folder, args.genome), shell=True)
	except CalledProcessError:
		raise SystemExit("ERROR: read_chr_bed didn't work - now we don't have a proteome.bed.dna file.\nAborting program.")
	# Commented the above out for speed - it's slow, so for current testing purposes I'm just copying it from elsewhere
	#shutil.copy('/ifs/data/proteomics/tcga/scripts/quilts/pyquilts/proteome.bed.dna', results_folder+"/log/")
	
	# Next, create a proteome.bed file containing only variants...probably.
	# I could probably combine them more prettily, but for now I'll just concatenate the files.
	if args.somatic:
		get_variants(args.somatic+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "S")
	if args.germline:
		get_variants(args.germline+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "G")
	write_to_status("Get variants completed, proteome.bed file written")	
	
	# Combine them, not very prettily:
	if args.somatic and args.germline:
		dest = open(results_folder+"/log/proteome.bed.var",'w')
		for f in [results_folder+"/log/proteome.bed.S.var",results_folder+"/log/proteome.bed.G.var"]:
			with open(f,'r') as src:
				shutil.copyfileobj(src, dest)
		dest.close()
	elif args.somatic:
		shutil.copy(results_folder+"/log/proteome.bed.S.var", results_folder+"/log/proteome.bed.var")
	elif args.germline:
		shutil.copy(results_folder+"/log/proteome.bed.G.var", results_folder+"/log/proteome.bed.var")
	write_to_status("Somatic and germline combined, proteome.bed.SG or whatever output")

	# Finish out the variants, if there are any	
	if args.somatic or args.germline:
		# Combine (some more) and sort variants.
		# Also variants are sorted by what they do to the sequence (add/remove stop, change AA, etc)
		# I will continue to do this, probably? I'll just ignore more later on
		sort_variants(results_folder+"/log/proteome.bed.dna", results_folder+"/log/proteome.bed.var", "proteome")
		write_to_status("Variants sorted")	

		# Translate the variant sequences into a fasta file.
		translate(results_folder+"/log/", "proteome.aa.var.bed.dna", logfile, 'aa')
		write_to_status("Translated SAAVs")	
		translate(results_folder+"/log/", "proteome.indel.var.bed.dna", logfile, 'indel')
		write_to_status("Translated indels")
		shutil.copy(results_folder+"/log/proteome.aa.var.bed.dna.fasta", results_folder+"/fasta/parts/proteome.aa.fasta")
		shutil.copy(results_folder+"/log/proteome.indel.var.bed.dna.fasta", results_folder+"/fasta/parts/proteome.indel.fasta")

		# Time for tryptic peptides! Remember: C-term (after) K/R residues
		#make_aa_peptide_fasta(results_folder+"/log/", args.no_missed_cleavage)
		#write_to_status("Tryptic peptides done.")
		#make_indel_peptide_fasta(results_folder+"/log/", logfile)
		#write_to_status("Indel tryptic peptides done.")
	
	# Let's move on to alternate splice junctions.
	# We're removing some of the earlier parts for now to expedite testing.
	if args.junction:
		write_to_status("Starting on junctions now.")
		junc_flag = None
		# If MapSplice was used instead of Tophat, copy its junctions.txt file to 
		# junctions.bed so the merge_junction_files function will pick up on it
		if not args.mapsplice:
			junc_flag = convert_tophat_to_mapsplice(args.junction)
			if junc_flag:
				args.junction = None
				write_to_log("Could not convert junctions.bed to junctions.txt - check junction bed file location", logfile)
				warnings.warn("Could not convert %s/junctions.bed to junctions.txt - check junction bed file location. Skipping..." % args.junction)
				quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
		# Merge junction files found in junction folder.
		if not junc_flag:
			junc_flag = merge_junction_files(args.junction, results_folder+'/log')
		if junc_flag:
			args.junction = None
			quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
		else:
			filter_known_transcripts(args.proteome+'/transcriptome.bed', results_folder+'/log', logfile)
			filter_alternative_splices(results_folder+'/log/', args.threshA, args.threshAN, args.threshN, logfile)
			write_to_status("Filtered alternative splices into the appropriate types.")
		
			# Make a fasta out of the alternative splices with conserved exon boundaries
			write_to_status("About to do a read_chr_bed")
			try:
				check_call("%s/read_chr_bed %s/log/merged-junctions.alt.filtered.bed %s" % (script_dir, results_folder, args.genome), shell=True)
				# Don't know why this copies instead of moving. If I never use merged-junctions.filter.A.bed.dna again, just move it or have read_chr_bed output the alternative.bed.dna file instead.
				shutil.copy(results_folder+'/log/merged-junctions.alt.filtered.bed.dna', results_folder+'/log/alternative.bed.dna')
			except CalledProcessError:
				warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.alt.filtered.bed.dna file. Will not have a fasta file of alternative splices with conserved exon boundaries.")
			write_to_status("Done with read_chr_bed to create alternative.bed.dna")
			
			# Make a fasta out of the alternative splices with conserved donor boundaries
			write_to_status("About to do a read_chr_bed")
			try:
				check_call("%s/read_chr_bed %s/log/merged-junctions.donor.filtered.bed %s" % (script_dir, results_folder, args.genome), shell=True)
				# Don't know why this copies instead of moving.
				shutil.copy(results_folder+'/log/merged-junctions.donor.filtered.bed.dna', results_folder+'/log/donor.bed.dna')
			except CalledProcessError:
				warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.donor.filtered.bed.dna file. Will not have a fasta file of alternative splices with conserved exon boundaries.")
			write_to_status("Done with read_chr_bed to create donor.bed.dna")
			
			# Now to tackle the novels...
			write_to_status("About to do a read_chr_bed")
			try:
				check_call("%s/read_chr_bed %s/log/merged-junctions.novel.filtered.bed %s" % (script_dir, results_folder, args.genome), shell=True)
				shutil.copy(results_folder+'/log/merged-junctions.novel.filtered.bed.dna', results_folder+'/log/novel.bed.dna')
			except CalledProcessError:
				warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.novel.filtered.bed.dna file. Will not have a fasta file of novel spliceforms.")
			write_to_status("Done with read_chr_bed to create novel.bed.dna")
		
			# Translating the junctions. Looks like it requires a slightly different function than the old translation function.
			# Basically, though, can use the indel translation function (keep the exon with the variant and everything that comes after) for all of them except the ones where the exon boundaries are both new, in which case we need to do all six reading frames. And we saved those...where? notA?
			# Single frame translations
			translate(results_folder+"/log/", "alternative.bed.dna", logfile, 'juncA')
			translate(results_folder+"/log/", "donor.bed.dna", logfile, 'juncAN')
			# Six-frame translations
			#shutil.copy('/ifs/data/proteomics/tcga/scripts/quilts/pyquilts/results_20171116.144632/log/novel_splices.bed.dna', results_folder+'/log/novel_splices.bed.dna')
			translate_novels(results_folder+"/log/", "novel.bed.dna", logfile)
		
			# Move junctions to results folder
			#shutil.copy(results_folder+"/log/alternative.bed.dna.fasta", results_folder+"/fasta/parts/proteome.alternative.fasta")
			#shutil.copy(results_folder+"/log/alternative_frameshift.bed.dna.fasta", results_folder+"/fasta/parts/proteome.alternative_frameshift.fasta")
			shutil.copy(results_folder+"/log/alternative.bed.dna.fasta", results_folder+"/fasta/parts/proteome.alternative_splices.fasta")
			with open(results_folder+"/log/donor.bed.dna.fasta", 'r') as src:
				with open(results_folder+"/fasta/parts/proteome.alternative_splices.fasta",'a') as dest:
					shutil.copyfileobj(src, dest)
			shutil.copy(results_folder+"/log/novel.bed.dna.fasta", results_folder+"/fasta/parts/proteome.novel_splices.fasta")
			write_to_status("Translated alternative-splice junctions")
	
	if args.fusion:
		write_to_status("Starting on fusions now.")
		fusion_flag = merge_and_filter_fusion_files(args.fusion, results_folder)
		if fusion_flag:
			args.fusion = None
			quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
		else:
			#fusion time
			create_fusion_bed(results_folder)
			write_to_status("About to do a read_chr_bed")
			try:
				check_call("%s/read_chr_bed %s/log/fusions.bed %s" % (script_dir, results_folder, args.genome), shell=True)
			except CalledProcessError:
				warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a fusions.bed.dna file. Will not be able to perform fusions.")
			write_to_status("Done with read_chr_bed to fusions.bed.dna")
			translate_fusions(results_folder)
			shutil.copy(results_folder+"/log/fusions.fasta", results_folder+"/fasta/parts/proteome.fusions.fasta")
		
	# Combine all of the proteome parts into one.
	combine_output_fastas(results_folder+'/fasta/')
	write_to_status("DONE")
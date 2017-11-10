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
	parser.add_argument('--output_dir', type=str, default="/ifs/data/proteomics/tcga/scripts/quilts/pyquilts",
		help="full path to output folder")
	parser.add_argument('--proteome', type=str, default="/ifs/data/proteomics/tcga/databases/refseq_human_20160914", help="full path to folder containing reference proteome")
	parser.add_argument('--genome', type=str, default="/ifs/data/proteomics/tcga/databases/genome_human", help="full path to folder containing reference genome")
	#parser.add_argument('--somatic', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna/vcf/TCGA-20130502-S", help="VCF file of somatic variants")
	#parser.add_argument('--germline', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna-germline/vcf/GATK26-G-Nature", help="VCF file of germline variants")
	#parser.add_argument('--junction', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/rna/tophat/junction_tophat208bowtie2_gMx", help="BED file of splice junctions [in progress]")
	parser.add_argument('--germline',type=str)
	parser.add_argument('--somatic',type=str)
	parser.add_argument('--junction', type=str)
	parser.add_argument('--mapsplice', action='store_true', default=False, help="Use this if your junctions were created by MapSplice rather than Tophat. Automatically converts the junctions.txt file at the junction location to a file called junctions.bed which can be found by this script.")
	parser.add_argument('--fusion', type=str, help="Fusion genes [currently unsupported]")
	parser.add_argument('--threshA', type=int, default=2)
	parser.add_argument('--threshAN', type=int, default=3)
	parser.add_argument('--threshN', type=int, default=3)
	parser.add_argument('--variant_quality_threshold', type=float, default=0.0, help="Quality threshold for variants")
	parser.add_argument('--no_missed_cleavage', action='store_true', default=False, help="Tryptic peptide fasta by default allows for a single missed cleavage; adding this argument will tell the virtual trypsinizer to assume perfect cleavage")
	
	# Pull out the arguments
	args = parser.parse_args()
	# Check that we have a somatic and/or germline and/or junction file. Abort if not.
	# Will add fusion to this check later.
	if not args.somatic and not args.germline and not args.junction:
		raise SystemExit("ERROR: Must have at least one variant file (somatic, germline, and/or junctions).\nAborting program.")
	
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
		raise SystemExit("ERROR: Output directory does not exist.\nAborting program.")
	
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
		shutil.copy(args.proteome+"/proteome.fasta",results_folder+"/fasta/parts/"+args.proteome.split("/")[-1]+".fasta")
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

	# Set up some tracking stats etc.
	qual_removed_all = 0
	total_variants_all = 0
	kept_variants_all = 0
	duplicates_all = 0
	header_written = False
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
			if len(spline) < 6:
				# Bummer, we don't have six tab-separated fields, this isn't a valid VCF line
				write_to_log("Error parsing: "+line, vcf_log_location)
				line = f.readline()
				continue
			if spline[4] == '.':
				# It's a "no variant". Ignore it.
				# This hasn't been tested too rigorously because it didn't appear in any of my test data,
				# but someone else had some of these so I had him insert this and it seemed to work for him.
				#write_to_log("No variant: "+line, vcf_log_location) # Uncomment this if you want to write these lines to the log for testing
				line = f.readline()
				continue
			if line[0] == '#':
				# It's a header line. Write it as-is to the file if it's the first file.
				if not header_written:
					w.write('\t'.join(spline[0:6])+'\n')
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
		header_written = True # Don't want to write a header unless it's the first file
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
			chr, start, name, lengths, offsets = spline[0].lstrip('chr'), int(spline[1]), spline[3], spline[-2], spline[-1]
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
		return changes, aa_substs
	
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
				if chunk_pos+len(to_delete) > chunk_length:
					write_to_status("Found a deletion that goes beyond the exon. %s\t%s" % (ex_head[0], indel))
				seq_chunk = orig_seq_chunk[:chunk_pos]+to_insert+orig_seq_chunk[chunk_pos+len(to_delete):]
				if chunk_pos+len(to_delete) <= chunk_length:
					new_length = chunk_length-len(to_delete)+len(to_insert)
				else:
					new_length = chunk_pos+1 # adds one for the base we're adding back in
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

def translate_seq(sequence, strand):
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
		translated += codon_map[sequence[i:i+3]]
	
	# Now this should return the translated sequence up to the first stop codon.	
	# That could either be at the actual end of the sequence, or at a variant that adds a stop.
	return translated.split('*')[0]

def translate(log_dir, bed_file):
	'''Reads in a bed.dna file and translates it to a protein .fasta file.'''
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
			try:
				if prev_gene != gene:
					prev_AA_subst = []
				if AA_subst not in prev_AA_subst:
					translated = translate_seq(sequence.upper(), strand)
					out_fasta.write(second_header)
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
	try:
		if prev_gene != gene or AA_subst != prev_AA_subst:
			translated = translate_seq(sequence.upper(), strand)
			out_fasta.write(second_header)
			out_fasta.write(translated+'\n')
	except UnboundLocalError:
		pass

	f.close()
	out_fasta.close()

### These functions are used to make the peptide fasta files.

def trypsinize(sequence, pos = 0):
	'''Virtually trypsinizes a protein and returns its tryptic peptides. 
	If a position is given, it only returns the peptides following that position (used for indels
	and splice junctions, where anything before the variant is unchanged).
	Right now, just chops it after any K or R that isn't followed by a P.'''
	tryptic_peptides = []
	seq = ''
	for i in range(len(sequence)):
		letter = sequence[i]
		seq += letter
		try:
			next_letter = sequence[i+1]
		except IndexError:
			next_letter = '*'
		if (letter == 'K' or letter == 'R') and next_letter != 'P':
			if pos <= i:
				tryptic_peptides.append(seq)
			seq = ''
	if seq != '':
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
	'''Finds all *junctions.bed files in the directory.'''
	files = os.listdir(junc_dir)
	junc_files = []
	for f in files:
		if f.endswith('junctions.bed'):
			junc_files.append(f)
	return junc_files

def merge_junction_files(junc_dir, log_dir):
	'''Merges junction files found in junc_dir. Returns 1 if unable to open any junction files or nothing if it runs to completion.'''
	
	# May or may not allow this to change. Not a thousand percent sure what it is.
	min_length = 8
	min_sum_length = 45
	
	num_observ_map = {}
	junction_map = {} # Is never used.
	
	# Tries to pull a list of junction files. Returns with a warning if unable to find any.
	if not os.path.isdir(junc_dir):
		# No junction files are going to be found. Gotta leave the function.
		warnings.warn("Unable to open junction folder %s. Giving up on finding splice junctions." % junc_dir)
		return 1
	junc_files = pull_junc_files(junc_dir)
	if len(junc_files) == 0:
		warnings.warn("Unable to find any *junctions.bed files in %s." % junc_dir)
		return 1
	
	for file in junc_files:
		label = ''
		if '-' in file:
			label = '-'+file.split('-')[0]
			
		f = open(junc_dir+'/'+file,'r')
		for line in f.readlines():
			spline = line.split()
			if len(spline) > 11:
				chrm, begin_ex_1, end_ex_2, junc_num, num_observ, len_exons, begin_exons = spline[0], int(spline[1]), int(spline[2]), spline[3], int(spline[4]), spline[10], spline[11]
				splen_exons = map(int,len_exons.split(','))
				spbegin_exons = map(int,begin_exons.split(','))
				begin_intron = begin_ex_1 + splen_exons[0] + 1
				end_intron = begin_ex_1 + spbegin_exons[-1] -1
				if 'chr' not in chrm:
					chrm = 'chr'+chrm
				if splen_exons[0] >= min_length and splen_exons[1] >= min_length and sum(splen_exons) >= min_sum_length:
					key = '%s#%d#%d' % (chrm,begin_intron,end_intron)
					# Why are these two (adding an int and appending to a list) so different? It's hideous. But this is as efficient as I can get with it right now.
					num_observ_map[key] = num_observ_map.get(key,0) + int(num_observ)
					junction_map.setdefault(key,[]).append(junc_num+label)
		f.close()
		
	# Time to write out our results
	w = open(log_dir+'/merged-junctions.bed','w')
	count = 0
	for key in sorted(num_observ_map.keys()):
		spkey = key.split('#')
		chr, beg_int, end_int = spkey[0], int(spkey[1]), int(spkey[2])
		begin_ex_1 = beg_int-50-1
		end_ex_2 = end_int+50+1
		start_ex_2 = end_int+1-begin_ex_1 
		w.write("%s\t%d\t%d\tj%d\t%d\t+\t%d\t%d\t0\t2\t50,50\t0,%d\n" % (chr,begin_ex_1,end_ex_2,count,num_observ_map[key],begin_ex_1,end_ex_2,start_ex_2))
		count += 1
	w.close()
							
def filter_known_transcripts(transcriptome_bed, results_folder, logfile):
 	'''Filters out spliceform transcripts that are known from those that aren't, and annotates the ones that can be annotated. Haven't done any serious testing on this yet, so...yeah, use at yr own risk'''

	junctions_model = {}
	count_model_exons = 0
	count_model_begin = 0
	count_model_end = 0
	
	try:
		f = open(transcriptome_bed,'r')
	except IOError:
		raise_warning("Could not find transcriptome file at %s!\n Treating all junctions as novel." % transcriptome_bed)
	
	if f:	
		line_number = 0 # Might never get used
		for line in f.readlines():
			spline = line.split('\t')
			chrm, begin_exon, protein, strand, block_count, block_sizes, block_starts = spline[0], int(spline[1]), spline[3], spline[5], int(spline[9]), spline[10], spline[11]
			sizes = map(int, block_sizes.split(','))
			starts = map(int, block_starts.split(','))
			start = begin_exon

			for i in range(0,block_count-1):
				begin_intron = start+sizes[i]+1
				end_intron = begin_exon+starts[i+1]-1
				key = "%s#%d#%d" % (chrm, begin_intron, end_intron)
				junctions_model[key] = "C"
				if i==0:
					key = "%s#%d" % (chrm, begin_intron)
					if strand == '+':
						junctions_model[key] = "B"
						count_model_begin += 1
					else:
						junctions_model[key] = "E"
						count_model_end += 1
				if i==block_count-1-1:
					key = "%s#%d" % (chrm, end_intron)
					if strand == '-':
						junctions_model[key] = "B"
						count_model_begin += 1
					else:
						junctions_model[key] = "E"
						count_model_end += 1
				start = begin_exon + starts[i+1]
				count_model_exons += 1
			line_number += 1
		f.close()
	
	f = open(results_folder+'/merged-junctions.bed','r')
	w = open(results_folder+'/merged-junctions.filter.bed','w')
	num_observ_max = 0 # number of maximum observations for a particular intron in RNAseq data
	count_junct_map = {'B': 0, 'E': 0, 'C': 0, '': 0}
	count_junctions = 0
	stat = {}
	
	for line in f.readlines():
		spline = line.split()
		if len(spline) > 11:
			chrm, begin_ex_1, end_ex_2, junc_num, num_observ, len_exons, begin_exons = spline[0], int(spline[1]), int(spline[2]), spline[3], int(spline[4]), spline[10], spline[11]
			splen_exons = map(int,len_exons.split(','))
			spbegin_exons = map(int,begin_exons.split(','))
			non_matches = 0
			junction_status = ''
			num_observ_max = max(num_observ, num_observ_max)
			if num_observ >= 1:	
				begin_intron = begin_ex_1 + splen_exons[0]+1
				end_intron = begin_ex_1 + spbegin_exons[-1]-1
				length_intron = end_intron-begin_intron
				# Figure out whether it already exists in the database and what type it is
				if "%s#%d#%d" % (chrm, begin_intron, end_intron) in junctions_model.keys():
					junction_status = 'C'
				elif "%s#%d" % (chrm, begin_intron) in junctions_model.keys():
					if junction_status == '':
						junction_status = junctions_model["%s#%d" % (chrm, begin_intron)]
				else:
					if "%s#%d" % (chrm, end_intron) in junctions_model.keys():
						if junction_status == '':
							junction_status = junctions_model["%s#%d" % (chrm, end_intron)]
				# Log it
				if junction_status == '':
					w.write(line)
				else:
					key = "%s#%d" % (junction_status, num_observ)
					stat[key] = stat.get(key,0) + 1
					count_junct_map[junction_status] += 1
				count_junctions += 1
	f.close()
	w.close()
	
	# Write out a stats file 
	w = open(results_folder+'/merged-junctions.filter.stats.txt','w')
	w.write("num_observ\tC_matches\tB_matches\tE_matches\n")
	for i in range(1,num_observ_max+1):
		w.write(str(i))
		for type in ['C','B','E','']:
			key = '%s#%d' % (type, i)
			if key not in stat.keys():
				stat[key] = 0
			w.write('\t%d' % (stat[key]))
		w.write('\n')
	w.close()
	
	write_to_log("\n---Transcript Filtering Results---",logfile)
	write_to_log("Model exons: %d" % count_model_exons, logfile)
	write_to_log("Model transcript begin: %d" % count_model_begin, logfile)
	write_to_log("Model transcript end: %d" % count_model_end, logfile)
	write_to_log("Total RNA-Seq introns: %d" % count_junctions, logfile)
	for type in ['C','B','E','']:
		write_to_log("RNA-Seq introns \"%s\": %d" % (type, count_junct_map[type]), logfile)

def filter_alternative_splices(log_dir, threshA, threshAN, threshN, logfile):
	# Leaving out "print details", which seems to just write out a file for each protein, which seems...unnecessary and excessive? If it needs to come back it can come back.
	# And if I do put that back, I need to make sure to find all the points in the original code where it's used (search print_details==1, OUT_DETAILS, should pick them up)
	
	### Keeping track of a bunch of things from the proteome.bed file
	f = open(log_dir+'/proteome.bed','r')
	
	models = {} # Stores tuple from proteome.bed for each protein
	junctions_model = {}
	junctions_model_A = {}
	junctions_model_begin_intron = {}
	junctions_model_end_intron = {}
	count_model_exons = 0
	
	line_number = 0
	line = f.readline()
	while line:
		spline = line.rstrip().split()
		chr, begin_exon, end_exon, protein, strand, block_count, block_sizes, block_starts = spline[0], int(spline[1]), int(spline[2]), spline[3], spline[5], int(spline[9]), spline[10], spline[11]
		models[protein] = line
		sizes = map(int,block_sizes.split(','))
		starts = map(int,block_starts.split(','))
		start = begin_exon
		for i in range(0,block_count-1):
			begin_intron = start + sizes[i] + 1
			end_intron = begin_exon + starts[i+1] - 1
			junctions_model.setdefault("%s#%d#%d" % (chr, begin_intron, end_intron),[]).append("%s,%d,%d" % (protein, i, i+1))
			junctions_model_begin_intron.setdefault("%s#%d" % (chr, begin_intron),[]).append("%s,%d" % (protein, i))
			junctions_model_end_intron.setdefault("%s#%d" % (chr, end_intron),[]).append("%s,%d" % (protein, i+1))
			
			# Looks to pick alternate blocks for possible weird splicing? is this really the most efficient way to do this? maybe it will make more sense later.
			for j in range(i+2,block_count):
				end_intron_tmp = begin_exon + starts[j] - 1
				junctions_model_A.setdefault("%s#%d#%d" % (chr, begin_intron, end_intron_tmp),[]).append("%s,%d,%d" % (protein, i, j))
			start = begin_exon + starts[i+1]
		count_model_exons += block_count # does this belong here? what is this even tracking
		line_number += 1
		line = f.readline()
	
	f.close()

	### Time to match our putative splice junctions to model introns.
	### ...and write out a whole ton of files. Are any of these useful?
	### I'll rename them when I find out what they are. Stop naming things with underscores!
	f = open(log_dir+'merged-junctions.filter.bed','r')
	outfile = open(log_dir+'merged-junctions.filter.A.bed','w')
	outfile_ = open(log_dir+'merged-junctions.filter.A_.bed','w')
	outfile_AN = open(log_dir+'merged-junctions.filter.AN.bed','w')
	outfile_AN_ = open(log_dir+'merged-junctions.filter.AN_.bed','w')
	outfile_NOT = open(log_dir+'merged-junctions.filter.notA.bed','w')
	outfile_EXTRA = open(log_dir+'merged-junctions.filter.A.extra.bed','w')
	outfile_EXTRA_ = open(log_dir+'merged-junctions.filter.A.extra_.bed','w')
	outfile_AN_EXTRA = open(log_dir+'merged-junctions.filter.AN.extra.bed','w')
	outfile_AN_EXTRA_ = open(log_dir+'merged-junctions.filter.AN.extra_.bed','w')
	outfile_NOT_EXTRA = open(log_dir+'merged-junctions.filter.notA.extra.bed','w')
	
	num_observ_max = 0
	junction_count = 0
	count_junctions_map = {}
	stat_map = {}
	
	line = f.readline()
	while line:
		spline = line.rstrip().split()
		chr, begin_ex_1, end_ex_2, junc_num, num_observ, len_exons, begin_exons = spline[0], int(spline[1]), int(spline[2]), spline[3], int(spline[4]), map(int,spline[10].split(',')), map(int,spline[11].split(','))
		size = [] # Will this get used? I bet it won't
		non_matches = 0
		junc_status = ""
		N = True # This is an indicator for something.
		num_observ_max = max(num_observ_max, num_observ)
		if num_observ >= 1: # should always be
			begin_intron = begin_ex_1 + len_exons[0] + 1
			end_intron = begin_ex_1 + begin_exons[-1] - 1
			stop_ex_1 = begin_intron - 1
			start_ex_2 = end_intron + 1
			len_intron = end_intron - begin_intron
			
			key = "%s#%d#%d" % (chr, begin_intron, end_intron)
			# This apparently never executes. But why?
			if key in junctions_model.keys():
				print "Actually it totally does execute, dorkus %s" % key
				junc_status = "C"
				N = False
			else:
				# I know these upcoming functions have a lot more parameters than should 
				# be strictly necessary, but I thought this function was getting lengthy
				# and I wanted to break it up a bit.
				if key in junctions_model_A.keys():
					# If the RNAseq intron matched any of the beginnings/ends of introns
					# from the gene, new gene created with exons that precede & follow
					# RNAseq intron. (Conserved exon boundaries? I want to say yes?)
					N = False
					junc_status = "A"
					#write_to_log("%s: %s" % (junc_status, key), logfile)
					gn_name = conserved_exon_boundaries(junctions_model_A[key], models, num_observ, start_ex_2, threshA, outfile, outfile_EXTRA, outfile_, outfile_EXTRA_)
				else:
					key = "%s#%d" % (chr, begin_intron)
					if key in junctions_model_begin_intron.keys():
						# If only the beginning of the RNAseq intron was matched to 
						# the existing in-database beginning of the DNA intron, we have
						# a truncated exon.
						N = False
						junc_status = "AN1"
						#write_to_log("%s: %s" % (junc_status, key), logfile)
						gn_name, outside = truncated_exon(junctions_model_begin_intron[key], models, end_intron, num_observ, start_ex_2, threshAN, outfile_AN, outfile_AN_EXTRA, outfile_AN_, outfile_AN_EXTRA_)
						if outside:
							junc_status = ""
					# Should this next chunk be in another else?
					key = "%s#%d" % (chr, end_intron)
					if key in junctions_model_end_intron.keys():
						# If only the end of the RNAseq intron was matched to 
						# the existing in-database end of the DNA intron, we have
						# the elongation of an exon within an intron.
						N = False
						junc_status = "AN2"
						#write_to_log("%s: %s" % (junc_status, key), logfile)
						gn_name, outside = elongated_exon(junctions_model_end_intron[key], models, begin_intron, num_observ, start_ex_2, threshAN, outfile_AN, outfile_AN_EXTRA, outfile_AN_, outfile_AN_EXTRA_)
						if outside:
							junc_status = ""
					if N or gn_name == "NO NAME":
						# There should be another "or" condition here, I guess, where 
						# something is returned from elongated/truncated weirdly?
						# For now I'm just going to keep this one and see what happens.
						# It all comes out the same.
						gn_name = "NO_GENE-N-%d-%s-%s-%d-%d" % (num_observ, chr, junc_num, begin_intron, end_intron)
		
			new_line = "%s\t%d\t%d\t%s\t%s\n" % (chr, begin_ex_1, end_ex_2, gn_name, '\t'.join(spline[4:]))
			if junc_status == "":
				if threshN <= num_observ:
					outfile_NOT.write(new_line)
				else:
					outfile_NOT_EXTRA.write(new_line)
			else:
				key = "%s#%d" % (junc_status, num_observ)
				stat_map[key] = stat_map.get(key,0) + 1
				count_junctions_map[junc_status] = count_junctions_map.get(junc_status, 0) + 1
				if junc_status[:2] == "AN":
					if threshAN <= num_observ:
						outfile_NOT.write(new_line)
					else:
						outfile_NOT_EXTRA.write(new_line)
			junction_count += 1
					
		line = f.readline()
	
	f.close()
	outfile.close()
	outfile_.close()
	outfile_AN.close()
	outfile_AN_.close()
	outfile_NOT.close()
	outfile_EXTRA.close()
	outfile_EXTRA_.close()
	outfile_AN_EXTRA.close()
	outfile_AN_EXTRA_.close()
	outfile_NOT_EXTRA.close()
	
	w = open(log_dir+"merged-junctions.filterA.stat.txt",'w')
	w.write('num_observ\tC_matches\tA_matches\tAN1_matches\tAN2_matches\n')
	for i in range(1,num_observ_max+1):
		w.write(str(i))
		for type in ["C","A","AN1","AN2"]:
			key = "%s#%d" % (type, i)
			w.write("\t%d" % stat_map.get(key, 0))
		w.write("\n")
	w.close()
	
	write_to_log("Model exons: %d\n" % count_model_exons, logfile)
	write_to_log("RNAseq introns: %d\n" % junction_count, logfile)
	for type in ["C","A","AN1","AN2"]:
		write_to_log("%s: %d\n" % (type, count_junctions_map.get(type, 0)), logfile)

def conserved_exon_boundaries(prot_info, models, num_observ, start_exon_2, threshA, out_file, out_extra, out_file_, out_extra_):
	# If I have a variable with an underscore at the end, that means it's just something I've formatted to be printed later and isn't used in any calculations or anything.
	for prot in prot_info:
		sprot = prot.split(',')
		protein, alt1, alt2 = sprot[0],int(sprot[1]),int(sprot[2])
		if protein in models:
			split_prot = models[protein].rstrip().split()
			chr, start, end, name, block_count, block_sizes, block_starts = split_prot[0], int(split_prot[1]), int(split_prot[2]), split_prot[3], int(split_prot[9]), map(int,split_prot[10].split(',')), map(int,split_prot[11].split(','))
			start_exon_1 = ""
			size_exon_1 = ""
			size_exon_2 = ""
			removed_length = 0
			
			block_count_ = 0
			block_sizes_ = []
			block_starts_ = []
			for i in range(0, block_count):
				if i <= alt1 or alt2 <= i:
					block_count_ += 1
					block_sizes_.append(str(block_sizes[i]))
					block_starts_.append(str(block_starts[i]))
				else:
					removed_length += block_sizes[i]
				if i == alt1:
					start_exon_1 = start+block_starts[i] # start location of exon 1 that bridges
					size_exon_1 = block_sizes[i]
				if i == alt2:
					# Do we use this?
					#stop_exon_2=start+block_starts[i]+block_sizes[i] #the end location of the exon2 that bridges
					size_exon_2 = block_sizes[i]
			# Is "chr" the same as "chromosome"? It's just the name, should be fine
			name_ = "%s-A-%d-%d-%d-%s-%d-%d-%d" % (name, alt1+1, alt2+1, num_observ, chr, start_exon_1+size_exon_1, start_exon_2, removed_length)
			to_print = "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t0,0,255\t%d\t%s\t%s\n" % (chr, start, end, name_, num_observ, split_prot[5], split_prot[6], split_prot[7], block_count_, ','.join(block_sizes_), ','.join(block_starts_))
			if removed_length % 3 == 0:
				if threshA <= num_observ:
					out_file.write(to_print)
				else:
					out_extra.write(to_print)
			else:
				if threshA <= num_observ:
					out_file_.write(to_print)
				else:
					out_extra_.write(to_print)
	try:
		return name_
	except UnboundLocalError:
		return "NO NAME"

def truncated_exon(prot_info, models, end_intron, num_observ, start_exon_2, threshAN, out_an, out_an_extra, out_an_, out_an_extra_):
	# If I have a variable with an underscore at the end, that means it's just something I've formatted to be printed later and isn't used in any calculations or anything.
	outside = True
	for prot in prot_info:
		sprot = prot.split(',')
		protein, alt1 = sprot[0], int(sprot[1])
		if protein in models.keys():
			split_prot = models[protein].rstrip().split()
			chr, start, end, name, strand, block_count, block_sizes, block_starts = split_prot[0], int(split_prot[1]), int(split_prot[2]), split_prot[3], split_prot[5], int(split_prot[9]), map(int,split_prot[10].split(',')), map(int,split_prot[11].split(','))
			start_exon_1 = ""
			size_exon_1 = ""
			size_exon_2 = ""
			
			if end_intron >= end:
				# Not within the same gene
				# I think nothing happens if this occurs.
				# Move to next protein.
				print "OH NO %d %d" % (end_intron, end)
				continue
			outside = False
			block_count_ = 0
			block_sizes_ = []
			block_starts_ = []
			modify = True
			for i in range(0, block_count):
				# Create new gene from exons that precede the beginning of the RNAseq intron
				# and the exons that come after the end of the RNAseq intron
				if i == alt1: 
					# We're at the right start intron
					start_exon_1 = start+block_starts[i] # start location of exon 1
					size_exon_1 = block_sizes[i]
				if i <= alt1 or start+block_starts[i-1]+block_sizes[i-1] > end_intron:
					# Add sizes and starts on both sides of the junction
					# The above is guaranteed to happen if this happens, is that cool?
					block_count_ += 1
					block_sizes_.append(str(block_sizes[i]))
					block_starts_.append(str(block_starts[i]))
				else:
					# I don't really know what's going on here, honestly
					if start+block_starts[i]+block_sizes[i] < end_intron:
						# If end of RNAseq intron between ends of [i]th & [i+1]th exons,
						# modifies [i+1] exon. We...don't know why? At least KVR doesn't
						# and I'll have to think about it.
						e1 = start+block_starts[i]
						e2 = start+block_starts[i]+block_sizes[i]
					else:
						if modify:
							i1 = start+block_starts[i]+block_sizes[i]-end_intron # new size
							i2 = end_intron-start+1 # new start
							size_exon_2 = i1
							block_count_ += 1
							block_sizes_.append(str(i1))
							block_starts_.append(str(i2))
							modify = False
				

			# Is "chr" the same as "chromosome"? Think so. It's just the name, should be fine
			name_ = "%s-AN1-%d-%d-%s-%d-%d" % (name, alt1+1, num_observ, chr, start_exon_1+size_exon_1, start_exon_2)
			to_print_beginning = "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t" % (chr, start, end, name_, num_observ, split_prot[5], split_prot[6], split_prot[7])
			to_print_end = "\t%d\t%s\t%s\n" % (block_count_, ','.join(block_sizes_), ','.join(block_starts_))
			if modify:
				# At some point move this to write_to_log if it seems important enough to keep
				raise_warning("Error: %s\n" % name)
			else:
				if strand == '+':
					to_print = "%s0,100,0%s" % (to_print_beginning, to_print_end)
					if threshAN <= num_observ:
						out_an.write(to_print)
					else:
						out_an_extra.write(to_print)
				else:
					to_print = "%s107,142,35%s" % (to_print_beginning, to_print_end)
					if threshAN <= num_observ:
						out_an_.write(to_print)
					else:
						out_an_extra_.write(to_print)
		#else:
		#	print "%s not in models, format should be %s" % (protein, models.keys()[0])
	try:
		return name_, outside
	except UnboundLocalError:
		return "NO NAME", outside

def elongated_exon(prot_info, models, begin_intron, num_observ, start_exon_2, threshAN, out_an, out_an_extra, out_an_, out_an_extra_):
	# If I have a variable with an underscore at the end, that means it's just something I've formatted to be printed later and isn't used in any calculations or anything.
	outside = True
	for prot in prot_info:
		sprot = prot.split(',')
		protein, alt2 = sprot[0], int(sprot[1])
		if protein in models.keys():
			split_prot = models[protein].rstrip().split()
			chr, start, end, name, strand, block_count, block_sizes, block_starts = split_prot[0], int(split_prot[1]), int(split_prot[2]), split_prot[3], split_prot[5], int(split_prot[9]), map(int,split_prot[10].split(',')), map(int,split_prot[11].split(','))
			start_exon_1 = ""
			size_exon_1 = ""
			size_exon_2 = ""
			
			if start >= begin_intron:
				# Not within the same gene
				# I think nothing happens if this occurs.
				# Move to next protein.
				#outside = False
				print "OH NO %s %d %d %d" % (protein, alt2, start, begin_intron)
				continue
			outside = False
			block_count_ = 0
			block_sizes_ = []
			block_starts_ = []
			modify = True
			offset = 0
			for i in range(0, block_count):
				# Create new gene from exons that precede the beginning of the RNAseq intron
				# and the exons that come after the end of the RNAseq intron
				if i == alt2: 
					# We're at the right start intron
					size_exon_2 = block_sizes[i] # end of exon 2
				if i == 0 and begin_intron < start:
					# This should never happen, right? since this originally was
					# inside an if block for start<begin_intron
					# If beginning of intron is before beginning of gene
					block_count_ += 1
					e = block_sizes[i] + 201 # why 201?
					offset = start-begin_intron+e
					block_sizes_.append(str(e))
					block_starts_.append(str(block_starts[i]))
				else:
					if i >= alt2 or start+block_starts[i+1]<begin_intron:
						# Create new gene from exons preceding start of intron and exons
						# following end of intron
						block_count_ += 1
						block_sizes_.append(str(block_sizes[i]))
						e = block_starts[i] + offset
						block_starts_.append(str(e))
					else:
						if start+block_starts[i] > begin_intron:
							# If beginning of RNAseq intron between beginnings
							# of [i]th & [i+1]th exons, modifies [i]th exon.
							e1 = start + block_starts[i]
							e2 = start + block_starts[i] + block_sizes[i]
						else:
							# What is the point of this? Why not just leave this stuff
							# where it was?
							if modify:
								i1 = begin_intron - (start + block_starts[i])
								i2 = block_starts[i] # exon start
								size_exon_1 = i1
								block_count_ += 1
								block_sizes_.append(str(i1))
								block_starts_.append(str(i2))
								modify = False
								
			# Is "chr" the same as "chromosome"? It's just the name, should be fine
			#name_ = "%s-AN2-%d-%d-%s-%d-%d" % (name, alt2+1, num_observ, chr, int(start_exon_1)+int(size_exon_1), start_exon_2)
			name_ = "%s-AN2-%d-%d-%s-%d-%d" % (name, alt2+1, num_observ, chr, begin_intron, start_exon_2)
			to_print_beginning = "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t" % (chr, start-offset, end, name_, num_observ, split_prot[5], split_prot[6], split_prot[7])
			to_print_end = "\t%d\t%s\t%s\n" % (block_count_, ','.join(block_sizes_), ','.join(block_starts_))
			if modify:
				# At some point move this to write_to_log if it seems important enough to keep
				raise_warning("Error: %s\n" % name)
			else:
				if strand == '-':
					to_print = "%s184,134,11%s" % (to_print_beginning, to_print_end)
					if threshAN <= num_observ:
						out_an.write(to_print)
					else:
						out_an_extra.write(to_print)
				else:
					to_print = "%s255,215,0%s" % (to_print_beginning, to_print_end)
					if threshAN <= num_observ:
						out_an_.write(to_print)
					else:
						out_an_extra_.write(to_print)
		#else:
		#	print "%s not in models, format should be %s" % (protein, models.keys()[0])	
	try:
		return name_, outside
	except UnboundLocalError:
		return "NO NAME", outside

### These functions are used to translate junction files.


def translate_junctions(input_file):
	'''Main function for adding non-frameshifted junctions to the tryptic peptide fasta.
	Basically, add the peptide that includes the junction and everything that comes after it.
	Later: maybe remove the parts that come after it if it isn't frameshifted? That might be complicated.
	NOT FINISHED'''
	# Grab the indel variants
	vars = {}
	f = open(input_file,'r')
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
		translate(results_folder+"/log/", "proteome.aa.var.bed.dna")
		write_to_status("Translated SAAVs")	
		translate(results_folder+"/log/", "proteome.indel.var.bed.dna")
		write_to_status("Translated indels")

		# Time for tryptic peptides! Remember: C-term (after) K/R residues
		make_aa_peptide_fasta(results_folder+"/log/", args.no_missed_cleavage)
		write_to_status("Tryptic peptides done.")
		make_indel_peptide_fasta(results_folder+"/log/", logfile)
		write_to_status("Indel tryptic peptides done.")
	
	# Let's move on to alternate splice junctions.
	# We're removing some of the earlier parts for now to expedite testing.
	if args.junction:
		write_to_status("Starting on junctions now.")
		# If MapSplice was used instead of Tophat, copy its junctions.txt file to 
		# junctions.bed so the merge_junction_files function will pick up on it
		if args.mapsplice:
			try:
				shutil.copy(args.junction+'/junctions.txt', args.junction+'/junctions.bed')
			except IOError:
				write_to_log("Could not copy junctions.txt to junctions.bed - check MapSplice junction file location", logfile)
				warnings.warn("Could not copy %s/junctions.txt to junctions.bed - check MapSplice junction file location. Skipping..." % args.junction)
		# Merge junction files found in junction folder.
		junc_flag = merge_junction_files(args.junction, results_folder+'/log')
		if junc_flag:
			args.junction = None
		quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
		filter_known_transcripts(args.proteome+'/transcriptome.bed', results_folder+'/log', logfile)
		shutil.copy('/ifs/data/proteomics/tcga/scripts/quilts/pyquilts/merged-junctions.filter.bed', results_folder+"/log/")
		filter_alternative_splices(results_folder+'/log/', args.threshA, args.threshAN, args.threshN, logfile)
		write_to_status("Filtered alternative splices into the appropriate types.")
		# Make a fasta out of the alternative splices with conserved exon boundaries
		write_to_status("About to do a read_chr_bed")
		try:
			check_call("%s/read_chr_bed %s/log/merged-junctions.filter.A.bed %s" % (script_dir, results_folder, args.genome), shell=True)
			# Don't know why this copies instead of moving. If I never use merged-junctions.filter.A.bed.dna again, just move it or have read_chr_bed output the alternative.bed.dna file instead.
			shutil.copy(results_folder+'/log/merged-junctions.filter.A.bed.dna', results_folder+'/log/alternative.bed.dna')
		except CalledProcessError:
			warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.filter.A.bed.dna file. Will not have a fasta file of alternative splices with conserved exon boundaries.")
		write_to_status("Done with read_chr_bed to create alternative.bed.dna")
				
		'''
		# Grabbing the variants that are in these splice junctions, I guess. This seems like fishing for dregs.
		# Removing this for now, just because it takes so long and I might want to do it differently.
		if args.somatic:
			get_variants(args.somatic+"/merged_pytest/merged.vcf", results_folder+"/log/merged-junctions.filter.A.bed", "S")
			if not args.germline: # Move the somatic variants to the merged .bed.var file, otherwise wait
				shutil.copy(results_folder+"/log/merged-junctions.filter.A.bed.S.var", results_folder+"/log/merged-junctions.filter.A.bed.var")
		if args.germline:
			get_variants(args.germline+"/merged_pytest/merged.vcf", results_folder+"/log/merged-junctions.filter.A.bed", "G")
			if args.somatic: # Put the germline variants in with the somatic ones
				with open(results_folder+"/log/merged-junctions.filter.A.bed.G.var", 'r') as src:
					with open(results_folder+"/log/merged-junctions.filter.A.bed.var",'w') as dest:
						shutil.copyfileobj(src, dest)
				dest.close()
				src.close()
			else:
				shutil.copy(results_folder+"/log/merged-junctions.filter.A.bed.var", results_folder+"/log/merged-junctions.filter.A.bed.G.var")
		if args.somatic or args.germline:
			sort_variants(results_folder+"/log/merged-junctions.filter.A.bed.dna", results_folder+"/log/merged-junctions.filter.A.bed.var", "merged-junctions.filter.A")
			# Move them into the file with the other junctions for translating
			with open(results_folder+"/log/merged-junctions.filter.A.aa.var.bed.dna", 'r') as src:
					with open(results_folder+"/log/alternative.bed.dna",'a') as dest:
						shutil.copyfileobj(src, dest)
			dest.close()
			src.close()
		'''
		
		# Okay, well, I'm going to try to just do the same thing on A_, AN, and AN_.
		write_to_status("About to do a read_chr_bed")
		shutil.copy(results_folder+"/log/merged-junctions.filter.A_.bed", results_folder+"/log/merged-junctions.filter.A_AN.bed")
		with open(results_folder+"/log/merged-junctions.filter.AN.bed", 'r') as src:
			with open(results_folder+"/log/merged-junctions.filter.A_AN.bed",'a') as dest:
				shutil.copyfileobj(src, dest) # Adds in AN

		# Does the read_chr_bed		
		try:
			check_call("%s/read_chr_bed %s/log/merged-junctions.filter.A_AN.bed %s" % (script_dir, results_folder, args.genome), shell=True)
			# Don't know why this copies instead of moving. If I never use merged-junctions.filter.A.bed.dna again, just move it or have read_chr_bed output the alternative.bed.dna file instead.
			shutil.copy(results_folder+'/log/merged-junctions.filter.A_AN.bed.dna', results_folder+'/log/alternative_frameshift.bed.dna')
		except CalledProcessError:
			warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.filter.A_AN.bed.dna file. Will not have a fasta file of alternative splices with conserved exon boundaries.")
		
		# Translating the junctions. Looks like it requires a slightly different function than the old translation function.
		# Basically, though, can use the indel translation function (keep the exon with the variant and everything that comes after) for all of them except the ones where the exon boundaries are both new, in which case we need to do all six reading frames. And we saved those...where? notA?
		# Single frame translations
		translate(results_folder+"/log/", "alternative.bed.dna")
		translate(results_folder+"/log/", "alternative_frameshift.bed.dna")
		# Six-frame translations
		#translate_six_frame()
		write_to_status("Translated alternative-splice junctions")
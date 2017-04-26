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

Am I ever going to get rid of duplicate tryptic peptides (peptides that show up in multiple proteins)? Maybe that can be a postprocessing step. Or maybe that's a search problem.

Emily Kawaler
'''

import argparse
import os
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
	parser.add_argument('--proteome', type=str, default="/ifs/data/proteomics/tcga/databases/refseq_human_20130727", help="full path to folder containing reference proteome")
	parser.add_argument('--genome', type=str, default="/ifs/data/proteomics/tcga/databases/genome_human", help="full path to folder containing reference genome")
	parser.add_argument('--somatic', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna/vcf/TCGA-20130502-S", help="VCF file of somatic variants")
	parser.add_argument('--germline', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna-germline/vcf/GATK26-G-Nature", help="VCF file of germline variants")
	#parser.add_argument('--germline', type=str, help="VCF file of germline variants")
	parser.add_argument('--junction', type=str, help="BED file of splice junctions [currently unsupported]")
	parser.add_argument('--fusion', type=str, help="Fusion genes [currently unsupported]")
	parser.add_argument('--threshA', type=int, default=2)
	parser.add_argument('--threshAN', type=int, default=3)
	parser.add_argument('--threshN', type=int, default=3)
	parser.add_argument('--variant_quality_threshold', type=float, default=0.0, help="Quality threshold for variants")
	parser.add_argument('--no_missed_cleavage', action='store_true', default=False, help="Tryptic peptide fasta by default allows for a single missed cleavage; adding this argument will tell the virtual trypsinizer to assume perfect cleavage")
	
	# Pull out the arguments
	args = parser.parse_args()
	# Check that we have a somatic and/or germline file. Abort if not.
	# Will add junction/fusion to this check later.
	if not args.somatic and not args.germline:
		raise SystemExit("ERROR: Must have at least one variant file (somatic and/or germline).\nAborting program.")
	
	return args

def set_up_output_dir(output_dir, ref_proteome):
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
		shutil.copy(ref_proteome+"/proteome.fasta",results_folder+"/fasta/parts/"+ref_proteome.split("/")[-1]+".fasta")
	except IOError:
		raise SystemExit("ERROR: Reference proteome .fasta file not found at "+ref_proteome+"/proteome.fasta.\nAborting program.")	
	try:
		shutil.copy(ref_proteome+"/proteome.bed",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome .bed file not found at "+ref_proteome+"/proteome.bed.\nAborting program.")
	# Not sure if these two should be quite so important.
	try:
		shutil.copy(ref_proteome+"/proteome-descriptions.txt",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome gene descriptions file not found at "+ref_proteome+"/proteome-descriptions.txt.\nAborting program.")
	try:
		shutil.copy(ref_proteome+"/proteome-genes.txt",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome gene names file not found at "+ref_proteome+"/proteome-genes.txt.\nAborting program.")
	
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
			chr, start, name, lengths, offsets = spline[0].lstrip('chr'), int(spline[1]), spline[3], spline[10], spline[11]
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
			# Maybe I shouldn't subtract one from pos here - did it when I made the vcf file in the first place
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
		if len(orig_nt) != 1 or len(new_nt) != 1:
			continue
		triplet_start = ((pos/3)*3) # start of the triplet containing pos. counting from 0, not 1
		triplet_orig = full_seq[triplet_start:(triplet_start+3)].upper()
		triplet_subst = var.split(':')[0].split(str(pos))[-1]
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
	return changes, aa_substs

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
		second_header[-1] = second_header_map_end[:-2]+' '+var+')'
		out_aa_bed_dna.write(' '.join(second_header)+(' (VAR:%s-%s%d%s:%s)\n' % (var[0], subst[0], subst[1], subst[2], var.split(':')[1]))) # Will write the AA position starting from 0, not 1. This is what old QUILTS did. We can change it.
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
		
def sort_variants(proteome_file, variant_file):
	'''Goes through the variants and sorts them by type, writing out a bunch of intermediate files in the process. Right now the types are "single-AA non-stop variant" and "not that", but eventually will have more.'''
	global codon_map
	
	# Setting this up here because this is the first time we'll use it. But it won't be the last.
	codon_map = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
	
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
	out_aa = open(file_base+"/proteome.bed.aa.var", 'w')
	# This one will basically just be that first line. Might want an extra line for each variant though, depends how it's used. Actually, almost certainly do.
	out_aa_bed = open(file_base+"/proteome.aa.var.bed", 'w')
	# proteome.bed.dna, but with the variants included. Def want an extra line for each variant.
	out_aa_bed_dna = open(file_base+"/proteome.aa.var.bed.dna", 'w')
	# Right now just proteome.bed.aa.var but with the non-AA variants. Doesn't do anything right now because I don't care about it very much.
	out_other = open(file_base+"/proteome.bed.other.var", 'w')
	
	f = open(proteome_file, 'r')
	line = f.readline()
	while line:
		# Line type one: First gene header line
		if line[:3] == 'chr':
			# Finish processing the previous gene - focus first on proteome.bed.aa.var
			# Using a try/except block here to deal with the possibility of this being the first line
			# There has to be a more graceful way to do that, right?
			try:
				changed_vars, aa_substs = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []))
				if changed_vars != []:
					write_out_aa(name, changed_vars, aa_substs, out_aa)
					write_out_aa_bed(header_line.split('\t'), changed_vars, out_aa_bed)
					write_out_aa_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, changed_vars, aa_substs, out_aa_bed_dna)
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
		changed_vars, aa_substs = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []))
		if changed_vars != []:
			write_out_aa(name, changed_vars, aa_substs, out_aa)
			write_out_aa_bed(header_line.split('\t'), changed_vars, out_aa_bed)
			write_out_aa_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, changed_vars, aa_substs, out_aa_bed_dna)
	except UnboundLocalError:
		pass

	# Close everything
	f.close()
	out_aa.close()
	out_aa_bed.close()
	out_aa_bed_dna.close()
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

def translate(log_dir):
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
	f = open(log_dir+"proteome.aa.var.bed.dna",'r')
	out_fasta = open(log_dir+"proteome.aa.var.bed.dna.fasta",'w') # The original QUILTS writes two other files but they're just duplicates of proteome.aa.var.bed and proteome.aa.var.bed.dna, it seems. For now I'm leaving them out.
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

def trypsinize(sequence):
	'''Virtually trypsinizes a protein and returns its tryptic peptides. Right now, just chops it after any K or R that isn't followed by a P.'''
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
		if 25 >= len(new_pep) >= 6:
			out_file.write('>%s START:%d END:%d VAR:%s\n' % (gene, peptide_start_pos, peptide_start_pos+len(new_pep)-1, ','.join(var_string)))
			out_file.write('%s\n' % new_pep)
		new_peptides.append([new_pep, var_string])
	return new_peptides

def write_missed_cleavage_peptides(gene, prev_peptides, cur_peptides, cur_start_pos, out_file):
	'''Allows for a single missed cleavage - basically does a Cartesian join between all variants of the previous peptide and all variants of the current peptide (including the case where one, but not both, has no variants)'''
	if prev_peptides != [] and cur_peptides != []:
		prev_peptide_start = cur_start_pos-len(prev_peptides[0][0])
		cur_peptide_end = cur_start_pos+len(cur_peptides[0][0])-1
		for i in prev_peptides:
			for j in cur_peptides:
				if i[1] != [] or j[1] != []:
					vars = i[1]+j[1]
					if 25 >= len(i[0]+j[0]) >= 6:
						out_file.write('>%s START:%d END:%d (missed cleavage after %d) VAR:%s\n' % (gene, prev_peptide_start, cur_peptide_end, cur_start_pos-1, ','.join(vars)))
						out_file.write('%s\n' % (i[0]+j[0]))

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
	for i in range(len(tryptic_peptides)):
		peptide = tryptic_peptides[i]
		vars_in_peptide = [var for var in vars.keys() if total_aas <= var-1 < total_aas+len(peptide)]
		cur_peptides = write_peptides(gene, {v: vars[v] for v in vars_in_peptide}, total_aas, peptide, out_file)
		if not no_missed_cleavage:
			# If we are allowing missed cleavages...
			write_missed_cleavage_peptides(gene, prev_peptides, cur_peptides, total_aas, out_file)
		prev_peptides = cur_peptides
		total_aas += len(peptide)

def make_peptide_fasta(log_dir, no_missed_cleavage):
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
					assign_variants(gene, tryptic_peptides, vars[gene], tryp_fasta, no_missed_cleavage)
			except UnboundLocalError:
				pass
			# Start processing the new gene
			sequence = ''
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
		line = f.readline()
	# And finish the last one
	try:
		if gene in vars:
			translated = translate_seq(sequence.upper(), strand) # no variants here yet
			tryptic_peptides = trypsinize(translated)
			assign_variants(gene, tryptic_peptides, vars[gene], tryp_fasta, no_missed_cleavage)
	except UnboundLocalError:
		pass
		
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
	
# Main function!
if __name__ == "__main__":
	# Parse input, make sure we have at least one variant file.
	args = parse_input_arguments()
	script_dir = os.path.dirname(os.path.realpath(__file__)) # can this really be the best way to do this!?

	# Set up log/status files
	output_dir = args.output_dir
	results_folder = set_up_output_dir(output_dir, args.proteome)
	write_to_status("Started")
	write_to_log("Version Python.0", logfile)
	write_to_log("Reference DB used: "+args.proteome.split("/")[-1], logfile)
	
	# Time to merge and quality-threshold the variant files!
	if args.somatic:
		som_flag = merge_and_qual_filter(args.somatic, args.variant_quality_threshold)
	if args.germline:
		germ_flag = merge_and_qual_filter(args.germline, args.variant_quality_threshold)
	if som_flag and germ_flag:
		# We can keep going with only one variant file, but if we find neither, we have to quit.
		raise SystemExit("ERROR: Couldn't find any variant files!\nAborting program.")
	write_to_status("Merge and qual filter finished")
		
	# Now let's remove everything in the somatic file that is duplicated in the germline file
	# since if it shows up in both, it's a germline variant.
	if args.somatic and args.germline and not (som_flag or germ_flag):
		remove_somatic_duplicates(args.germline, args.somatic)
	write_to_status("Somatic duplicates removed")

	'''# Call read_chr_bed.c, which takes the reference genome and proteome.bed file as input and produces
	# a fasta file of exomes.
	# Possible but unlikely future work: rewrite the C file (still in C though) so it's more efficient?
	# I dunno, it seems fine for now.
	try:
		check_call("%s/read_chr_bed %s/log/proteome.bed %s" % (script_dir, results_folder, args.genome), shell=True)
	except CalledProcessError:
		raise SystemExit("ERROR: read_chr_bed didn't work - now we don't have a proteome.bed.dna file.\nAborting program.")'''
	# Commented the above out for speed - it's slow, so for current testing purposes I'm just copying it from elsewhere
	shutil.copy('/ifs/data/proteomics/tcga/scripts/quilts/pyquilts/proteome.bed.dna', results_folder+"/log/")
	
	# Next, create a proteome.bed file containing only variants...probably.
	# I could probably combine them more prettily, but for now I'll just concatenate the files.
	if args.somatic and not som_flag:
		get_variants(args.somatic+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "S")
	if args.germline and not germ_flag:
		get_variants(args.germline+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "G")
	write_to_status("Get variants completed, proteome.bed file written")	
	
	# Combine them, not very prettily:
	if args.somatic and args.germline and not (som_flag or germ_flag):
		dest = open(results_folder+"/log/proteome.bed.var",'w')
		for f in [results_folder+"/log/proteome.bed.S.var",results_folder+"/log/proteome.bed.G.var"]:
			with open(f,'r') as src:
				shutil.copyfileobj(src, dest)
		dest.close()
	elif args.somatic and not som_flag:
		shutil.copy(results_folder+"/log/proteome.bed.S.var", results_folder+"/log/proteome.bed.var")
	elif args.germline and not germ_flag:
		shutil.copy(results_folder+"/log/proteome.bed.G.var", results_folder+"/log/proteome.bed.var")
	write_to_status("Somatic and germline combined, proteome.bed.SG or whatever output")
	
	# Combine (some more) and sort variants.
	# Also variants are sorted by what they do to the sequence (add/remove stop, change AA, etc)
	# I will continue to do this, probably? I'll just ignore more later on
	sort_variants(results_folder+"/log/proteome.bed.dna", results_folder+"/log/proteome.bed.var")
	write_to_status("Variants sorted")	

	# Translate the variant sequences into a fasta file.
	translate(results_folder+"/log/")
	write_to_status("Translated")	

	# Time for tryptic peptides! Remember: C-term (after) K/R residues
	#make_peptide_fasta(results_folder+"/log/", args.no_missed_cleavage)
	#write_to_status("Tryptic peptides done. Should be finished.")

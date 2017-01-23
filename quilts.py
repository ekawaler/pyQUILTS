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
	# The only one I found that has both somatic and germline
	parser.add_argument('--somatic', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna/vcf/TCGA-20130502-S", help="VCF file of somatic variants")
	parser.add_argument('--germline', type=str, default="/ifs/data/proteomics/tcga/samples/breast/TCGA-E2-A15A/dna-germline/vcf/GATK26-G-Nature", help="VCF file of germline variants")
	parser.add_argument('--junction', type=str, help="BED file of splice junctions [currently unsupported]")
	parser.add_argument('--fusion', type=str, help="Fusion genes [currently unsupported]")
	parser.add_argument('--threshA', type=int, default=2)
	parser.add_argument('--threshAN', type=int, default=3)
	parser.add_argument('--threshN', type=int, default=3)
	parser.add_argument('--variant_quality_threshold', type=float, default=0.0, help="Quality threshold for variants")
	
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
			try:
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
	global codon_map
	
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
		triplet_new = triplet_orig[:subst_pos] + triplet_subst + triplet_orig[subst_pos+1:] # This may change it in both.
		# If it's the reverse strand, check the reversed and translated triplet instead
		# I guess the .vcf files are based on proteome.bed.dna? So the position and nucleotides refer
		# to the negative strand if proteome.bed.dna only has the negative strand.
		if reverse_flag:
			triplet_orig = triplet_orig[::-1].translate(translate_table)
			triplet_new = triplet_new[::-1].translate(translate_table)
		AA_old = codon_map[triplet_orig]
		AA_new = codon_map[triplet_new]
		# If neither the old nor the new codon is a stop codon and the AA changes,
		# count it!
		if AA_old != AA_new and AA_new != '*' and AA_old != '*':
			changes.append(var)
	return changes

def sort_variants(proteome_file, variant_file):
	# sort_variants.pl proteome.bed.dna proteome.bed.var
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
	'''w = open(variant_file+".SG.combined", 'w')
	for var in variants_map:
		w.write("%s\t%s\n" % (var, ','.join(variants_map[var])))
	w.close()'''
	
	# Open some files to write to. Right now, only looking at single AA changes.
	# I hate these filenames. Do they actually mean anything?
	out_aa = open(file_base+"/proteome.bed.aa.var", 'w')
	# This one will basically just be that first line
	out_aa_bed = open(file_base+"/proteome.aa.var.bed", 'w')
	out_aa_bed_dna = open(file_base+"/proteome.aa.var.bed.dna", 'w')
	out_other = open(file_base+"/proteome.bed.other.var", 'w')
	
	f = open(proteome_file, 'r')
	line = f.readline()
	while line:
		# Line type one: First gene header line
		if line[:3] == 'chr':
			# Finish processing the previous gene - focus first on proteome.bed.aa.var
			# Using a try/except block here to deal with the possibility of this being the first line
			# There has to be a more graceful way to do this, right?
			try:
				changed_vars = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []))
				if changed_vars != []:
					out_aa.write("%s\t%s\t%s\n" % (header_line.split('\t')[3], ','.join(changed_vars), strand))
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
		# I think we can ignore the pre- and post-100 for now, since they're mostly used for stop codon stuff
		# Eventually, though, we'll need them.
		line = f.readline()
		
	# Do the last one
	try:
		changed_vars = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []))
		if changed_vars != []:
			out_aa.write("%s\t%s\n" % (header_line.split('\t')[3], ','.join(changed_vars)))
	except UnboundLocalError:
		pass

	# Close everything
	f.close()
	out_aa.close()
	out_aa_bed.close()
	out_aa_bed_dna.close()
	out_other.close()
	
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
		
	# Now let's remove everything in the somatic file that is duplicated in the germline file
	# since if it shows up in both, it's a germline variant.
	if args.somatic and args.germline and not (som_flag or germ_flag):
		remove_somatic_duplicates(args.germline, args.somatic)

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
	# perl $script_root/get_variants.pl "$result_dir/log/proteome.bed" $somatic/merged/merged.vcf S
	# I could probably combine them more prettily, but for now I'll just concatenate the files.
	if args.somatic and not som_flag:
		get_variants(args.somatic+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "S")
	if args.germline and not germ_flag:
		get_variants(args.germline+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "G")
		
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
		
	# Combine (some more) and sort variants. (sort_variants.pl proteome.bed.dna proteome.bed.var)
	# It appears that this is the point in the original when overlapping G/S were removed (around line 203)
	# Also variants are sorted by what they do to the sequence (add/remove stop, change AA, etc)
	# I will continue to do this, probably? I'll just ignore more later on
	sort_variants(results_folder+"/log/proteome.bed.dna", results_folder+"/log/proteome.bed.var")
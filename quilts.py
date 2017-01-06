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

# ahhhh look at these hideous global variables
global logfile
global statusfile
global results_folder

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
	parser.add_argument('--variant_quality_threshold', type=float, default=15.0, help="Quality threshold for variants")
	
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

			# If we passed those checks, let's say it's good...
			total_variants += 1
			# Not sure why we subtract one from pos, but it happened in the original.
			try:
				chr, pos, id, old, new, qual = spline[0].lstrip('chr'), int(spline[1])-1, spline[2], spline[3], spline[4], float(spline[5])
			except ValueError:
				write_to_log("Error parsing: "+line.rstrip(), vcf_log_location)
				line = f.readline()
				continue	
					
			# Quality check!
			if qual < quality_threshold:
				qual_removed += 1
				line = f.readline()
				continue

			# Is there a higher-quality version of this variant already found?
			line_map_key = "%s#%d#%s#%s" % (chr, pos, old, new)
			if line_map_key in existing_variants:
				if existing_variants[line_map_key][5] > qual:
					duplicates_removed += 1
					line = f.readline()
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

### These functions are used to create bed files of only variants that exist in exons.

def get_variants(vcf_file, proteome_file, type):
	# perl $script_root/get_variants.pl "$result_dir/log/proteome.bed" $somatic/merged/merged.vcf S
	w = open(proteome_file+".var", 'w')
	f = open(proteome_file, 'r')
	
	# Basically goes through proteome.bed and finds all positions in exons
	# Then goes through the .vcf file and checks if they're in any exons
	# I think he saves literally every number. I'm going to just have it check if it's between the start and end, should be more efficient
	# Output line: NP_XX (name) \t type-old(position in genome)new:(quality),(repeat)\t type-old(position in exome)new:(quality),(repeat) Looks like instead of quality it might have taken the wrong field at some point? Sometimes it looks like instead of quality we just get a dot, even when it shouldn't be missing. I'll leave quality in for now and see if we use it later. Output all genes, even if they have no variants.
	
	start_map = {}
	lengths_map = {}
	offsets_map = {}
	pos_map = {}
	est = ExonSearchTree()
	line = f.readline() # Going to assume no headers for now
	
	# I kind of want to make a class for variants but that would probably be overly aggressive. So I won't.
	# But let the record show that I wanted to.
	while line:
		spline = line.rstrip().split('\t')
		try:
			chr, start, name, lengths, offsets = spline[0].lstrip('chr'), int(spline[1]), spline[3], spline[10], spline[11]
		except ValueError:
			# Maybe write this to a log somewhere?
			warnings.warn("Failed to parse %s" % line)
			line = f.readline()
			continue
		start_map[name] = start
		lengths_map[name] = lengths
		offsets_map[name] = offsets
		splengths = lengths.split(',')
		spoffsets = offsets.split(',')
		for i in range(len(spoffsets)):
			# The original saved a map with every position that ever appears in an exon
			# and the names of the genes it goes with. This feels slow and inefficient, so instead
			# I made some wacked-out tree thing for storage and search. See if you like it!
			est.add_exon(chr, int(spoffsets[i])+start, int(spoffsets[i])+start+int(splengths[i]), name)
			
		line = f.readline()

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
	if args.somatic and not som_flag:
		get_variants(args.somatic+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "S")
	#if args.germline and not germ_flag:
		#get_variants(args.germline+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "G")
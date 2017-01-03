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
from subprocess import call
import warnings

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
	day_time_string = str(today.year)+str(today.month)+str(today.day)+'.'+str(today.hour)+str(today.minute)+str(today.second)
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
		raise SystemExit("ERROR: Reference proteome not found at "+ref_proteome+"/proteome.fasta.\nAborting program.")

### These functions are all dicked. Will be deleted gradually once I've scavenged them for anything useful.

def check_for_vcf_file(vcf_dir, type):
	''' Hunt down the VCF file. Raise a warning if not found. Helper function for set_up_vcf_single and set_up_vcf_both.'''
	# Right now, only accepts files with <type> (somatic or germline) in the name. can change that if necessary
	files = os.listdir(vcf_dir)
	vcf_filename = None
	for f in files:
		if f.endswith('.vcf') and f.find(type) != -1:
			vcf_filename = f
	if not vcf_filename:
		warnings.warn("Warning: Couldn't find variant file in %s!" % vcf_dir)
	return vcf_filename
		
def set_up_vcf_dir(vcf_dir):
	'''Make the directory we'll put our VCF files in.'''
	if not os.path.isdir(vcf_dir):
		# For now, this is the only thing we care about. If there aren't any VCF files, there's no program.
		# So we abort. Later, we can just raise a non-kill exception and write a complaint to the log.
		raise SystemExit("VCF directory not found at "+vcf_dir+".\nAborting program.")
	try:
		# Change to merged later, but for now...
		os.makedirs(vcf_dir+"/merged_pytest")
	except OSError:
		#raise SystemExit("VCF directory "+vcf_dir+"/merged already exists!\nAborting program.")
		#raise_warning("VCF directory "+vcf_dir+"/merged already exists!")
		warnings.warn("VCF directory "+vcf_dir+"/merged already exists!")

def set_up_vcf_single(quality_threshold, vcf_dir, type):
	'''If we have either somatic or germline but not both, all we do here is the quality threshold and subtraction of one from the position.'''
	
	# Hunt down the VCF file. Abort if not found.
	# (It'll be a warning later when we have fusions/junctions, but for now we need variants.)
	# Right now, only accepts files with "somatic" or "germline" in the name. can change that if necessary
	vcf_filename = check_for_vcf_file(vcf_dir, type)
	if not vcf_filename:
		raise SystemExit("ERROR: Couldn't find variant file at %s!\nAborting program." % vcf_dir)	
	
	# Cool, we have a valid file, let's open our log, output and input files.
	vcf_log_location = vcf_dir+'/merged_pytest/merged.log'
	vcf_log = open(vcf_log_location,'w')
	write_to_log(vcf_dir+'/'+vcf_filename, vcf_log_location)
	f = open(vcf_dir+'/'+vcf_filename,'r')
	w = open(vcf_dir+'/merged_pytest/merged.vcf','w')
	
	# I read one line at a time instead of going with f.readlines() because
	# if the file is long, it's bad to hold all of the lines in memory.
	# In case you were wondering, which you almost certainly weren't.
	line = f.readline()
	qual_removed = 0
	total_variants = 0
	kept_variants = 0
	while line:
		spline = line.split('\t')
		if len(spline) != 6:
			# Bummer, we don't have six tab-separated fields, this isn't a valid VCF line
			write_to_log("Error parsing: "+line, vcf_log_location)
			line = f.readline()
			continue
		if line[0] == '#':
			# It's a header line. Write it as-is to the file.
			w.write(line)
			line = f.readline()
			continue
		
		# If we passed those checks, let's say it's good...
		total_variants += 1
		# Not sure why we subtract one from pos, but it happened in the original.
		try:
			chr, pos, id, old, new, qual = spline[0], int(spline[1])-1, spline[2], spline[3], spline[4], float(spline[5])
		except ValueError:
			write_to_log("Error parsing: "+line.rstrip(), vcf_log_location)
			line = f.readline()
			continue	
					
		# Quality check!
		if qual < quality_threshold:
			qual_removed += 1
			line = f.readline()
			continue
		
		# Passed the quality check? Cool, write it out to the file and continue.
		kept_variants += 1
		w.write(line)
		line = f.readline()
	
	# Write a summary to the log and close all files. We're done here.
	write_to_log("",vcf_log_location)
	write_to_log("%d total variants, %d failed quality check at threshold %f, %d variants kept in final version" % (total_variants, qual_removed, quality_threshold, kept_variants), vcf_log_location)
	
	f.close()
	w.close()
	vcf_log.close()
	
def set_up_vcf_both(quality_threshold, somatic_dir, germline_dir):
	'''Sets up the variant file when we have both somatic and germline. As in single, includes a quality check. Unlike in single, makes sure to exclude any germline variants that already exist in the somatic file.'''

	# We don't care that the VCF file isn't in any kind of order, right?
	# Also: I only save somatic variants to be checked against germline if they pass the quality check. This cool? I think it's cool.
	
	# Hunt down the VCF files. If one isn't found, raise a warning, then go to
	# set_up_vcf_single with the one that is found. If neither are found, return with an error.
	# (It'll be a warning later, when we have junctions/fusions, but for now we need a variant file.)
	somatic_filename = check_for_vcf_file(somatic_dir, "somatic")
	germline_filename = check_for_vcf_file(germline_dir, "germline")
	if not somatic_filename and not germline_filename:
		raise SystemExit("ERROR: Couldn't find either variant file!\nAborting program.")
	elif not somatic_filename:
		warnings.warn("Could not find the somatic variant file. Continuing with germline only.")
		set_up_vcf_single(quality_threshold, germline_dir, "germline")
		return
	elif not germline_filename:
		warnings.warn("Could not find the germline variant file. Continuing with somatic only.")
		set_up_vcf_single(quality_threshold, somatic_dir, "somatic")
		return
	
	# Cool, we have two valid files, let's open our log, somatic and output files.
	vcf_log_location = somatic_dir+'/merged_pytest/merged.log'
	vcf_log = open(vcf_log_location,'w')
	write_to_log("Somatic variant file: "+somatic_dir+'/'+somatic_filename, vcf_log_location)	
	f = open(somatic_dir+'/'+somatic_filename,'r')
	w = open(somatic_dir+'/merged_pytest/merged.vcf','w')
	
	# I copy/pasted code and I'm not proud of myself. Gonna do it again, too.
	# If I end up having to make extensive edits to this code, I'll put it in a function.
	# But I think there are just enough small changes that I need to do it this way :(
	line = f.readline()
	somatic_variants = set([]) # Should be stored as a hash, so quicker lookups than in arrays?
	qual_removed_som = 0
	total_variants_som = 0
	kept_variants_som = 0
	while line:
		spline = line.split('\t')
		if len(spline) != 6:
			# Bummer, we don't have six tab-separated fields, this isn't a valid VCF line
			write_to_log("Error parsing: "+line, vcf_log_location)
			line = f.readline()
			continue
		if line[0] == '#':
			# It's a header line. Write it as-is to the file.
			w.write(line)
			line = f.readline()
			continue
		
		# If we passed those checks, let's say it's good...
		total_variants_som += 1
		# Not sure why we subtract one from pos, but it happened in the original.
		try:
			chr, pos, id, old, new, qual = spline[0], int(spline[1])-1, spline[2], spline[3], spline[4], float(spline[5])
		except ValueError:
			write_to_log("Error parsing: "+line.rstrip(), vcf_log_location)
			line = f.readline()
			continue		

		# Quality check!
		if qual < quality_threshold:
			qual_removed_som += 1
			line = f.readline()
			continue
		
		# Passed the quality check? Cool, write it out to the file, save it to the list and continue.
		kept_variants_som += 1
		somatic_variants.add("%s#%s" % (chr, pos))
		w.write(line)
		line = f.readline()
	
	write_to_log("", vcf_log_location)
	f.close()
	
	# Now do the same thing with germline variants, but in addition to the quality check,
	# remove any that overlap with somatic variants.
	f = open(germline_dir+'/'+germline_filename,'r')
	write_to_log("Germline variant file: "+germline_dir+'/'+germline_filename, vcf_log_location)
	
	line = f.readline()
	qual_removed_germ = 0
	duplicate_removed_germ = 0
	total_variants_germ = 0
	kept_variants_germ = 0
	
	while line:
		spline = line.split('\t')
		if len(spline) != 6:
			# Bummer, we don't have six tab-separated fields, this isn't a valid VCF line
			write_to_log("Error parsing: "+line.rstrip(), vcf_log_location)
			line = f.readline()
			continue
		if line[0] == '#':
			# It's a header line. Skip - should have used the one from the somatic variants..
			line = f.readline()
			continue
		
		# If we passed those checks, let's say it's good...
		total_variants_germ += 1
		# Not sure why we subtract one from pos, but it happened in the original.
		try:
			chr, pos, id, old, new, qual = spline[0], int(spline[1])-1, spline[2], spline[3], spline[4], float(spline[5])
		except ValueError:
			write_to_log("Error parsing: "+line.rstrip(), vcf_log_location)
			line = f.readline()
			continue
		
		# Does it already exist in the somatic file? If so, dump it!
		if chr+"#"+str(pos) in somatic_variants:
			duplicate_removed_germ += 1
			write_to_log("Duplicate found: "+line.rstrip(), vcf_log_location)
			line = f.readline()
			continue
		
		# Quality check!
		if qual < quality_threshold:
			qual_removed_germ += 1
			line = f.readline()
			continue

		# Passed all the checks? Cool, write it out to the file and continue.
		kept_variants_germ += 1
		w.write(line)
		line = f.readline()
					
	# Write a summary to the log and close all files. We're done here.
	write_to_log("",vcf_log_location)
	write_to_log("Somatic stats: %d total variants, %d failed quality check at threshold %f, %d variants kept in final version" % (total_variants_som, qual_removed_som, quality_threshold, kept_variants_som), vcf_log_location)
	write_to_log("Germline stats: %d total variants, %d were duplicates of somatic variants, %d failed quality check at threshold %f, %d variants kept in final version" % (total_variants_germ, duplicate_removed_germ, qual_removed_germ, quality_threshold, kept_variants_germ), vcf_log_location)
	write_to_log("Total stats: %d total variants, %d germline duplicate variants removed, %d failed quality check at threshold %f, %d variants kept in final version" % ((total_variants_germ+total_variants_som), duplicate_removed_germ, (qual_removed_germ+qual_removed_som), quality_threshold, (kept_variants_som+kept_variants_germ)), vcf_log_location)

	f.close()
	w.close()
	vcf_log.close()

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
	
	# Set up log/status files
	output_dir = args.output_dir
	set_up_output_dir(output_dir, args.proteome)
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

'''
	if args.somatic and args.germline:
		set_up_vcf_both(args.variant_quality_threshold, args.somatic, args.germline)
	else:
		if args.somatic:
			set_up_vcf_single(args.variant_quality_threshold, args.somatic, "somatic")
		else:
			set_up_vcf_single(args.variant_quality_threshold, args.germline, "germline")
'''
	
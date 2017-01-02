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
	parser.add_argument('--proteome', type=str, default="/ifs/data/proteomics/tcga/databases/ensembl_human_37.70", help="full path to folder containing reference proteome")
	# The only one I found that has both somatic and germline
	parser.add_argument('--somatic', type=str, default="/ifs/data/proteomics/tcga/samples/breast-xenografts/whim11/rna/vcf", help="VCF file of somatic variants")
	parser.add_argument('--germline', type=str, default="/ifs/data/proteomics/tcga/samples/breast-xenografts/whim11/rna/vcf", help="VCF file of germline variants")
	parser.add_argument('--junction', type=str, help="BED file of splice junctions [currently unsupported]")
	parser.add_argument('--fusion', type=str, help="Fusion genes [currently unsupported]")
	parser.add_argument('--threshA', type=int, default=2)
	parser.add_argument('--threshAN', type=int, default=3)
	parser.add_argument('--threshN', type=int, default=3)
	parser.add_argument('--variant_quality_threshold', type=float, default=1.0, help="Quality threshold for variants")
	
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

### These functions are used while getting variants.

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
	'''If we have either somatic or germline but not both, all we do here is the quality threshold.'''
	
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
		chr, pos, id, old, new, qual = spline[0], int(spline[1])-1, spline[2], spline[3], spline[4], float(spline[5])
		
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
	# Zeroth: Make sure both files are openable.
	# First: Do somatic like we did in set_up_vcf_single, except we now save
	# the chromosomes/positions.
	# Then: Do germline like we did in set_up_vcf_single, except we now don't save
	# the variants that were in the somatic file.
	
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
	write_to_log(somatic_dir+'/'+somatic_filename, vcf_log_location)	
	f = open(somatic_dir+'/'+somatic_filename,'r')
	w = open(somatic_dir+'/merged_pytest/merged.vcf','w')
		
	f.close()
	w.close()
	vcf_log.close()

### These functions are used everywhere.

def write_to_log(message, log_file):
	'''Writes a message to the output log. It's only one line, but I made it its own function
	because it will be much easier to understand in the body of the code.'''
	# Error message can only be one line long, since adding '\n' messes up the call
	# Calling this for each line you want to error-write is a workaround, since it automatically
	# adds one.
	# I think this is more efficient than opening and writing to the end of the file with the 
	# Python I/O tools, but if not it would probably be more convenient to use those.
	call("echo "+message+" >> "+log_file, shell=True)

def write_to_status(message):
	'''Writes a message to the status log. It's only one line, but I made it its own function
	because it will be much easier to understand in the body of the code.'''
	# Same comment as in write_to_log
	call("echo "+message+" >> "+statusfile, shell=True)

def raise_warning(warn_message):
	'''The default warning is ugly! I'm making a better one. Okay I'm not, this is a waste of time right now.'''
	my_warning = warnings.warn(warn_message)
	print my_warning
	
# Main function!
if __name__ == "__main__":
	# Parse input.
	args = parse_input_arguments()
	
	# Set up log/status files
	output_dir = args.output_dir
	set_up_output_dir(output_dir, args.proteome)
	write_to_status("Started")
	write_to_log("Version Python.0", logfile)
	write_to_log("Reference DB used: "+args.proteome.split("/")[-1], logfile)
	
	# Time to get some variants!
	# merge_vcf_orig.pl twice: with somatic first, then germline as argument
	# We are still cool with removing germline variants that also appear in somatic.
	# So...merge_vcf.pl seems like it's really only used to remove repeats. We don't want to do that.
	# If there's only one file, we'll skip this step and just copy the VCF to the working directory.	
	# If there are two, we need to merge them.
	# Kind of clunky right now, will fix when i have a better idea what I'm doing.
	# Seems like the original allows for multiple files. Why? Is this necessary?
	if not args.somatic:
		set_up_vcf_dir(args.germline)
	else:
		set_up_vcf_dir(args.somatic)

	if args.somatic and args.germline:
		set_up_vcf_both(args.variant_quality_threshold, args.somatic, args.germline)
	else:
		if args.somatic:
			set_up_vcf_single(args.variant_quality_threshold, args.somatic, "somatic")
		else:
			set_up_vcf_single(args.variant_quality_threshold, args.germline, "germline")
			pass
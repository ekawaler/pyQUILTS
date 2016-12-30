'''
quilts.py

Remaking QUILTS in Python.

It's in its early stages, so half the lines are just me asking myself questions
that I will forget to answer in the future.

Emily Kawaler
'''

import argparse
import os
import datetime
from subprocess import call

global logfile
global statusfile
global results_folder

def parse_input_arguments():
	# Set up the argument parser
	# Do I want to have some sort of root directory so they don't have to enter full paths? Probably not - 
	# full paths are a pain but they're more flexible.
	parser = argparse.ArgumentParser(description="QUILTS") # What even is this description for, anyway? Maybe I'll remove it eventually
	parser.add_argument('--output_dir', type=str, default="/ifs/data/proteomics/tcga/scripts/quilts/quilts_py/test",
		help="full path to output folder") # Make not-optional at some point
	parser.add_argument('--proteome', type=str, default="ensembl_human_37.70", help="folder containing reference proteome") # Decide whether to keep defaults later. If not, make sure to check them in this function.
	parser.add_argument('--somatic', type=str, help="VCF file of somatic variants")
	parser.add_argument('--germline', type=str, help="VCF file of germline variants")
	parser.add_argument('--junction', type=str, help="BED file of splice junctions [currently unsupported]")
	parser.add_argument('--fusion', type=str, help="Fusion genes [currently unsupported]")
	parser.add_argument('--threshA', type=int, default=2)
	parser.add_argument('--threshAN', type=int, default=3)
	parser.add_argument('--threshN', type=int, default=3)
	
	# Pull out the arguments
	args = parser.parse_args()
	# Check that we have a somatic and/or germline file. Abort if not.
	# Will add junction/fusion to this check later.
	if not args.somatic and not args.germline:
		raise SystemExit("ERROR: Must have at least one variant file (somatic and/or germline).\nAborting program.")
	
	return args

def set_up_output_dir(output_dir):
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
	today = datetime.datetime.today()
	day_time_string = str(today.year)+str(today.month)+str(today.day)+'.'+str(today.hour)+str(today.minute)+str(today.second)
	results_folder = output_dir+'/results_'+day_time_string
	os.makedirs(results_folder)
	
	# Gives us addresses for the logfile and statusfile,
	# and writes the starting date and time to them.
	logfile = results_folder+'/log.txt'	
	statusfile = results_folder+'/status.txt'
	write_to_log("Logfile created: "+str(today))
	write_to_status("Status file created: "+str(today))

def write_to_log(message):
	# Error message can only be one line long, since adding '\n' messes up the call
	# Calling this for each line you want to error-write is a workaround, since it automatically
	# adds one.
	# I think this is more efficient than opening and writing to the end of the file with the 
	# Python I/O tools, but if not it would probably be more convenient to use those.
	call("echo "+message+" >> "+logfile, shell=True)

def write_to_status(message):
	# Same comment as in write_to_log
	call("echo "+message+" >> "+statusfile, shell=True)

if __name__ == "__main__":
	# Parse input.
	args = parse_input_arguments()
	
	# Set up log/status files
	output_dir = args.output_dir
	set_up_output_dir(output_dir)


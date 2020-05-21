import argparse
import sys
import os
import warnings

def parse_input_arguments():
	''' Sets up an argument parser and checks a couple of the inputs for correctness.'''
	
	# Set up the argument parser
	parser = argparse.ArgumentParser(description="QUILTS")
	parser.add_argument('--sample_list', type=str, help="list of samples to run QUILTS on")
	parser.add_argument('--output_dir', type=str, default='.', help="full path to output folder")
	parser.add_argument('--proteome', type=str, help="full path to folder containing reference proteome")
	parser.add_argument('--genome', type=str, help="full path to folder containing reference genome")
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
	
	if not args.sample_list:
		raise SystemExit("ERROR: Must have input list of samples.\nAborting program.")
	
	# Check that we have a somatic and/or germline and/or junction file. Abort if not.
	# Will add fusion to this check later.
	if not args.somatic and not args.germline and not args.junction:
		raise SystemExit("ERROR: Must have at least one variant file (somatic, germline, and/or junctions).\nAborting program.")
	
	if not args.proteome:
		raise SystemExit("ERROR: Must include a directory containing the reference proteome. See README in prepare_proteome folder for details.\nAborting program.")
	
	if not args.genome:
		raise SystemExit("ERROR: Must include a directory containing the reference genome. See README in prepare_genome folder for details.\nAborting program.")
	
	if not args.output_dir:
		print warnings.warn("No output directory given - defaulting to CWD")
	else:
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
	
	return args

#def raise_warning(warn_message):
#	'''Raising a warning.'''
#	my_warning = warnings.warn(warn_message)
#	print my_warning

args = parse_input_arguments()

f = open(args.sample_list,'r')

for f in f.readlines():

	f = f.rstrip()
	cur_dir = os.path.dirname(os.path.realpath(__file__))
		
	out_sh = "%s/%s/quilts.%s.sh" % (args.output_dir, f, f)
	
	w = open(out_sh,'w')

	w.write("python %s/quilts.py --genome %s --proteome %s --output_dir %s" % (cur_dir, args.genome, args.proteome, args.output_dir))

	if args.junction:
		w.write(" --junction %s" % args.junction)

	if args.germline:
		w.write(" --germline %s/%s/germline" % (args.germline, f))

	if args.somatic:
		w.write(" --somatic %s/%s/somatic" % (args.somatic, f))

	if args.threshA:
		w.write(" --threshA" % (args.threshA))

	if args.threshN:
		w.write(" --threshN" % (args.threshN))

	if args.threshAN:
		w.write(" --threshAN" % (args.threshAN))
	
	w.close()
	
	os.system("chmod g+wrx %s/%s/quilts.%s.sh" % (args.output_dir, f, f))
	#os.system("qsub -cwd -S /bin/bash %s/%s/quilts.%s.sh" % (args.output_dir, f, f))
		
	
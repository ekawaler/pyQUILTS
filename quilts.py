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
#from sets import Set
from itertools import product
import sqlite3

# Sorry about the global variables!
global logfile
global referrorfile
global statusfile
global results_folder
global codon_map

class Sample:
    # I think I'll just use the position in the header as the primary key for the samples.
    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.variants = {}

    def add_variant(self, var_pk, zygosity):
        self.variants[var_pk] = zygosity

class Variant:
    def __init__(self, prot_acc,genvar,inprotvar,snpid,genotypes,var_allele_no):
        self.chr = genvar.split('-')[0]
        self.chr_pos = int(re.findall(r'\d+', genvar.split(':')[0].split('-')[1])[0])
        self.ref_nt = genvar.split(':')[0].split('-')[1].split(str(self.chr_pos))[0]
        self.alt_nt = genvar.split(':')[0].split('-')[1].split(str(self.chr_pos))[1]
        self.qual = genvar.split(':')[1]
        self.prot_pos = int(re.findall(r'\d+', inprotvar.split(':')[0].split('-')[1])[0])
        self.prot_accession = prot_acc
        self.snpid = snpid
        self.genotypes = genotypes.split(':')
        self.var_allele_no = var_allele_no
        self.aa_pos = None
        self.ref_aa = None
        self.alt_aa = None
    
    def print_variant(self):
        print_string = "chr%s:%s%d%s" % (self.chr,self.ref_nt, self.chr_pos, self.alt_nt)
        return print_string

    def varsimple_output(self):
        vsstring = "(%d|%s|%s)" % (self.aa_pos, self.alt_aa, self.ref_aa)
        return vsstring
    
    def snp_output(self):
        #snpstring = "(chr%s|%d|%d|%s|%s|Novel)" % (self.chr, self.chr_pos+1, self.chr_pos+1,self.alt_nt,self.ref_nt)
        snpstring = "(chr%s|%d|%d|%s|%s|%s)" % (self.chr, self.chr_pos+1, self.chr_pos+1,self.alt_nt,self.ref_nt,self.snpid)
        return snpstring
    
    def set_variant_aa(self, ref_aa, aa_pos, alt_aa):
        self.ref_aa = ref_aa
        self.aa_pos = aa_pos
        self.alt_aa = alt_aa

class Protein:
    def __init__(self, primary_key,seq):
        self.primary_key = primary_key
        self.sequence = seq
        self.accession = None
        self.prot_name = None
        self.gene_name = None
        self.gene_accession = None
        self.tax_id = 9606 # Shouldn't be hardcoded
        self.tax_name = "Homo_Sapiens" # Shouldn't be hardcoded
        self.transcript_accession = None
        self.seq_length = len(seq)
    
    def process_ensembl_header(self,header):
        #>ENSP00000494411.1 pep chromosome:GRCh38:3:41194741:41239949:1
        #  gene:ENSG00000168036.18 transcript:ENST00000643541.1
        #  gene_biotype:protein_coding transcript_biotype:protein_coding
        #  gene_symbol:CTNNB1 description:catenin beta 1 [Source:HGNC Symbol;Acc:HGNC:2514]
        self.accession = header[1:].split('.')[0]
        self.gene_name = header.split('gene_symbol:')[1].split()[0]
        self.gene_accession = header.split('gene:')[1].split('.')[0]
        self.transcript_accession = header.split('transcript:')[1].split('.')[0]
        if "description" in header:
            self.prot_name = header.split('description:')[1].rstrip()
        else:
            self.prot_name = self.gene_name

    def process_refseq_header(self,header):
        pass

    def process_uniprot_header(self,header):
        pass

def parse_input_arguments():
    ''' Sets up an argument parser and checks a couple of the inputs for correctness.'''
    parser = argparse.ArgumentParser(description="QUILTS") # It's QUILTS!
    parser.add_argument('--output_dir', type=str, default=".", help="full path to output folder (defaults to .)")
    parser.add_argument('--reference_proteome', choices=['ensembl','refseq','uniprot'], default="ensembl", help="Reference proteome source (defaults to ensembl). Accepts ensembl, refseq, or uniprot")
    parser.add_argument('--proteome', type=str, default=".", help="full path to folder containing reference proteome (defaults to .)")
    parser.add_argument('--genome', type=str, default=".", help="full path to folder containing reference genome (defaults to .)")
    parser.add_argument('--somatic', type=str, help="full path to somatic variant VCF file")
    parser.add_argument('--germline', type=str, help="full path to germline variant VCF file")
    parser.add_argument('--variant_quality_threshold', type=float, default=0.0, help="Quality threshold for variants (defaults to 0.0)")
    parser.add_argument('--dbname', type=str, default='index_database.db', help="Name of index database (defaults to index_database.db)")
    
    # Pull out the arguments
    args = parser.parse_args()
    # Check that we have a somatic and/or germline and/or junction file. Abort if not.
    if not args.somatic:
        raise SystemExit("ERROR: Must have a somatic variant file.\nAborting program.")

    if args.somatic.rsplit('/',1)[-1] == 'filtered.vcf':
        raise SystemExit("ERROR: 'filtered.vcf' is a protected name. Please rename your VCF and try again. \nAborting program.")

    if args.dbname[-3:] != '.db':
        args.dbname = args.dbname+'.db'
        
    return args

def set_up_output_dir(output_dir, args):
    ''' Sets up the results folder (creates it, sets up a log file and status file).
        Because this will probably confuse people (me) later, "output directory" refers
        to the directory the user specifies in which the "results folder" will reside. The
        results folder contains the actual results.'''
        
    # I hate this global variable stuff but is there another way to do it
    # without having to pass the logfile/statusfile to all the functions so they can be used?
    global logfile
    global referrorfile
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
    referrorfile = results_folder+'/reference_mismatches.txt'
    write_to_log("Logfile created: "+str(today), logfile)
    write_to_status("Status file created: "+str(today))
    write_to_log("Reference mismatches file created: "+str(today), referrorfile)
    
    # Creates some folders within the results folder
    # Currently just copying what's in the Perl version.
    # Likely to change as I start building out other functionality.
    os.makedirs(results_folder+"/log")
    os.makedirs(results_folder+"/peff")
    os.makedirs(results_folder+"/fasta")
    #os.makedirs(results_folder+"/fasta/parts")
    
    # Moves reference proteome to the working area.
    try:
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

    return results_folder

def save_ref_prot(ref_prot_loc, ref_db_type):
    '''This function saves the proteome as a map of Protein objects.'''
    ref_prot = {}
    f = open(ref_prot_loc,'r')
    header = f.readline()
    seq = ''
    line = f.readline()
    count = 0
    while line:
        if line[0] == '>':
            prot = Protein(count,seq)
            if ref_db_type == 'ensembl':
                prot.process_ensembl_header(header)
            elif ref_db_type == 'refseq':
                prot.process_refseq_header(header)
            elif ref_db_type == 'uniprot':
                prot.process_uniprot_header(header)
            ref_prot[prot.accession] = prot
            header = line
            seq = ''
            count += 1
            line = f.readline()
        else:
            seq += line.rstrip()
            line = f.readline()
    prot = Protein(count,seq)
    if ref_db_type == 'ensembl':
        prot.process_ensembl_header(header)
    elif ref_db_type == 'refseq':
        prot.process_refseq_header(header)
    elif ref_db_type == 'uniprot':
        prot.process_uniprot_header(header)
    ref_prot[prot.accession] = prot # the last one
    return ref_prot

def set_up_db(dbloc):
    conn = sqlite3.connect(dbloc)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    c.execute('''CREATE TABLE proteins
             (PKProteins integer primary key, 
             ProteinAccession text, 
             ProteinName text, 
             GeneName text, 
             GeneAccession text,
             TranscriptAccession text,
             TaxID integer,
             TaxName text,
             PeptideSeq text,
             Length integer)''')
    c.execute('''CREATE TABLE variants
             (PKVariants integer primary key, 
             FKProteins integer, 
             SAAVPos integer, 
             SAAVRef text, 
             SAAVAlt text,
             Chromosome text, 
             SNPPos integer, 
             SNPRef text, 
             SNPAlt text, 
             SNPAccession text, 
             foreign key (FKProteins) references proteins(PKProteins))''')
    c.execute('''CREATE TABLE samples
             (PKSamples integer primary key, 
             SampleName text)''')
    c.execute('''CREATE TABLE samples_to_variants
             (PKSampToVar integer primary key,
             FKSamples integer, 
             FKVariants integer,
             Zygosity text,
             foreign key (FKSamples) references samples(PKSamples),
             foreign key (FKVariants) references variants(PKVariants))''')
    conn.commit()
    return conn

### These functions are for merging and quality checking the variant files.

def set_up_vcf_output_dir(vcf_dir):
    '''Make the directory we'll put our merged VCF files in. Returns with an error if directory doesn't exist.'''
    if not os.path.isdir(vcf_dir):
        # No VCF files are going to be found. Gotta leave the function.
        return 1
    try:
        os.makedirs(vcf_dir+"/merged")
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

def qual_filter(vcf_dir, quality_threshold):
    '''Returns None if successful, or 1 if unable to find variant files. If duplicates appear, pick the one with the highest quality score.'''
    
    # Set up the output directory, if possible. If not found, return with a warning.
    #output_dir_flag = set_up_vcf_output_dir(vcf_dir)
    #if output_dir_flag:
    #    warnings.warn("VCF directory not found at %s" % vcf_dir)
    #    return 1

    # Pull the .vcf file from the directory. If not found, return with a warning.
    #vcf_files = pull_vcf_files(vcf_dir)
    #if len(vcf_files) == 0:
    #    warnings.warn("Unable to find any .vcf files in %s" % vcf_dir)
    #    return 1
    if not os.path.isfile(vcf_dir):
        warnings.warn("Unable to find the .vcf file at %s" % vcf_dir)
        return 1

    # Cool, we have a valid file, let's open our log and output files.
    vcf_directory = vcf_dir.rsplit('/',1)[0]
    vcf_log_location = vcf_directory+'/qualfilter.log'
    vcf_log = open(vcf_log_location,'w')
    w = open(vcf_directory+'/filtered.vcf','w')
    #w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")    

    # Set up some tracking stats etc.
    qual_removed = 0
    total_variants = 0
    kept_variants = 0
    duplicates_removed = 0 # Duplicates in this file that were lower-quality than older duplicates
    old_duplicates_removed = 0 # Duplicates in this file that were higher-quality than older duplicates
    existing_variants = {}

    f = open(vcf_dir, 'r')
    write_to_log(vcf_dir, vcf_log_location)
    
    line = f.readline()
    while line:
        spline = line.split('\t')
        if line[:6] == '#CHROM':
            # It's the important header line! Check for sample names
            w.write(line)
            line = f.readline()
            continue
        if line[0] == '#':
            # It's a header line. Ignore it.
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
            if spline[5].rstrip() == '.':
                spline[5] = str(quality_threshold+.01)
            chr, pos, id, old, new_array, qual = spline[0].lstrip('chr'), int(spline[1])-1, spline[2], spline[3], spline[4].split(','), eval(spline[5])
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
                if existing_variants[line_map_key].split('\t')[5]!= '.' and eval(existing_variants[line_map_key].split('\t')[5]) > qual:
                    duplicates_removed += 1
                    continue
                else:
                    old_duplicates_removed += 1

            # Passed the other checks? Cool, store the line in our map and continue.
            # Overwrite the duplicate if it exists.
            kept_variants += 1
            existing_variants[line_map_key] = line

        line = f.readline()

    # Add everything up and write a status line to the log
    #header_written = True # Don't want to write a header unless it's the first file
    write_to_log("\nFrom variant file %s: %d total variants, %d failed quality check at threshold %f, %d duplicates found \(%d overwritten\), %d variants kept in final version" % (vcf_dir, total_variants, qual_removed, quality_threshold, (duplicates_removed+old_duplicates_removed), old_duplicates_removed, kept_variants), vcf_log_location)
    write_to_log("------------------------", vcf_log_location)
    f.close()
        
    # Write all variants to filtered.vcf
    for key in sorted(existing_variants.keys()):
        #w.write('chr'+'\t'.join(str(i) for i in existing_variants[key])+'\n')
        w.write(existing_variants[key])

    # Write final status line to the log and close all files, we're done.
    #write_to_log("\nTotal for all variant files: %d total variants, %d failed quality check at threshold %f, %d duplicates removed, %d variants kept in final version" % (total_variants_all, qual_removed_all, quality_threshold, duplicates_all, len(existing_variants.keys())), vcf_log_location)
    w.close()
    vcf_log.close()
    return None

def remove_somatic_duplicates(germ_dir, som_dir):
    '''Removes duplicate germline/somatic variants from the somatic file.'''
    # Shouldn't have to check and make sure the files are there - that should have happened earlier
    # Open the germ file, save all of the chr/pos/old/new in a hash table for easy lookup
    # Then open the somatic file, go through each line, if one matches just don't write it back out
    # Afterwards, overwrite the old somatic merged.vcf with the new one
    
    global logfile
    
    # Open the germline file, save all of the information in a hash table for easy lookup, close it
    f = open(germ_dir+'/filtered.vcf','r')
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
    
    w = open(som_dir+'/temp.vcf','w')
    f = open(som_dir+'/filtered.vcf','r')
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
    shutil.move(som_dir+'/temp.vcf', som_dir+'/filtered.vcf')

def samples_to_database(samples, conn):
    c = conn.cursor()

    samp_pk = 1

    for samp in samples:
        c.execute('INSERT INTO samples VALUES (?,?)', 
        (samp_pk,samp))
    
        samp_pk += 1

    conn.commit()

def get_variants(vcf_file, proteome_file, type, conn):
    '''This function is used to create .bed files of only variants that exist in exons.'''
    f = open(proteome_file, 'r')
    
    # Basically goes through proteome.bed and finds all positions in exons
    # Then goes through the .vcf file and checks if they're in any exons
    # Output line: NP_XX (name) \t type-old(position in genome)new:(quality),(repeat)\t type-old(position in exome)new:(quality),(repeat) 
    
    est = ExonSearchTree()
    all_names = set([])
    line = f.readline()
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
    f.close()

    f = open(vcf_file, 'r')
    
    line = f.readline()
    genotypes = False
    header = line.rstrip().split('\t')
    if len(header) > 8:
        # We have genotype info
        genotypes = True
        samples_to_database(header[9:], conn)
    else:
        samples_to_database(["Sample"], conn)

    line = f.readline()
    all_genes = {}
    while line:
        spline = line.rstrip().split('\t')
        try:
            chr, pos, id, old, new, qual = spline[0].lstrip('chr'), int(spline[1]), spline[2], spline[3], spline[4], spline[5]
        except ValueError:
            # Maybe write this to a log somewhere?
            warnings.warn("Failed to parse %s" % line)
            line = f.readline()
            continue
        if not valid_nucleotides(new) or not valid_nucleotides(old):
            warnings.warn("Failed to parse %s" % line)
            line = f.readline()
            continue
        if old == '.' or old == '-':
            old = ''
        if new == '.' or new == '-':
            new = ''
        if id == '.' or id == '':
            id = "Novel"
        if genotypes:
            gtypes = ''
            for samp_no in range(9,len(header)):
                # According to VCF v4.2, these columns are samples and the first pre-colon field is the genotype
                gtypes = gtypes+spline[samp_no].split(':')[0]+':'
            gtypes = gtypes.rstrip(':')
        else:
            gtypes = "1"
        exon = est.find_exon(chr,pos)
        if exon != []:
            # Save pos in chr and pos in gene
            # Each exon returned is a [name, position in gene] pair
            for ex in exon:
                in_chr = "%s-%s%d%s:%s" % (chr, old, pos, new, qual) 
                in_gene = "%s-%s%d%s:%s" % (type, old, ex[1], new, qual) 
                try:
                    all_genes[ex[0]].append([in_chr, in_gene, id, gtypes])
                except KeyError:
                    all_genes[ex[0]] = []
                    all_genes[ex[0]].append([in_chr, in_gene, id, gtypes])
        line = f.readline()
    f.close()

    # Write variants out to file.
    # Need to let it write genes with no variants, also. [Why, though?]
    w = open(proteome_file+".var", 'w')    
    for key in all_names:
        in_chr = []
        in_gene = []
        ids = []
        gtypes = []
        try:
            vars = all_genes[key]
        except KeyError:
            vars = []
        for var in vars:
            in_chr.append(var[0])
            in_gene.append(var[1])
            ids.append(var[2])
            gtypes.append(var[3])
        w.write("%s\t%s\t%s\t%s\t%s\n" % (key, ','.join(in_chr), ','.join(in_gene), ','.join(ids), ','.join(gtypes)))
    w.close()

def valid_nucleotides(strg, search=re.compile(r'[^ACGTU\.\-.]').search):
    '''Thanks, Stack Overflow user! Use this to make sure the variant we're given contains
    only nucleotide letters - found a MAF file that sometimes put "TRUE" in there because
    people are jerks and you can only trust yourself. Also, at some point I might add
    support for the other IUPAC base notations (like N), but currently just throws those out.'''
    return not bool(search(strg))

def process_gene(header_line, second_header, exon_headers, exon_seqs, variants, ref_prot):
    '''Checks all variants in a gene for whether or not they cause a single-AA non-stop substitution. 
    Returns a list of the ones that do, and the AA substitutions they cause.'''
    global codon_map
    global logfile
    global referrorfile
    
    changes = []
    aa_substs = []
    indels = []
    indel_substs = []

    if header_line.split('\t')[3] not in ref_prot.keys():
        write_to_log("Protein %s not found in reference, skipping..." % header_line.split('\t')[3], referrorfile)
        return changes

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

    # With hg38, we get a weird thing where some of the genes in proteome.bed.dna have no sequence associated with them.
    # So...I guess I'll just skip those for now until I figure out why
    # I think this is fixed, and it was a problem in a different script. Keeping this in anyway, just in case.
    if '0' in full_seq  or 'N' in full_seq:
        write_to_log(header_line+full_seq, logfile)
        return changes
    
    prev_start = -1 # Previous codon start position
    prev_subst = [-1, ''] # Previous substitution position and nucleotide
    for var in variants:
        # Remove the ones with multiple nt in the variant
        # I'm going to remove the ones where both sides have the same number but that number is greater than 1
        # Just for now.
        pos = var.prot_pos # using the -1 so we're counting from 0, not 1
        orig_nt = var.ref_nt
        new_nt = var.alt_nt
        
        # Here we removed the ones where it wasn't a 1 to 1 substitution.
        # Now we automatically add it to the changes.
        if len(orig_nt) != 1 or len(new_nt) != 1:
            '''
            if reverse_flag:
                #pos = len(full_seq) - pos
                #orig_nt = orig_nt[::-1].translate(translate_table)
                #new_nt = new_nt[::-1].translate(translate_table)
                qual = var.qual
                prefix = var.split('-')[0]
                indels.append('S-%s%d%s:%s' % (orig_nt, pos, new_nt, qual))
            else:
                indels.append(var)
            '''
            continue
        
        triplet_start = ((pos/3)*3) # start of the triplet containing pos. counting from 0, not 1
        triplet_orig = full_seq[triplet_start:(triplet_start+3)].upper()
        triplet_subst = new_nt
        if triplet_orig == '':
            print('\n'+'Empty codon: '+var.print_variant(), full_seq, len(full_seq), triplet_start, header_line)
            continue
        if len(triplet_orig) != 3:
            print('\n'+'Triplet with length <3 (length of sequence likely not divisible by 3, for whatever reason): '+var.print_variant(), full_seq, len(full_seq), triplet_start, header_line)
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
            write_to_log("Start codon variant found: %s" % (header_line.split('\t')[3]+'\t'+var.print_variant()+'\t'+triplet_orig+'\t'+triplet_new+'\t'+full_seq[subst_pos-3:subst_pos+6]),logfile)
            continue
        # If neither the old nor the new codon is a stop codon and the AA changes,
        # count it!
        #if AA_old != AA_new and AA_new != '*' and AA_old != '*':
        # Now we're just checking whether the old and new AAs are different. Leaving stop codon things in...
        if AA_old != AA_new:
            if reverse_flag:
                total_AA = len(full_seq)/3
                position = total_AA-(pos/3)
            else:
                position = pos/3+1
            
            try:
                orig_seq = ref_prot[header_line.split('\t')[3]].sequence
                orig_aa = orig_seq[position-1]
                if orig_aa == AA_old:
                    var.set_variant_aa(AA_old, position, AA_new)
                    changes.append(var)
                elif (AA_old=='*') and (orig_aa in 'EOQUW') and (orig_aa != AA_new):
                    var.set_variant_aa(orig_aa, position, AA_new)
                    changes.append(var)
                else:
                    write_to_log("Original AA at position %d in %s not what we expected [expected %s, found %s], skipping..." % (position, header_line.split('\t')[3], AA_old, orig_seq[position-1]), referrorfile)
            #except KeyError:
            #    write_to_log("Protein %s not found in reference, skipping..." % header_line.split('\t')[3], referrorfile)
            except IndexError:
                orig_seq = ref_prot[header_line.split('\t')[3]].sequence
                if (AA_old=='*') and (position==len(orig_seq)+1):
                    var.set_variant_aa(AA_old, position, AA_new)
                    changes.append(var)
                else:
                    write_to_log("Protein %s not expected length, skipping..." % header_line.split('\t')[3], referrorfile)

        # If the new codon is a stop codon, 
        # uh...I found this comment here and it looks like I decided not to do whatever I was planning to do.
        prev_start = triplet_start
        prev_subst = [subst_pos, triplet_subst]
    return changes    

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
            start_sec += snp_pos-tot_len+1
            return "%s%d%s" % (orig_nt,start_sec,new_nt)
        else:
            tot_len += lengths[i]

def translate_seq(sequence, strand, return_all = False):
    '''Translates a DNA sequence to an AA sequence'''
    global codon_map
    
    if strand == '-':
        translate_table = maketrans("ACGTacgt","TGCAtgca")
        sequence = sequence[::-1].translate(translate_table)
    
    translated = ''
    
    for i in range(0,len(sequence)-2,3):
        translated += codon_map.get(sequence[i:i+3],'X')
    
    # Now this should return the translated sequence up to the first stop codon.    
    # That could either be at the actual end of the sequence, or at a variant that adds a stop.
    if return_all:
        # In some cases, I want to return the whole sequence, including internal stops.
        return translated
    else:
        return translated.split('*')[0]

def write_out_peff(variants,ref_prot):

    # Set up the PEFF
    output_peff = open(results_folder+"/peff/variant_proteome.peff",'w')

    output_peff.write("# PEFF 1.0\n")
    output_peff.write("# DbName=ENSEMBL38.100\n")
    output_peff.write("# Prefix=db\n")    
    output_peff.write("# DbDescription=ENSEMBL38.100\n")
    #output_peff.write('# CustomKeyDef=(KeyName=SNP|Description="Single nucleotide polymorphism information"| RegExp="((chr)?[A-Za-z]?[0-9]{0,2})\| ([0-9]+)\|[ACGTUNX]*\|[ACGTUNX]*)"|FieldNames=Chromosome,StartPosition,NewNTSequence,OriginalNTSequence|FieldTypes=string,integer,string,string)\n')
    output_peff.write('# CustomKeyDef=(KeyName=PKProtein|Description="Primary key for the protein in the index database"|RegExp="[0-9]+"|FieldNames=PrimaryKey|FieldTypes=integer)\n')
    output_peff.write("# DbVersion=38.100\n")
    output_peff.write("# DbSource=ENSEMBL, CPTAC\n")
    output_peff.write("# NumberOfEntries=%d\n" % len(ref_prot.keys()))
    output_peff.write("# SequenceType=AA\n")
    output_peff.write("# //\n")

    for prot in ref_prot.keys():
        protein = ref_prot[prot]
        varsimp_string = ''
        snp_string = ''
        vars_written = []

        output_peff.write('>db:%s \\PName=%s \\GName=%s \\NcbiTaxId=%d \\TaxName=%s \\Length=%d \\PKProtein=%d' % (prot,protein.prot_name, protein.gene_name, protein.tax_id, protein.tax_name, protein.seq_length, protein.primary_key))
        if prot in variants.keys():
            varsimp_string = ' \VariantSimple='
            snp_string = ' \SNP='
            vars = variants[prot]
            for v in vars:
                if v.varsimple_output() not in vars_written:
                    varsimp_string+=v.varsimple_output()
                    snp_string+=v.snp_output()
                    vars_written.append(v.varsimple_output())
            #output_peff.write('%s%s' % (varsimp_string, snp_string))
            output_peff.write('%s' % (varsimp_string)) # SNP string is in the index DB now

        output_peff.write('\n%s\n' % protein.sequence)
    
    output_peff.close()

def fill_index_db(variants, ref_prot, conn):

    c = conn.cursor()

    var_pk = 1
    s_to_v_pk = 1

    # For each protein...

    for prot in ref_prot.keys():
        protein = ref_prot[prot]
        c.execute('INSERT INTO proteins VALUES (?,?,?,?,?,?,?,?,?,?)', 
        (protein.primary_key,prot,protein.prot_name,protein.gene_name,protein.gene_accession,protein.transcript_accession,protein.tax_id, protein.tax_name, protein.sequence,protein.seq_length))
    
        if prot in variants.keys():
            vars = variants[prot]
            for var in vars:
                c.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?)', 
                (var_pk,protein.primary_key,var.aa_pos,var.ref_aa,var.alt_aa,var.chr,var.chr_pos,var.ref_nt,var.alt_nt,var.snpid))
                for v in range(len(var.genotypes)):
                    c.execute('INSERT INTO samples_to_variants VALUES (?,?,?,?)', (s_to_v_pk, v+1,var_pk,var.genotypes[v]))
                    s_to_v_pk += 1
                var_pk += 1

    conn.commit()

def sort_variants(proteome_file, variant_file, output_prefix, ref_prot):
    '''Goes through the variants and sorts them by type, writing out a bunch of intermediate files in the process.'''
    global codon_map
    
    # Grab all the variants and their locations within their genes
    variants_map = {}
    f = open(variant_file, 'r')
    line = f.readline()
    while line:
        spline = line.rstrip().split('\t')
        if len(spline) > 1:
            split_variants = spline[2].split(',')
            split_genomevar = spline[1].split(',')
            split_names = spline[3].split(',')
            split_samples = spline[4].split(',')
            for v in range(len(split_variants)):
                saav = Variant(spline[0],split_genomevar[v],split_variants[v],split_names[v],split_samples[v], v+1)
                variants_map.setdefault(spline[0], []).append(saav)
        line = f.readline()
    f.close()
    
    # Open some files to write to. Right now, only looking at single AA changes.
    # I hate these filenames. Do they actually mean anything?
    #out_aa = open(file_base+"/"+output_prefix+".bed.aa.var", 'w')
    # This one will basically just be that first line. Might want an extra line for each variant though, depends how it's used. Actually, almost certainly do.
    #out_aa_bed = open(file_base+"/"+output_prefix+".aa.var.bed", 'w')
    # proteome.bed.dna, but with the variants included. Def want an extra line for each variant.
    #out_aa_bed_dna = open(file_base+"/"+output_prefix+".aa.var.bed.dna", 'w')
    # Right now just proteome.bed.aa.var but with the non-AA variants. Doesn't do anything right now because I don't care about it very much.
    #out_indel = open(file_base+"/"+output_prefix+".bed.indel.var", 'w')
    #out_indel_bed = open(file_base+"/"+output_prefix+".indel.var.bed", 'w')
    #out_indel_bed_dna = open(file_base+"/"+output_prefix+".indel.var.bed.dna", 'w')
    #out_other = open(file_base+"/"+output_prefix+".bed.other.var", 'w')
    #out_peff = open(file_base+"/"+output_prefix+".peff",'w')
    
    changed_var_map = {}

    f = open(proteome_file, 'r')
    line = f.readline()
    while line:
        # Line type one: First gene header line
        if line[:3] == 'chr':
            # Finish processing the previous gene - focus first on proteome.bed.aa.var
            # Using a try/except block here to deal with the possibility of this being the first line
            # There has to be a more graceful way to do that, right?
            try:
                changed_vars = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []), ref_prot)
                if changed_vars != []:
                    changed_var_map[name] = changed_vars
                #if indels != []:
                    # Make sure none of them are just removing intron stuff
                    #indels = test_indels(exon_headers, indels)
                    # If there are any left...
                    #if indels != []:
                        #out_indel.write("%s\t%s\n" % (name, ','.join(indels)))
                        #out_indel_bed.write("chr%s\t%d\t%d\t%s-indel\t1000\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (chr, start, end, name, strand, start, end, spline[-4], exon_count, ','.join(exon_lengths), ','.join(exon_offsets)))
                        #write_out_indel_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, indels, out_indel_bed_dna)
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
        changed_vars = process_gene(header_line, second_header, exon_headers, exon_seqs, variants_map.get(name, []), ref_prot)
        if changed_vars != []:
            changed_var_map[name] = changed_vars
        #if indels != []:
            # Make sure none of them are just removing intron stuff
            #indels = test_indels(exon_headers, indels)
            # If there are any left...
            #if indels != []:
                #out_indel.write("%s\t%s\n" % (name, ','.join(indels)))
                #out_indel_bed.write("chr%s\t%d\t%d\t%s-indel\t1000\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (chr, start, end, name, strand, start, end, spline[-4], exon_count, ','.join(exon_lengths), ','.join(exon_offsets)))
                #write_out_indel_bed_dna(header_line.split('\t'), second_header.split(), exon_headers, indels, out_indel_bed_dna)
    except UnboundLocalError:
        pass

    # Close everything
    f.close()
    #out_aa.close()
    #out_aa_bed.close()
    #out_aa_bed_dna.close()
    #out_indel.close()
    #out_indel_bed.close()
    #out_indel_bed_dna.close()
    #out_other.close()

    return changed_var_map

def tabulate_ref_mismatches():
    '''Opens up the reference mismatch folder and counts them up to add to the logfile'''
    global logfile
    global referrorfile

    missing_prots = set([])
    mismatch_prots = set([])
    rerf = open(referrorfile,'r')
    for line in rerf.readlines():
        if line.split()[0] == 'Protein':
            missing_prots.add(line.split()[1])
        elif line.split()[0] == "Original":
            mismatch_prots.add(line.split()[6])
    rerf.close()
    
    write_to_log("%d SAAV proteins not found in reference" % (len(missing_prots)),logfile)
    write_to_log("%d SAAV proteins do not match between the transcriptome and the proteome reference" % (len(mismatch_prots)),logfile)

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
    '''Raises a slightly less hideous warning than the default. Only used once, might delete later.'''
    my_warning = warnings.warn(warn_message)
    print(my_warning)

def quit_if_no_variant_files(args):
    if not args.somatic:
        raise SystemExit("ERROR: Couldn't find somatic variant file!\nAborting program.")

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
    
    # Open up index table
    conn = set_up_db(results_folder+'/'+args.dbname)

    # Set up codon map
    codon_map = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S",
    "TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T",
    "ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*",
    "TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R",
    "CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G",
    "GGA":"G","GGG":"G"}
    
    # Prep a map of the reference proteome for quality checks
    ref_prot = save_ref_prot(results_folder+"/fasta/reference_proteome.fasta",args.reference_proteome)

    # Time to merge and quality-threshold the variant files!
    if args.somatic:
        som_flag = qual_filter(args.somatic, args.variant_quality_threshold)
        if som_flag:
            args.somatic = None
    if args.germline:
        germ_flag = qual_filter(args.germline, args.variant_quality_threshold)
        if germ_flag:
            args.germline = None
    quit_if_no_variant_files(args) # Check to make sure we have our somatic variant file
    write_to_status("Merge and qual filter finished")
        
    # Now let's remove everything in the somatic file that is duplicated in the germline file
    # since if it shows up in both, it's probably actually a germline variant.
    if args.somatic and args.germline:
        remove_somatic_duplicates(args.germline.rsplit('/',1)[0], args.somatic.rsplit('/',1)[0])
    write_to_status("Probable germline variants filtered out of somatic variant file")

    # Call read_chr_bed.c, which takes the reference genome and proteome.bed file as input and produces
    # a fasta file of exomes.
    # Possible but unlikely future work: rewrite the C file (still in C though) so it's more efficient?
    # I dunno, it seems fine for now.
    try:
        check_call("%s/read_chr_bed %s/log/proteome.bed %s" % (script_dir, results_folder, args.genome), shell=True)
    except CalledProcessError:
        raise SystemExit("ERROR: read_chr_bed didn't work - now we don't have a proteome.bed.dna file. Try recompiling read_chr_bed.c.\nAborting program.")
    
    # Next, create a proteome.bed file containing only variants...probably.
    # I could probably combine them more prettily, but for now I'll just concatenate the files.
    # Now includes the name of the variant, if it has one
    get_variants(args.somatic.rsplit('/',1)[0]+"/filtered.vcf", results_folder+"/log/proteome.bed", "S", conn)
    write_to_status("Get variants completed, proteome.bed file written")    

    # Combine (some more) and sort variants.
    # Also variants are sorted by what they do to the sequence (add/remove stop, change AA, etc)
    final_vars = sort_variants(results_folder+"/log/proteome.bed.dna", results_folder+"/log/proteome.bed.var", "proteome", ref_prot)
    write_to_status("Variants sorted")    

    # Translate the variant sequences into a fasta file.
    #translate_saavs(results_folder+"/log/", "proteome.aa.var.bed.dna", logfile, ref_prot)
    #write_to_status("Translated SAAVs")    
    #translate(results_folder+"/log/", "proteome.indel.var.bed.dna", logfile, 'indel')
    #write_to_status("Translated indels")
    #shutil.copy(results_folder+"/log/proteome.aa.var.bed.dna.fasta", results_folder+"/fasta/variant_proteome.fasta")
    #shutil.copy(results_folder+"/log/proteome.indel.var.bed.dna.fasta", results_folder+"/fasta/parts/proteome.indel.fasta")

    # Write out the PEFF file
    write_out_peff(final_vars, ref_prot)

    # Write out the index database
    fill_index_db(final_vars, ref_prot, conn)

    # Check to see how many reference mismatches we have
    tabulate_ref_mismatches()
        
    # Combine all of the proteome parts into one.
    write_to_status("DONE")
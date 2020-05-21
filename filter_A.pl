#!/usr/local/bin/perl
#
# Step#1: ensembl_gtf_to_bed.pl
# Step#2: merge_junctions.pl
# Step#3: filter_known_transcripts.pl
# Step#4: filter_A.pl
# Step#5: protein_seq_from_genome_using_bed.pl
# Step#6: protein_seq_from_genome_using_junctions.pl

# This script takes as input 2 *.bed files: 
# 1. Filtered junctions from RNA-Seq data aligned to the genome (*-merged-junctions.filter.bed) and
# 2. A collection of gene models (e.g. Ensembl) where each tuple for a particular transcript has the structure described here: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
# On the second run of the script 3rd intput is "1" in order to create an output .bed file for each protein.

# Output files:
# 1. *-merged-junctions.filter.bed.A.bed: alternative splicing models;
# 2. *-merged-junctions.filter.bed.AN.bed: alternative splicing/novel expression;
# 3. *-merged-junctions.filter.bed.AN_.bed: alternative splicing/novel expression (different RGB for the new genes);
# 4. *-merged-junctions.filter.bed.notA.bed: novel expressions. 
use strict;
use Data::Dumper;

my $error=0;
my $filename1="";
my $filename2="";
my $print_details=0;
my $threshold_A=2;
my $threshold_AN=3;
my $threshold_N=3;

if ($ARGV[0]=~/\w/) { $filename1=$ARGV[0];} else { $error=1; } # RNAseq
if ($ARGV[1]=~/\w/) { $filename2=$ARGV[1];} else { $error=1; } # gene models
if ($ARGV[2]=~/\w/) { $print_details=$ARGV[2];} else { $print_details=0; } # gene models
if ($ARGV[3]=~/\w/) { $threshold_A=$ARGV[3];} else { $threshold_A=2 } 
if ($ARGV[4]=~/\w/) { $threshold_AN=$ARGV[4];} else { $threshold_AN=3; } 
if ($ARGV[5]=~/\w/) { $threshold_N=$ARGV[5];} else { $threshold_N=3; } 
print qq!filter_A.pl: $threshold_A $threshold_AN $threshold_N\n!;
$filename1=~s/\\/\//g;
$filename2=~s/\\/\//g;
my $filename1_=$filename1;
$filename1_=~s/\.bed//i;
if ($print_details==1) { mkdir($filename1_); }
my $dir=$filename1;
if ($dir!~/\/([\/]+)$/) { $dir="."; }
my $line = "";
my $count_model_exons=0;
my $count_model_begin=0;
my $count_model_end=0;
my $chromosome = "";
my %stat = ();
my $num_observ_max=0;
my %junctions_model=();
my %junctions_model_A=();
my %junctions_model_begin_intron=();
my %junctions_model_end_intron=();
my %models=();
my $count_junctions=0;
my %count_junctions=();

if ($error==0)
{
	# goes over the DNA database file(s) and creates 5 hashes:
	# "models" for each protein (key) stores the tuple from DNA database for that protein;
	# "junction_model_begin_intron", "junction_model_end_intron", "junction_model", and
	# "junction_model_A" that has alternative splicing models for each protein;
	if(open (IN_DNA, qq!$filename2!))
	{
		my $line_number = 0;
		while($line = <IN_DNA>)
		{
			chomp($line);
			if ($line=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				$chromosome = $1;
				my $begin_exon = $2;
				my $end_exon = $3;
				my $protein = $4;
				my $strand = $6;
				my $block_count = $10;
				my $block_sizes = $11;
				my $block_starts = $12; 
				if ($print_details==1)
				{
					if (open(OUT_DETAILS,">$filename1_/$protein.bed"))
					{
						print OUT_DETAILS qq!$line\n!;
						close(OUT_DETAILS);
					}
				}
				$models{$protein}=$line;
				my @sizes = split(',',$block_sizes);
				my @starts = split(',',$block_starts);
				my $start = $begin_exon;
				for (my $i = 0; $i < $block_count-1; $i++)
				{
					my $j=$i+1;
					my $begin_intron = $start + $sizes[$i] + 1;
					my $end_intron = $begin_exon + $starts[$j] - 1; 
					$junctions_model{qq!$chromosome#$begin_intron#$end_intron!}.="#$protein,$i,$j#";
					$junctions_model_begin_intron{qq!$chromosome#$begin_intron!}.="#$protein,$i#";
					$junctions_model_end_intron{qq!$chromosome#$end_intron!}.="#$protein,$j#";
					#corrected $block_count-1 to $block_count;
					for (my $j = $i+2; $j < $block_count; $j++)
					{
						my $end_intron_ = $begin_exon + $starts[$j] - 1; 
						$junctions_model_A{qq!$chromosome#$begin_intron#$end_intron_!}.="#$protein,$i,$j#";
					}
					$start = $begin_exon + $starts[$i+1];
					
				}
				$count_model_exons+=$block_count;
			}
			$line_number++;
		}
		close(IN_DNA);
	}
	
	# takes one RNAseq intron at a time from .bed files and compares it to the model introns;
	# prints matched RNAseq introns into files;
	open (OUT, qq!>$filename1.A.bed!);
	open (OUT_, qq!>$filename1.A_.bed!);
	open (OUT_AN, qq!>$filename1.AN.bed!);
	open (OUT_AN_, qq!>$filename1.AN_.bed!);
	open (OUT_NOT, qq!>$filename1.notA.bed!);
	open (OUT_EXTRA, qq!>$filename1.A.extra.bed!);
	open (OUT_EXTRA_, qq!>$filename1.A_.extra.bed!);
	open (OUT_AN_EXTRA, qq!>$filename1.AN.extra.bed!);
	open (OUT_AN_EXTRA_, qq!>$filename1.AN_.extra.bed!);
	open (OUT_NOT_EXTRA, qq!>$filename1.notA.extra.bed!);
	if(open (IN_RNA, qq!$filename1!))
	{
		while($line = <IN_RNA>)
		{
			chomp($line);
			if ($line=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				$chromosome= $1;
				my $begin_exon_1= $2;
				my $end_exon_2= $3;
				my $junc_number = $4;
				my $num_observ = $5; 
				my $length_exons = $11;
				my $begin_exons = $12;
				my @length_exon = split(',', $length_exons,2);
				my @begin_exon = split(',', $begin_exons,2);
				my $size = ();
				my $non_matches = 0;
				my $junction_status="";
				my $junction_name_mod="";
				my $name_2="";
				my $N=1; 
				if ($num_observ_max<$num_observ) { $num_observ_max=$num_observ; }
				if ($num_observ >= 1) 
				{	
					my $begin_intron = $begin_exon_1 + $length_exon[0] + 1;
					my $end_intron = $begin_exon_1 + $begin_exon[-1] - 1;
					my $stop_exon1=$begin_intron-1;
					my $start_exon2=$end_intron+1; 
					my $length_intron = $end_intron - $begin_intron;
					# never executes
					if ($junctions_model{qq!$chromosome#$begin_intron#$end_intron!}=~/\w/)
					{
						$junction_status="C";
						$N=0; 
					}
					else
					{	# if the RNAseq intron was matched to any of the beginnings and ends of the introns from the gene,
						# the new gene is created with the exons that preceed and follow RNAseq intron;
						if ($junctions_model_A{qq!$chromosome#$begin_intron#$end_intron!}=~/\w/)
						{
							$junction_status="A";
							$N=0; 
							my $temp=$junctions_model_A{qq!$chromosome#$begin_intron#$end_intron!};
							while ($temp=~s/^#([^\,]+)\,([^\,]+)\,([^#]+)#//)
							{
								my $protein=$1;
								my $alt1=$2; #exon at beginning
								my $alt2=$3; #exon at end
								#print qq!$chromosome#$begin_intron#$end_intron: $protein $alt1 $alt2\n!;
								if ($models{$protein}=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
								{
									my $chr = $1;
									my $start=$2;
									my $end=$3;
									my $beginning = "$chr\t$start\t$end"; 
									my $name = $4;
									my $middle = "$6\t$7\t$8";
									my $block_count = $10;
									my $block_sizes = $11;
									my $block_starts = $12; 
									my $block_count_ = 0;
									my $block_sizes_ = "";
									my $block_starts_ = ""; 
									my @sizes = split(',',$block_sizes);
									my @starts = split(',',$block_starts);
									my $start_exon1=""; #KVR
									my $size_exon1 = "";
									my $size_exon2=""; 
									my $removed_length=0;
									#my $stop_exon2="";  #KVR
									for (my $i = 0; $i < $block_count; $i++)
									{
										if ($i<=$alt1 or $alt2<=$i)
										{
											$block_count_++;
											$block_sizes_.="$sizes[$i],";
											$block_starts_.="$starts[$i],";
										}
										else
										{
											$removed_length+=$sizes[$i];
										}
										if ($i==($alt1)) #KVR
										{
											$start_exon1=$start+$starts[$i]; #the start location of the exon 1 that bridges
											$size_exon1=$sizes[$i];
										}
										if ($i==($alt2))#KVR
										{
											#$stop_exon2=$start+$starts[$i]+$sizes[$i]; #the end location of the exon2 that bridges
											$size_exon2=$sizes[$i];
										}
									}
									#my $name_.="$name-A-$num_observ-$chromosome-$start_exon1-$begin_exon_1-$end_exon_2-$stop_exon2"; #KVR
									#my $name_.="$name-A-$alt1-$alt2-$num_observ-$junc_number-$chromosome-$start_exon1-$size_exon1-$start_exon2-$size_exon2"; #KVR
									my $name_="$name-A-";
									$name_.=$alt1+1;
									$name_.="-";
									$name_.=$alt2+1;
									$name_.="-$num_observ-$chromosome-";
									$name_.=$start_exon1+$size_exon1;
									$name_.="-$start_exon2-$removed_length"; 
									$block_sizes_=~s/\,$//;
									$block_starts_=~s/\,$//;
									if ($removed_length % 3 == 0)
									{
										if ($threshold_A<=$num_observ)
										{
											print OUT qq!$beginning\t$name_\t$num_observ\t$middle\t0,0,255\t$block_count_\t$block_sizes_\t$block_starts_\n!;
										}
										else
										{
											print OUT_EXTRA qq!$beginning\t$name_\t$num_observ\t$middle\t0,0,255\t$block_count_\t$block_sizes_\t$block_starts_\n!;
										}
									}
									else
									{
										if ($threshold_A<=$num_observ)
										{
											print OUT_ qq!$beginning\t$name_\t$num_observ\t$middle\t0,0,255\t$block_count_\t$block_sizes_\t$block_starts_\n!;
										}
										else
										{
											print OUT_EXTRA_ qq!$beginning\t$name_\t$num_observ\t$middle\t0,0,255\t$block_count_\t$block_sizes_\t$block_starts_\n!;
										}
									}
									if ($print_details==1)
									{
										if (open(OUT_DETAILS,">>$filename1_/$name.bed"))
										{
											print OUT_DETAILS qq!$line\n!;
											print OUT_DETAILS qq!$beginning\t$junc_number-$name_\t$num_observ\t$middle\t0,0,255\t$block_count_\t$block_sizes_\t$block_starts_\n!;
											close(OUT_DETAILS);
										}
										$name_2=$name_; 
									}
								}
							}
						}
						else
						{	# if only the biginning of the RNAseq intron was matched to the existing in database beginning of the DNA intron
							if ($junctions_model_begin_intron{qq!$chromosome#$begin_intron!}=~/\w/)
							{
								$junction_status="AN1";
								$N=0; 
								my $outside=1;
								my $temp=$junctions_model_begin_intron{qq!$chromosome#$begin_intron!};
								while ($temp=~s/^#([^\,]+)\,([^#]+)#//)
								{
									my $protein=$1;
									my $alt1=$2;
									#print qq!AN1: $chromosome#$begin_intron#$end_intron: $protein $alt1\n!;
									if ($models{$protein}=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
									{
										my $chr = $1;
										my $start = $2;
										my $end = $3;
										my $name = $4;
										my $middle = "$6\t$7\t$8";
										my $strand=$6;
										my $block_count = $10;
										my $block_sizes = $11;
										my $block_starts = $12;
										my $start_exon1=""; #KVR
										my $size_exon1 = ""; #KVR
										my $size_exon2="";  #KVR
										#my $stop_exon2="";  #KVR
										if ($end_intron<$end) #within the same gene 
										{
											$outside=0;
											my $block_count_ = 0;
											my $block_sizes_ = "";
											my $block_starts_ = ""; 
											my @sizes = split(',',$block_sizes);
											my @starts = split(',',$block_starts);
											my $modify=1;
											for (my $i = 0; $i < $block_count; $i++)
											{
												# creates new gene from the exons that preceed the beginning of the RNAseq intron and
												# exons that follow after the end of the RNAseq intron;
												if ($i==$alt1) #KVR when we are at the right start intron
												{
													$start_exon1=$start+$starts[$i]; #KVR - start for the exon1
													$size_exon1=$sizes[$i];
												}
												if ($i<=$alt1 or $start+$starts[$i-1]+$sizes[$i-1]>$end_intron) #on both sides of the junction add in the sizes and starts 
												{
													$block_count_++;
													$block_sizes_.="$sizes[$i],";
													$block_starts_.="$starts[$i],";
												}
												else
												{	# if the end of the RNAseq intron is inbetween the ends of the [i]th and [i+1]th exons, modifies [i+1] exon
													# to fit the boundaries of the RNAseq intron; 
													if ($start+$starts[$i]+$sizes[$i]<$end_intron) #KVR WHAT DOES THIS DO??
													{
														my $e1=$start+$starts[$i];
														my $e2=$start+$starts[$i]+$sizes[$i];
														#print qq!Skipped $i. $e1-$e2 $end_intron\n!;
													}
													
													else
													{
														if ($modify==1)
														{
															my $i1=($start+$starts[$i]+$sizes[$i])-$end_intron; #new size 
															my $i2=$end_intron-$start+1; #new start
															#$stop_exon2=$start+$i2+$i1; #KVR end of exon2
															$size_exon2=$i1; 
															$block_count_++;
															$block_sizes_.="$i1,";
															$block_starts_.="$i2,";
															$modify=0;
														}
													}
												}
											}
											if ($junction_name_mod!~/\w/)
											{
												$junction_name_mod="$name-AN1-$alt1";
											}
											#my $name_.="$name-AN1-$junc_number-$num_observ-$chromosome-$start_exon1-$begin_exon_1-$end_exon_2-$stop_exon2";
											#my $name_.="$name-AN1-$num_observ-$junc_number-$chromosome-$start_exon1-$size_exon1-$start_exon2-$size_exon2"; #KVR
											my $name_="$name-AN1-";
											$name_.=$alt1+1;
											$name_.="-$num_observ-$chromosome-";
											$name_.=$start_exon1+$size_exon1;
											$name_.="-$start_exon2"; 
											$block_sizes_=~s/\,$//;
											$block_starts_=~s/\,$//;
											if ($modify==1) { print qq!Error: $name_\n!; }
											else
											{
												if ($strand=~/^\+$/)
												{
													if ($threshold_AN<=$num_observ)
													{
														print OUT_AN qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t0,100,0\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
													else
													{
														print OUT_AN_EXTRA qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t0,100,0\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
												}
												else
												{
													if ($threshold_AN<=$num_observ)
													{
														print OUT_AN_ qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t107,142,35\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
													else
													{
														print OUT_AN_EXTRA_ qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t107,142,35\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
												}
												if ($print_details==1)
												{
													if (open(OUT_DETAILS,">>$filename1_/$name.bed"))
													{
														print OUT_DETAILS qq!$line\n!;
														if ($strand=~/^\+$/)
														{
															print OUT_DETAILS qq!$chr\t$start\t$end\t$junc_number-$name_\t$num_observ\t$middle\t0,100,0\t$block_count_\t$block_sizes_\t$block_starts_\n!;
														}
														else
														{
															print OUT_DETAILS qq!$chr\t$start\t$end\t$junc_number-$name_\t$num_observ\t$middle\t107,142,35\t$block_count_\t$block_sizes_\t$block_starts_\n!;
														}
														close(OUT_DETAILS);
													}
												}
												$name_2=$name_; 
											}
										}
									}
								}
								if ($outside==1) { $junction_status=""; }
							}
							# if only the end of the RNAseq intron was matched to the existing in database end of the DNA intron
							if ($junctions_model_end_intron{qq!$chromosome#$end_intron!}=~/\w/)
							{
								$junction_status="AN2";
								$N=0; 
								my $outside=1;
								my $temp=$junctions_model_end_intron{qq!$chromosome#$end_intron!};
								while ($temp=~s/^#([^\,]+)\,([^#]+)#//)
								{
									my $protein=$1;
									my $alt2=$2;
									#print qq!AN2: $chromosome#$begin_intron#$end_intron: $protein $alt2\n!;
									if ($models{$protein}=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
									{
										my $chr = $1;
										my $start = $2;
										my $end = $3;
										my $name = $4;
										my $middle = "$6\t$7\t$8";
										my $strand=$6;
										my $block_count = $10;
										my $block_sizes = $11;
										my $block_starts = $12;
										my $start_exon1=""; #KVR
										my $size_exon1 = ""; #KVR
										my $size_exon2="";  #KVR
										#my $stop_exon2="";  #KVR
										if ($start<$begin_intron)
										{
											$outside=0;
											my $block_count_ = 0;
											my $block_sizes_ = "";
											my $block_starts_ = ""; 
											my @sizes = split(',',$block_sizes);
											my @starts = split(',',$block_starts);
											my $modify=1;
											my $offset=0;
											for (my $i = 0; $i < $block_count; $i++)
											{
												if ($i==$alt2) #KVR when we are at the right start intron
												{
													$size_exon2=$sizes[$i]; #KVR - end of exon 2
												}
												# checks if the beginning of the RNAseq intron is before the beginning of the DNA gene;
												if($i==0 and $begin_intron<$start)
												{
													$block_count_++;
													my $e=$sizes[$i]+201;
													$offset=$start-$begin_intron+$e;
													$block_sizes_.="$e,";
													$block_starts_.="$starts[$i],";
													#print qq!First $i.\n!;
												}
												else
												{
													# creates new gene from the exons that preceed the beginning of the RNAseq intron and
													# exons that follow after the end of the RNAseq intron;
													if ($i>=$alt2 or $start+$starts[$i+1]<$begin_intron)
													{
														$block_count_++;
														$block_sizes_.="$sizes[$i],";
														my $e=$starts[$i]+$offset;
														$block_starts_.="$e,";
														#print qq!Same $i. $i>$alt2 or $start+$starts[$i-1]<$begin_intron\n!;
													}
													else
													{	# if the beginning of the RNAseq intron is inbetween the beginnings of the [i]th and [i+1]th exons, modifies [i]th exon
														# to fit the boundaries of the RNAseq intron; 
														if ($start+$starts[$i]>$begin_intron)
														{
															my $e1=$start+$starts[$i];
															my $e2=$start+$starts[$i]+$sizes[$i];
															#print qq!Skipped $i. $e1-$e2 $begin_intron\n!;
														}
														else
														{
															if ($modify==1)
															{
																my $i1=$begin_intron-($start+$starts[$i]);
																my $i2=$starts[$i];#exon start
																$start_exon1=$i2 + $start;
																$size_exon1=$i1; 
																$block_count_++;
																$block_sizes_.="$i1,";
																$block_starts_.="$i2,";
																$modify=0;
																#print qq!Modified $i.\n!;
															}
														}
													}
												}
											}
											$start-=$offset;
											if ($junction_name_mod!~/\w/)
											{
												$junction_name_mod="$name-AN2-$alt2";
											}
											#my $name_="$name-AN2-$junc_number-$num_observ-$chromosome-$start_exon1-$begin_exon_1-$end_exon_2-$stop_exon2";
											#my $name_.="$name-AN2-$num_observ-$junc_number-$chromosome-$start_exon1-$size_exon1-$start_exon2-$size_exon2"; #KVR
											my $name_="$name-AN2-";
											$name_.=$alt2+1;
											$name_.="-$num_observ-$chromosome-";
											$name_.=$start_exon1+$size_exon1;
											$name_.="-$start_exon2"; 
											$block_sizes_=~s/\,$//;
											$block_starts_=~s/\,$//;
											$block_sizes_=~s/\,$//;
											$block_starts_=~s/\,$//;
											if ($modify==1) { print qq!Error: $name_\n!; }
											else
											{
												if ($strand=~/^\-$/)
												{
													if ($threshold_AN<=$num_observ)
													{
														print OUT_AN qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t184,134,11\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
													else
													{
														print OUT_AN_EXTRA qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t184,134,11\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
												}
												else
												{
													if ($threshold_AN<=$num_observ)
													{
														print OUT_AN_ qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t255,215,0\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
													else
													{
														print OUT_AN_EXTRA_ qq!$chr\t$start\t$end\t$name_\t$num_observ\t$middle\t255,215,0\t$block_count_\t$block_sizes_\t$block_starts_\n!;
													}
												}
												if ($print_details==1)
												{
													if (open(OUT_DETAILS,">>$filename1_/$name.bed"))
													{
														print OUT_DETAILS qq!$line\n!;
														if ($strand=~/^\-$/)
														{
															print OUT_DETAILS qq!$chr\t$start\t$end\t$junc_number-$name_\t$num_observ\t$middle\t184,134,11\t$block_count_\t$block_sizes_\t$block_starts_\n!;
														}
														else
														{
															print OUT_DETAILS qq!$chr\t$start\t$end\t$junc_number-$name_\t$num_observ\t$middle\t255,215,0\t$block_count_\t$block_sizes_\t$block_starts_\n!;
														}
														close(OUT_DETAILS);
													}
												}
												$name_2=$name_; 
											}
										}
									}
								}
								if ($outside==1) { $junction_status=""; }
							}
							if (($N==1) or ($name_2!~/^([^\-]+)\-([A-Z]+)/))
							{
								$name_2= "NO_GENE-N-$num_observ-$chromosome-$junc_number-$begin_intron-$end_intron"; 
							}
						}
					}

					my $new_line = ""; 
					if ($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.*)$/)
					{
						$new_line = "$1\t$2\t$3\t$name_2\t$5";
					}
					if ($junction_status!~/\w/) 
					{
						if ($threshold_N<=$num_observ)
						{
							print OUT_NOT qq!$new_line\n!;
						}
						else
						{
							print OUT_NOT_EXTRA qq!$new_line\n!;
						}
					}
					else
					{
						$stat{"$junction_status#$num_observ"}++;
						$count_junctions{"$junction_status"}++;
						if ($junction_status=~/^AN/) 
						{
							#my $line_=$line;
							#$line_=~s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+)/$1-$junction_name_mod/;
							if ($threshold_AN<=$num_observ)
							{
								print OUT_NOT qq!$new_line\n!;
							}
							else
							{
								print OUT_NOT_EXTRA qq!$new_line\n!;
							}
						}
					}
					$count_junctions++;
				}
			}
		}
		close(IN_RNA);
	}
	close(OUT);
	close(OUT_);
	close(OUT_AN);
	close(OUT_AN_);
	close(OUT_NOT);
	
	if (open(OUT, qq!>$filename1.filterA.stat.txt!))
	{
		print OUT qq!num_observ\tC_matches\tA_matches\tAN1_matches\tAN2_matches\n!;
		for(my $num_observ=1;$num_observ<=$num_observ_max;$num_observ++)
		{
			print OUT qq!$num_observ!;
			foreach my $type ("C","A","AN1","AN2")
			{
				if ($stat{"$type#$num_observ"}!~/\w/) { $stat{"$type#$num_observ"}=0;}
				print OUT qq!\t$stat{"$type#$num_observ"}!;
			}
			print OUT qq!\n!;
		}
		close(OUT);
	}
	print qq!Model exons: $count_model_exons\n!;
	print qq!RNA-Seq introns: $count_junctions\n!;
	foreach my $type ("C","A","AN1","AN2")
	{
		print qq!$type: $count_junctions{$type}\n!;
	}
}


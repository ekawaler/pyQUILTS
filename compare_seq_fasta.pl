#!/usr/local/bin/perl

use strict;

my $error=0;
my $filename1="";
my $filename2="";

if ($ARGV[0]=~/\w/) { $filename1=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $filename2=$ARGV[1];} else { $error=1; }

my $total_pep_count=0;
my %prot=();
my %desc=();
my %seq=();
my %proteins_count={};
my @stat=();
my $line="";

if ($error==0)
{
	if (open (LOG,">$filename1.compare_seq_fasta.log"))
	{
		my $bed_found=0;
		my $filename_=$filename1;
		$filename_=~s/\.fasta$//;
		my %bed=();
		if (open (IN,"$filename_.bed"))
		{
			while ($line=<IN>)
			{
				chomp($line);
				if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
				{
					my $name_=$4;
					$name_=~s/\s.*$//;
					#$name_=~s/\-.*$//;
					#if ($name_=~s/^gi\|[0-9]+\|ref\|(.*)\|$/$1/)
					if ($name_=~s/^ref\|(.*)\|$/$1/)
					{
						$name_=~s/\.[0-9]+$//;
					}
					$bed{$name_}=$line;
					$bed_found=1;
				} else { print LOG qq!Error parsing ($filename_.bed): $line!; }
			}
			close(IN);
		}
		my $filename=$filename1;
		my $i=1;
		if (open (IN,"$filename"))
		{
			#print qq!$filename\n!;
			my $name="";
			while ($line=<IN>)
			{
				chomp($line);
				if ($line=~/^>\s*(\S+)\s?(.*)$/)
				{
					my $name_=$1;
					my $description_=$2;
					$name_=~s/\s.*$//;
					#$name_=~s/\-.*$//;
					#if ($name_=~s/^gi\|[0-9]+\|ref\|(.*)\|$/$1/)
					if ($name_=~s/^ref\|(.*)\|$/$1/)
					{
						$name_=~s/\.[0-9]+$//;
					}
					if ($name=~/\w/ and $seq{"$i#$name"}=~/\w/)
					{
						$proteins_count{$filename}++;
						#print "$proteins_count. $name\n";
					}
					$name=$name_;
					$desc{"$i#$name"}=$description_;
					$seq{"$i#$name"}="";
					$prot{$name}=1;
				}
				else
				{
					$line=~s/\*\s*$//;
					$seq{"$i#$name"}.="$line";
				}
			}	
			if ($name=~/\w/ and $seq{"$i#$name"}=~/\w/)
			{
						$proteins_count{$filename}++;
						#print "$proteins_count. $name\n";
			}
			close(IN);
			print qq!$filename $proteins_count{$filename}\n!;
			print LOG qq!$filename $proteins_count{$filename}\n!;
		}
		my $filename=$filename2;
		my $i=2;
		if (open (IN,"$filename"))
		{
			#print qq!$filename\n!;
			my $name="";
			my $sequence="";
			while ($line=<IN>)
			{
				chomp($line);
				if ($line=~/^>\s*(\S+)\s?(.*)$/)
				{
					my $name_=$1;
					my $description_=$2;
					$name_=~s/\s.*$//;
					#$name_=~s/\-.*$//;
					#if ($name_=~s/^gi\|[0-9]+\|ref\|(.*)\|$/$1/)
					if ($name_=~s/^ref\|(.*)\|$/$1/)
					{
						$name_=~s/\.[0-9]+$//;
					}
					if ($name=~/\w/ and $seq{"$i#$name"}=~/\w/)
					{
						$proteins_count{$filename}++;
						#print "$proteins_count. $name\n";
					}
					$name=$name_;
					$desc{"$i#$name"}=$description_;
					$seq{"$i#$name"}="";
					$prot{$name}=1;
				}
				else
				{
					$line=~s/\*\s*$//;
					$seq{"$i#$name"}.="$line";
				}
			}	
			if ($name=~/\w/ and $seq{"$i#$name"}=~/\w/)
			{
						$proteins_count{$filename}++;
						#print "$proteins_count. $name\n";
			}
			close(IN);
			print qq!$filename $proteins_count{$filename}\n!;
			print LOG qq!$filename $proteins_count{$filename}\n!;
		}
		my $both=0;
		my $equal=0;
		my $star1=0;
		my $star2=0;
		my $longer=0;
		my $longer_fixed=0;
		my $shorter=0;
		my $different=0;
		if ($bed_found==1) { open (OUT_BED,">$filename_-corrected.bed"); }

		open (OUT1,">$filename1.extra.fasta");
		open (OUT2,">$filename2.extra.fasta");
		foreach my $name (keys %prot)
		{
			if ($seq{"1#$name"}=~/\w/ and $seq{"2#$name"}=~/\w/)
			{
				if ($seq{"1#$name"}=~/\*/)
				{
					print OUT_BED qq!$bed{$name}\n!;
					$star1++;
				}
				else
				{
					if ($seq{"2#$name"}=~/\*/)
					{
						print OUT_BED qq!$bed{$name}\n!;
						$star2++;
					}
					else
					{
						$both++;
						if ($seq{"1#$name"}=~/^$seq{"2#$name"}$/)
						{
							print OUT_BED qq!$bed{$name}\n!;
							$equal++;;
						}
						else 
						{
							$different++;
							if ($seq{"1#$name"}=~/^$seq{"2#$name"}/)
							{
								my $diff=substr $seq{"1#$name"},length($seq{"2#$name"});
								print LOG qq!\nDifferent (Longer: $diff->): $name\n1: $seq{"1#$name"}\n2: $seq{"2#$name"}\n!;
								$stat[length($diff)]++;
								$longer++;
								if (length($diff)==1)
								{
									if ($bed{$name}=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
									{
										my $chr_=$1;
										my $start_=$2;
										my $end_=$3;
										my $name_=$4;
										my $score_=$5;
										my $strand_=$6;
										my $copy_="$7\t$8\t$9";
										my $segment_count_=$10;
										my $segment_lengths_="$11,";
										my $segment_starts_="$12,";
										if ($strand_=~/\+/)
										{
											if ($segment_lengths_=~s/([0-9]+)\,$//) 
											{ 
												my $temp=$1; 
												$temp-=3; 
												$segment_lengths_.="$temp,"; 
												$end_-=3;
											}
										}
										else
										{
											if ($segment_lengths_=~s/^([0-9]+)\,//) 
											{ 
												my $temp=$1; 
												$temp-=3; 
												$segment_lengths_="$temp,$segment_lengths_,"; 
												$start_+=3;
												$temp=$segment_starts_;
												$segment_starts_="";
												if ($temp=~s/^([0-9]+)\,//) { $segment_starts_="$1,"; }
												while ($temp=~s/^([0-9]+)\,//) { my $temp=$1; $temp-=3; $segment_starts_.="$temp,"; }
											}
										}
										$segment_lengths_=~s/\,+$//;
										$segment_starts_=~s/\,+$//;
										print OUT_BED qq!$chr_\t$start_\t$end_\t$name_\t$score_\t$strand_\t$start_\t$end_\t0\t$segment_count_\t$segment_lengths_\t$segment_starts_\n!; #####
										$longer_fixed++;
									}
								}
							}
							else
							{
								if ($seq{"2#$name"}=~/^$seq{"1#$name"}/)
								{
									my $diff=substr $seq{"2#$name"},length($seq{"1#$name"});
									print LOG qq!\nDifferent (Shorter: ->$diff): $name\n1: $seq{"1#$name"}\n2: $seq{"2#$name"}\n!;
									$shorter++;
								}
								else
								{
									print LOG qq!\nDifferent: $name\n1: $seq{"1#$name"}\n2: $seq{"2#$name"}\n!;
								}
							}
						}
					}
				}
			}
			else
			{
				if ($seq{"1#$name"}!~/\w/)
				{
					print LOG qq!1: not found: $name\n!;
					print OUT2 qq!>$name $desc{"2#$name"}\n$seq{"2#$name"}\n!;
				}
				if ($seq{"2#$name"}!~/\w/)
				{
					print LOG qq!2: not found: $name\n!;
					print OUT1 qq!>$name $desc{"1#$name"}\n$seq{"1#$name"}\n!;
					print BED_OUT qq!$bed{$name}\n!;
				}
			}
		}
		close(OUT1);
		close(OUT2);
		print qq!Both: $both\n!;
		print qq!Equal: $equal\n!;
		print qq!Different: $different\n!;
		print qq!Longer: $longer (0: $stat[0], 1: $stat[1] 2: $stat[2], 3: $stat[3]) {!;
		print join /,/,@stat;
		print "}\n";
		print qq!Longer (BED file fixed): $longer_fixed\n!;
		print qq!Shorter: $shorter\n!;
		print qq!Stars: 1:$star1 2:$star2\n!;
		print LOG qq!Both: $both\n!;
		print LOG qq!Equal: $equal\n!;
		print LOG qq!Different: $different\n!;
		print LOG qq!Longer: $longer (0: $stat[0], 1: $stat[1] 2: $stat[2], 3: $stat[3]) {!;
		print LOG join /,/,@stat;
		print LOG "}\n";
		print LOG qq!Longer (BED file fixed): $longer_fixed\n!;
		print LOG qq!Shorter: $shorter\n!;
		print LOG qq!Stars: 1:$star1 2:$star2\n!;
		close(LOG);
		if ($bed_found==1) { close(OUT_BED); }
	}
}
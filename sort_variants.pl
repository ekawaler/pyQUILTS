#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $filename="";
my $filename_variants="";

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $filename_variants=$ARGV[1];} else { $error=1; }

my $filename__=$filename;
if ($filename__!~s/\.bed.dna$//) { $error=1; }
my $filename_=$filename_variants;
if ($filename_!~s/\.var$//) { $error=1; }

my %mapping = (	"TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
				"CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
				"ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
				"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
				
				"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
				"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
				"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
				"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
				
				"TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
				"CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
				"AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
				"GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
				
				"TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
				"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
				"AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
				"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G");
				
if ($error==0)
{
	my $line="";
	my %variants=();
	if (open (IN,"$filename_variants"))
	{
		print qq!$filename_variants\n!;
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]+)\t([^\t]*)\t([^\t]*)/)
			{
				my $name=$1;
				my $var=$3;
				$var=~s/[\r\n]+$//;
				$variants{$name}.=$var;
				print qq!$name $variants{$name}\n!;
			}
		}
		close(IN);
	}
	if (open (IN,"$filename"))
	{
		open (OUT,">$filename_.aa.var");
		open (OUT_BED,">$filename__.aa.var.bed");
		open (OUT_BED_DNA,">$filename__.aa.var.bed.dna");
		open (OUT_NO,">$filename_.noaa.var");
		open (OUT_STOP,">$filename_.stop.var");
		open (OUT_STOP_BED,">$filename__.stop.var.bed");
		open (OUT_STOP_BED_DNA,">$filename__.stop.var.bed.dna");
		open (OUT_STOP_,">$filename_.stop-removed.var");
		open (OUT_STOP_BED_,">$filename__.stop-removed.var.bed");
		open (OUT_STOP_BED_DNA_,">$filename__.stop-removed.var.bed.dna");
		open (OUT_INDEL,">$filename_.indel.var");
		open (OUT_INDEL_BED,">$filename__.indel.var.bed");
		open (OUT_INDEL_BED_DNA,">$filename__.indel.var.bed.dna");
		open (OUT_ERROR,">$filename_.error.var");
		my $proteins_count=0;
		my $name="";
		my $description="";
		my $bed="";
		my $map="";
		my $map_variants="";
		my $var_protein="";
		my $seq="";
		my $seq_before="";
		my $seq_after="";
		my $seq_line_count=0;
		my $seq_line_count_=0;
		my @seq_line=();
		my $seq_line_leftover="";
		my $seq_line_before="";
		my $seq_line_after="";
		my $chr="";
		my $start="";
		my $end="";
		my $strand="";
		my $lengths="";
		my $offsets="";
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n\r]*)/)
			{
				$chr=$1;
				$start=$2;
				$end=$3;
				$name=$4;
				$strand=$6;
				$lengths=$11;
				$offsets=$12;
				$description="";
				$bed=$line;
				$seq="";
				$map="";
				$map_variants="";
				$var_protein="";
				$seq_before="";
				$seq_after="";
				$seq_line_count=0;
				$seq_line_count_=0;
				@seq_line=();
				$seq_line_leftover="";
				$seq_line_before="";
				$seq_line_after="";
			} 
			else 
			{
				if ($line=~/^>(\S+)\s+(.*)\s*\(MAP:/)
				{
					$description=$2;
					if ($line=~/\(MAP:\s*([^\)]+)\)/)
					{
						$map="(MAP:$1)";
						if ($map=~/^\S+\s\S+\s\S+\s(\S+)\)$/)
						{
							$map_variants=$1;
						}
					}
					if ($line=~/\(VAR:\s*([^\)]+)\)/)
					{
						$var_protein=$1;
					}
					if ($description=~/\w/ and $description!~/\s$/) { $description.=" "; }
				} 
				else 
				{
					if ($line=~/^([^\t]*)\t([0-9]+)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
					{
						$seq.="\U$8";
						my $line_number=$2;
						$seq_line[$line_number]=$line;
						if ($line_number!=$seq_line_count)
						{
							print OUT_ERROR qq!$line_number\!=$seq_line_count\n!; 
						}
						$seq_line_count++;
					} 
					else 
					{ 
						if ($line=~/^([^\t]*)\t(\-1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
						{
							$seq_before.="\U$8";
							$seq_line_before=$line;
						} 
						else 
						{ 
							if ($line=~/^([^\t]*)\t(\+1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
							{
								$seq_after.="\U$8";
								$seq_line_after=$line;
								my $protein="";
								my $var_aa="";
								my $var_noaa="";
								my $var_stop="";
								my $var_stop_="";
								my $var_indel="";
								my $var_aa_string="";
								my $var_noaa_string="";
								my $var_stop_string="";
								my $var_stop_string_="";
								my $var_indel_string="";
								my $temp=$variants{$name};
								my $temp_="";
								while ($temp=~s/^([GS])\-([A-Z\.]+),?([A-Z\.]*)([0-9]+)([A-Z\.\-]+),?([A-Z\.\-]*):([0-9\.edED\+\-]+),//)
								{
									my $type=$1;
									my $old=$2;
									my $old_extra=$3;
									my $pos=$4;
									my $new=$5;
									my $new_extra=$6;
									my $qual=$7;
									if ($old_extra=~/\w/ or $new_extra=~/\w/) { print OUT_ERROR qq!Ignored: $line: $type-$old(,$old_extra)$pos$new(,$new_extra):$qual\n!; }
									if ($type=~/^S$/)
									{
										if ($temp_!~/([GS])\-[A-Z\.]+$pos([A-Z\.]+)\:([0-9\-\+\.edED]+),/)
										{
											$temp_.="$type-$old$pos$new:$qual,";
										}
									}
									else
									{
										if ($temp_!~s/S\-[A-Z\.]+$pos([A-Z\.]+)\:([0-9\-\+\.edED]+),/$type-$old$pos$new:$qual,/)
										{
											if ($temp_!~/([GS])\-[A-Z\.]+$pos([A-Z\.]+)\:([0-9\-\+\.edED]+),/)
											{
												$temp_.="$type-$old$pos$new:$qual,";
											}
										}
									}
								}
								if ($temp=~/\w/) { print OUT_ERROR qq!Error: $line: $temp\n!; }
								$variants{$name}=$temp_;
								$temp=$variants{$name};
								my $seq_new=$seq;
								$temp_="";
								while ($temp=~s/^([GS])\-([A-Z\.\-]+)([0-9]+)([A-Z\.\-]+):([0-9\.edED\+\-]+),//)
								{
									my $type=$1;
									my $old=$2;
									my $pos=$3;
									my $new=$4;
									my $qual=$5;
									my $old_ = substr $seq,$pos,length($old);
									if (not($old=~/\./ or $old=~/^$old_$/))
									{
										print OUT_ERROR qq!$name $variants{$name} $pos:$old\!=$old_\n!; 
									}
									else
									{
										if(length($new)!=length($old) or $old=~/\./ or $new=~/\./ or $old=~/\-/ or $new=~/\-/)
										{
											$var_indel.="$type-$old$pos$new:$qual,";
										}
										else
										{
											substr $seq_new,$pos,length($old),$new;
											$temp_.="$type-$old$pos$new:$qual,";
										}
									}
								}
								if ($temp=~/\w/) { print OUT_ERROR qq!Error: $line: $temp\n!; }
								$variants{$name}=$temp_;
								$temp=$variants{$name};
								while ($temp=~s/^([GS])\-([A-Z]+)([0-9]+)([A-Z]+):([0-9\.edED\+\-]+),//)
								{
									my $type=$1;
									my $old=$2;
									my $pos=$3;
									my $new=$4;
									my $qual=$5;
									my $stop=0;
									my $stop_=0;
									my $noaa=1;
									if ($strand=~/\-/)
									{
										my $seq_rev = reverse $seq;
										$seq_rev=~tr/ATCG/TAGC/;
										my $seq_rev_new = reverse $seq_new;
										$seq_rev_new=~tr/ATCG/TAGC/;
										my $pos_=length($seq_rev_new)-1-$pos;
										for(my $pos_aa=int($pos_/3);$pos_aa<=int($pos_/3)+int(length($new)/3);$pos_aa++)
										{
											my $pos_aa_=$pos_aa+1;
											my $triplet_old = substr($seq_rev, 3*$pos_aa, 3);
											my $triplet = substr($seq_rev_new, 3*$pos_aa, 3);
											#print qq!$type-$old$pos$new:$qual: $triplet_old->$triplet $pos_aa:$mapping{$triplet_old}->$mapping{$triplet}\n!;
											if ($mapping{$triplet_old}!~/^\*$/ and $mapping{$triplet}=~/^\*$/) 
											{
												$stop=1; 
												if ($var_stop_string!~/$mapping{$triplet_old}$pos_aa_\*\:/)
												{
													$var_stop_string.="$type\-$mapping{$triplet_old}$pos_aa_$mapping{$triplet}:$qual,";
												}
											}
											else
											{
												if ($mapping{$triplet_old}=~/^\*$/ and $mapping{$triplet}!~/^\*$/) 
												{ 
													$stop_=1;
													if ($var_stop_string_!~/\*$pos_aa_$mapping{$triplet}\:/)
													{
														$var_stop_string_.="$type\-$mapping{$triplet_old}$pos_aa_$mapping{$triplet}:$qual,";
													}
												}
												else
												{
													if ($mapping{$triplet_old} ne $mapping{$triplet}) 
													{
														$noaa=0;
														if ($var_aa_string!~/$mapping{$triplet_old}$pos_aa_$mapping{$triplet}\:/)
														{
															$var_aa_string.="$type\-$mapping{$triplet_old}$pos_aa_$mapping{$triplet}:$qual,";
														}
													}
												}
											}
										}
									}
									else
									{
										for(my $pos_aa=int($pos/3);$pos_aa<=int($pos/3)+int(length($new)/3);$pos_aa++)
										{
											my $pos_aa_=$pos_aa+1;
											my $triplet_old = substr($seq, 3*$pos_aa, 3);
											my $triplet = substr($seq_new, 3*$pos_aa, 3);
											#print qq!$type-$old$pos$new:$qual: $triplet_old->$triplet $pos_aa:$mapping{$triplet_old}->$mapping{$triplet}\n!;
											if ($mapping{$triplet_old}!~/^\*$/ and $mapping{$triplet}=~/^\*$/) 
											{
												$stop=1; 
												if ($var_stop_string!~/$mapping{$triplet_old}$pos_aa_\*\:/)
												{
													$var_stop_string.="$type-$mapping{$triplet_old}$pos_aa_$mapping{$triplet}:$qual,";
												}
											}
											else
											{
												if ($mapping{$triplet_old}=~/^\*$/ and $mapping{$triplet}!~/^\*$/) 
												{ 
													$stop_=1;
													if ($var_stop_string_!~/\*$pos_aa_$mapping{$triplet}\:/)
													{
														$var_stop_string_.="$type-$mapping{$triplet_old}$pos_aa_$mapping{$triplet}:$qual,";
													}
												}
												else
												{
													if ($mapping{$triplet_old} ne $mapping{$triplet}) 
													{
														$noaa=0;
														if ($var_aa_string!~/$mapping{$triplet_old}$pos_aa_$mapping{$triplet}\:/)
														{
															$var_aa_string.="$type-$mapping{$triplet_old}$pos_aa_$mapping{$triplet}:$qual,";
														}
													}
												}
											}
										}
									}
									if ($stop==1)
									{
										$var_stop.="$type-$old$pos$new:$qual,";
									}
									else
									{
										if ($stop_==1)
										{
											$var_stop_.="$type-$old$pos$new:$qual,";
										}
										else
										{
											if ($noaa==1)
											{
												$var_noaa.="$type-$old$pos$new:$qual,";
											}
											else
											{
												$var_aa.="$type-$old$pos$new:$qual,";
											}
										}
									}
								}
								if ($temp=~/\w/) { print OUT_ERROR qq!Error: $line: $temp\n!; }
								if ($var_aa=~/\w/)
								{
									my @seq_line_=@seq_line;
									my $name__=$name;
									$name__.="-variant";
									print OUT qq!$name\t$var_aa\n!;
									$map="(MAP:$chr:$start$strand $lengths $offsets $map_variants";
									my $temp=$var_aa;
									while ($temp=~s/^([GS])\-([A-Z]+)([0-9]+)([A-Z]+):([0-9\.edED\+\-]+),//)
									{
										my $type=$1;
										my $old=$2;
										my $pos=$3;
										my $new=$4;
										my $qual=$5;
										$map.="$type-$old$pos$new:$qual,";
										my $sum=0;
										for(my $k=0;$k<$seq_line_count;$k++)
										{
											if ($seq_line_[$k]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/)
											{
												my $length=$5;
												if ($sum<=$pos and $pos<$sum+$length)
												{
													#print qq!---$type-$old$pos$new:$qual\n!;
													#print qq!1. $seq_line_[$k]\n!;
													substr $seq_line_[$k],length($seq_line_[$k])-$length+$pos-$sum,length($old),$new;
													#print qq!2. $seq_line_[$k]\n!;
												}
												$sum+=$length;
											}
										}
									}
									if ($temp=~/\w/) { print OUT_ERROR qq!Error: $var_aa: $temp\n!; }
									$map.=")";
									if ($var_protein=~/\w/ or $var_aa_string=~/\w/) { $map.=" (VAR:$var_protein$var_aa_string)"; }
									$bed=~s/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$2\t$3\t$name__\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12/;
									print OUT_BED qq!$bed\n!;
									print OUT_BED_DNA qq!$bed\n!;
									print OUT_BED_DNA qq!>$name__ $description$map\n!;
									print OUT_BED_DNA qq!$seq_line_before\n!;
									my $sum=0;
									for(my $k=0;$k<$seq_line_count;$k++)
									{
										print OUT_BED_DNA qq!$seq_line_[$k]\n!;
									}
									print OUT_BED_DNA qq!$seq_line_after\n!;
								}
								if ($var_noaa=~/\w/)
								{
									print OUT_NO qq!$name\t$var_noaa\n!;
								}
								if ($var_stop=~/\w/)
								{
									my @seq_line_=@seq_line;
									my $name__=$name;
									$name__.="-stop-introduced";
									print OUT_STOP qq!$name\t$var_stop\n!;
									$map="(MAP:$chr:$start$strand $lengths $offsets $map_variants";
									my $temp=$var_stop;
									while ($temp=~s/^([GS])\-([A-Z]+)([0-9]+)([A-Z]+):([0-9\.edED\+\-]+),//)
									{
										my $type=$1;
										my $old=$2;
										my $pos=$3;
										my $new=$4;
										my $qual=$5;
										$map.="$type-$old$pos$new:$qual,";
										my $sum=0;
										for(my $k=0;$k<$seq_line_count;$k++)
										{
											if ($seq_line_[$k]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/)
											{
												my $length=$5;
												if ($sum<=$pos and $pos<$sum+$length)
												{
													#print qq!---$type-$old$pos$new:$qual\n!;
													#print qq!1. $seq_line_[$k]\n!;
													substr $seq_line_[$k],length($seq_line_[$k])-$length+$pos-$sum,length($old),$new;
													#print qq!2. $seq_line_[$k]\n!;
												}
												$sum+=$length;
											}
										}
									}
									if ($temp=~/\w/) { print OUT_ERROR qq!Error: $var_stop: $temp\n!; }
									$map.=")";
									if ($var_protein=~/\w/ or $var_stop_string=~/\w/) { $map.=" (VAR:$var_protein$var_stop_string)"; }
									$bed=~s/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$2\t$3\t$name__\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12/;
									print OUT_STOP_BED qq!$bed\n!;
									print OUT_STOP_BED_DNA qq!$bed\n!;
									print OUT_STOP_BED_DNA qq!>$name__ $description$map\n!;
									print OUT_STOP_BED_DNA qq!$seq_line_before\n!;
									my $sum=0;
									for(my $k=0;$k<$seq_line_count;$k++)
									{
										print OUT_STOP_BED_DNA qq!$seq_line_[$k]\n!;
									}
									print OUT_STOP_BED_DNA qq!$seq_line_after\n!;
								}
								if ($var_stop_=~/\w/)
								{
									my @seq_line_=@seq_line;
									my $name__=$name;
									$name__.="-stop-removed";
									print OUT_STOP_ qq!$name\t$var_stop_\n!;
									$map="(MAP:$chr:$start$strand $lengths $offsets $map_variants";
									my $temp=$var_stop_;
									while ($temp=~s/^([GS])\-([A-Z]+)([0-9]+)([A-Z]+):([0-9\.edED\+\-]+),//)
									{
										my $type=$1;
										my $old=$2;
										my $pos=$3;
										my $new=$4;
										my $qual=$5;
										$map.="$type-$old$pos$new:$qual,";
										my $sum=0;
										for(my $k=0;$k<$seq_line_count;$k++)
										{
											if ($seq_line_[$k]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/)
											{
												my $length=$5;
												if ($sum<=$pos and $pos<$sum+$length)
												{
													#print qq!---$type-$old$pos$new:$qual\n!;
													#print qq!1. $seq_line_[$k]\n!;
													substr $seq_line_[$k],length($seq_line_[$k])-$length+$pos-$sum,length($old),$new;
													#print qq!2. $seq_line_[$k]\n!;
												}
												$sum+=$length;
											}
										}
									}
									if ($temp=~/\w/) { print OUT_ERROR qq!Error: $var_stop_: $temp\n!; }
									$map.=")";
									if ($var_protein=~/\w/ or $var_stop_string_=~/\w/) { $map.=" (VAR:$var_protein$var_stop_string_)"; }
									$bed=~s/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$2\t$3\t$name__\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12/;
									print OUT_STOP_BED_ qq!$bed\n!;
									print OUT_STOP_BED_DNA_ qq!$bed\n!;
									print OUT_STOP_BED_DNA_ qq!>$name__ $description$map\n!;
									print OUT_STOP_BED_DNA_ qq!$seq_line_before\n!;
									my $sum=0;
									for(my $k=0;$k<$seq_line_count;$k++)
									{
										print OUT_STOP_BED_DNA_ qq!$seq_line_[$k]\n!;
									}
									print OUT_STOP_BED_DNA_ qq!$seq_line_after\n!;
								}
								if ($var_indel=~/\w/)
								{
									print OUT_INDEL qq!$name\t$var_indel\n!;
									my @seq_line_=@seq_line;
									my $name__=$name;
									$name__.="-indel";
									$map="(MAP:$chr:$start$strand $lengths $offsets $map_variants";
									my $temp=$var_indel;
									while ($temp=~s/^([GS])\-([A-Z\.\-]+)([0-9]+)([A-Z\.\-]+):([0-9\.edED\+\-]+),//)
									{
										my $type=$1;
										my $old=$2;
										my $pos=$3;
										my $new=$4;
										my $qual=$5;
										$map.="$type-$old$pos$new:$qual,";
										if ($old=~/^\.$/) { $old=""; }
										if ($new=~/^\.$/) { $new=""; }
										my $sum=0;
										for(my $k=0;$k<$seq_line_count;$k++)
										{
											if ($seq_line_[$k]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/)
											{
												my $start=$4;
												my $length_original=$5;
												my $length_actual=length($8);
												if ($sum<=$pos and $pos<$sum+$length_original)
												{
													#print qq!---$type-$old$pos$new:$qual\n!;
													#print qq!1. $seq_line_[$k]\n!;
													my $seq_line_new=substr $seq_line_[$k],0,length($seq_line_[$k])-$length_actual;
													$seq_line_new.=substr $seq_line_[$k],length($seq_line_[$k])-$length_actual+$pos-$sum+length($old);
													my $seq_line_new_start=$start+$pos-$sum+length($old);
													my $seq_line_new_length=length(substr $seq_line_[$k],length($seq_line_[$k])-$length_actual+$pos-$sum+length($old));
													$seq_line_new=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$2\t$3\t$seq_line_new_start\t$seq_line_new_length\t$6\t$7\t$8/;
													$seq_line_[$k]=substr $seq_line_[$k],0,length($seq_line_[$k])-$length_actual+$pos-$sum+length($old);
													substr $seq_line_[$k],length($seq_line_[$k])-length($old),length($old),$new;
													my $seq_line_old_length=$length_original-$seq_line_new_length;
													$seq_line_[$k]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$2\t$3\t$4\t$seq_line_old_length\t$6\t$7\t$8/;
													#print qq!2. $seq_line_[$k]\n!;
													#print qq!3. $seq_line_new\n!;
													for(my $l=$seq_line_count-1;$l>$k;$l--)
													{
														my $l_=$l+1;
														$seq_line_[$l_]=$seq_line_[$l];
														$seq_line_[$l_]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$l_\t$3\t$4\t$5\t$6\t$7\t$8/;
													}
													my $l_=$k+1;
													$seq_line_[$l_]=$seq_line_new;
													$seq_line_[$l_]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$l_\t$3\t$4\t$5\t$6\t$7\t$8/;
													$seq_line_count++;
												}
												$sum+=$length_original;
											}
										}
									}
									$map.=")";
									if ($var_protein=~/\w/ or $var_indel_string=~/\w/) { $map.=" (VAR:$var_protein$var_indel_string)"; }
									$bed=~s/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$/$1\t$2\t$3\t$name__\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12/;
									print OUT_INDEL_BED qq!$bed\n!;
									print OUT_INDEL_BED_DNA qq!$bed\n!;
									print OUT_INDEL_BED_DNA qq!>$name__ $description$map\n!;
									print OUT_INDEL_BED_DNA qq!$seq_line_before\n!;
									my $sum=0;
									for(my $k=0;$k<$seq_line_count;$k++)
									{
										print OUT_INDEL_BED_DNA qq!$seq_line_[$k]\n!;
									}
									print OUT_INDEL_BED_DNA qq!$seq_line_after\n!;
								}
							} 
							else 
							{ 
								print OUT_ERROR qq!Error parsing: $line\n!; 
							}
						}
					}
				}
			}
		}
		close(IN);
	}	
	close(OUT);
	close(OUT_STOP);
	close(OUT_STOP_);
	close(OUT_INDEL);
	close(OUT_BED);
	close(OUT_STOP_BED);
	close(OUT_STOP_BED_);
	close(OUT_INDEL_BED);
	close(OUT_BED_DNA);
	close(OUT_STOP_BED_DNA);
	close(OUT_STOP_BED_DNA_);
	close(OUT_INDEL_BED_DNA);
	close(OUT_ERROR);
}
#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $filename="";
my $continue="";

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $continue=$ARGV[1];} else { $continue=0; }

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
	my $filename_=$filename;
	$filename_=~s/\.bed\.dna$//;
	my $line="";
	my $line_number=0;
	my %descriptions=();
	if (open (IN,"$filename_-descriptions.txt"))
	{
		while ($line=<IN>)
		{
			chomp($line); $line=~s/([\n\r]+)$//;
			if ($line=~/^([^\t]+)\t([^\t]+)/)
			{
				my $name=$1;
				my $desc=$2;
				$name=~s/^([^\-]+)\-.*$/$1/;
				$descriptions{$name}=$desc;
			}
		}
		close(IN);
	}
	my %genes=();
	if (open (IN,"$filename_-genes.txt"))
	{
		while ($line=<IN>)
		{
			chomp($line); $line=~s/([\n\r]+)$//;
			if ($line=~/^([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]+)/)
			{
				my $name=$1;
				my $gene=$4;
				$name=~s/^([^\-]+)\-.*$/$1/;
				$genes{$name}=$gene;
			}
		}
		close(IN);
	}
	if (open (IN,"$filename"))
	{
		open (OUT,">$filename.fasta");
		open (OUT_BED,">$filename.fasta.bed");
		open (OUT_BED_DNA,">$filename.fasta.bed.dna");
		my $proteins_count=0;
		my $name="";
		my $name_old="";
		my $bed="";
		my $map="";
		my $map_variants="";
		my $var_protein="";
		my $description="";
		my $seq="";
		my $seq_length="";
		my $seq_before="";
		my $seq_after="";
		my $seq_line_count=0;
		my $seq_line_count_=0;
		my @seq_line=();
		my @seq_line_=();
		my $seq_line_leftover="";
		my $seq_line_before="";
		my $seq_line_after="";
		my $chr="";
		my $start="";
		my $end="";
		my $strand="";
		my $block=0;
		while ($line=<IN>)
		{
			chomp($line); $line=~s/([\n\r]+)$//;
			$line_number++;
			if ($line=~/^(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
			{
				$name_old=$name;
				$chr=$1;
				$start=$2;
				$end=$3;
				$name=$4;
				$strand=$6;
				$bed=$line;
				$map="(MAP:$chr:$start$strand $11 $12)";
				$map_variants="";
				$var_protein="";
				$description="";
				$seq="";
				$seq_length="";
				$seq_before="";
				$seq_after="";
				$seq_line_count=0;
				$seq_line_count_=0;
				@seq_line=();
				@seq_line_=();
				$seq_line_leftover="";
				$seq_line_before="";
				$seq_line_after="";
				if ($block!=0) { print qq!translate.pl:$line_number:$filename: Warning: end of block not found: $name_old\n!;}
				$block=1;
			} 
			else 
			{
				if ($line=~/^>(\S+) (.*) *(\(MAP:.*\))/)
				{
					$description=$2;
					if ($line=~/\(MAP:\s*([^\)]+)\)/)
					{
						$map="(MAP:$1)";
						if ($map=~/^\S+\s\S+\s\S+\s(\S+)\)$/)
						{
							$map_variants=$1;
						}
						if ($map_variants=~/\w/) { $map_variants=" $map_variants"; }
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
						$seq_line[$line_number]="$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8";
						#print qq!$seq_line[$line_number]\n!;
						if ($line_number!=$seq_line_count)
						{
							print qq!$line_number\!=$seq_line_count\n!; 
						}
						$seq_line_count++;
					} 
					else 
					{ 
						if ($line=~/^([^\t]*)\t(\-1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
						{
							$seq_before="\U$8";
							$seq_line_before="$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8";
						} 
						else 
						{ 
							if ($line=~/^([^\t]*)\t(\+1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
							{
								$block=0;
								$seq_after="\U$8";
								$seq_line_after="$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8";
								$seq_length=length($seq);
								#print qq!$seq_line_after\n!;
								if ($seq!~/Error Reading Sequence/i)
								{
									if ($continue==1)
									{
										if ($strand=~/\-/)
										{
											$seq="$seq_before$seq";
											$seq_before="";
											if ($seq_line_before=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/)
											{
												my $start_before_=$4;
												my $length_before_=$5;
												my $seq_before_=$8;
												if ($seq_line[0]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/)
												{
													my $length_=$5+$length_before_;
													my $seq_="$seq_before_$8";
													$seq_line[0]="$1\t$2\t$3\t$4\t$length_\t$6\t$7\t$seq_";
													$seq_line_before="$1\t-1\t$3\t$start_before_\t0\t$6\t0\t";
													#print "-0-$map\n";
													if ($map=~/^\(MAP:chr[0-9XYMxym]+\:[0-9]+[\+\-] (\S+) (\S+)(.*)\)/)
													{
														my $lengths="$1";
														my $offsets_="$2,";
														$start-=$length_before_;
														$lengths=~s/^([0-9]+)/$length_/;
														my $offsets="";
														while($offsets_=~s/^([0-9]+),//)
														{
															my $offset_=$1;
															$offset_+=$length_before_;
															if ($offsets=~/\w/)
															{
																$offsets.="$offset_,";
															}
															else
															{
																$offsets.="0,";
															}
														}
														$offsets=~s/\,$//;
														$map="(MAP:$chr:$start$strand $lengths $offsets$map_variants)";
														#print "---$map\n";
													}
												}
											}
										}
										else
										{
											$seq.=$seq_after;
											$seq_after="";
											if ($seq_line_after=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/)
											{
												my $start_after_=$4;
												my $length_after_=$5;
												my $seq_after_=$8;
												if ($seq_line[$seq_line_count-1]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/)
												{
													my $length_=$5+$length_after_;
													my $seq_="$8$seq_after_";
													$seq_line[$seq_line_count-1]="$1\t$2\t$3\t$4\t$length_\t$6\t$7\t$seq_";
													$start_after_+=$length_after_;
													$seq_line_after="$1\t+1\t$3\t$start_after_\t0\t$6\t0\t";
													#print "-0-$map\n$seq_line[$seq_line_count-1]\n";
													if ($map=~/^\(MAP:chr[0-9XYMxym]+\:[0-9]+[\+\-] (\S+) (\S+)(.*)\)$/)
													{
														my $lengths="$1";
														my $offsets="$2";
														$lengths=~s/,([0-9]+)$/,$length_/;
														$map="(MAP:$chr:$start$strand $lengths $offsets$map_variants)";
														#print "---$map\n";
													}
												}
											}
										}
									}
									my $protein="";
									if ($strand=~/\-/)
									{
										$seq = reverse $seq;
										$seq=~tr/ATCG/TAGC/;
									}
									for(my $n=0;$n<length($seq);$n=$n+3)
									{
										my $triplet = substr($seq, $n, 3);
										if (length($triplet)==3)
										{
											if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
											$protein.=$mapping{$triplet}; 
										}
									}
									my $name_=$name;
									$name_=~s/^([^\-]+)\-.*$/$1/;
									if ($description!~/\w/) 
									{ 
										$description="$descriptions{$name_}"; 
										if ($genes{$name_}=~/\w/) 
										{
											if ($description=~/\w/ and $description!~/\s$/) { $description.=" "; }
											$description.="GN=$genes{$name_}"; 
										}
									}
									if ($description=~/\w/ and $description!~/\s$/) { $description.=" "; }
									my $name__=$name;
									my $stop = index($protein, '*')+1;
									my $protein_=$protein;
									my $extra="";
									if ($continue==1 and $protein!~/\*/) { print qq!translate.pl:$line_number:$filename: Warning: $name__ stop not found\n!; }
									if ($protein=~s/\*(.*)$//)
									{
										$extra=$1;
										if ($extra=~/\w/)
										{
											if ($name__!~/\-stop$stop$/ and $name__!~/\-stop$stop\-/) 
											{ 
												$name__.="-stop$stop";
												my $temp=$map_variants;
												$map_variants="";
												while($temp=~s/^\s*([GS])\-([A-Z\*]+)([0-9]+)([A-Z\*]+)\:([^\,]+)\,//)
												{
													my $type=$1;
													my $old=$2;
													my $pos=$3;
													my $new=$4;
													my $qual=$5;
													if ($strand=~/\+/)
													{
														if ($pos<=$stop*3) { $map_variants.="$type-$old$pos$new:$qual,"; }
													}
													else
													{
														if ($seq_length-$pos<=$stop*3) { $map_variants.="$type-$old$pos$new:$qual,"; }
													}
												}
												if ($map_variants=~/\w/) { $map_variants=" $map_variants"; }
												my $temp=$var_protein;
												$var_protein="";
												while($temp=~s/^\s*([GS])\-([A-Z\*]+)([0-9]+)([A-Z\*]+)\:([^\,]+)\,//)
												{
													my $type=$1;
													my $old=$2;
													my $pos=$3;
													my $new=$4;
													my $qual=$5;
													if ($pos<=$stop) { $var_protein.="$type-$old$pos$new:$qual,"; }
												}
											}
											if ($map=~/^\(MAP:chr[0-9XYMxym]+\:[0-9]+[\+\-] (\S+) (\S+)(.*)\)$/)
											{
												my $lengths_="$1,";
												my $offsets_="$2,";
												my $lengths="";
												my $offsets="";
												my $sum=0;
												if ($strand=~/\+/)
												{
													$seq_line_before=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t$2\t$3\t$4\t$5\t$6\t$7\t$8/;
													my $k=0;
													while($lengths_=~s/^([^\,]+)\,//)
													{
														my $length=$1;
														if($offsets_=~s/^([^\,]+)\,//)
														{
															my $offset=$1;
															#print qq!$offset $length $sum $stop\n!;
															if ($sum+$length<3*$stop)
															{
																#print qq!1\n!;
																$lengths.="$length,";
																$offsets.="$offset,";
																$end=$start+$offset+$length;
																$seq_line_[$seq_line_count_]=$seq_line[$k];
																$seq_line_[$seq_line_count_]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t$2\t$3\t$4\t$5\t$6\t$7\t$8/;
																$seq_line_count_++;
															}
															else
															{
																if ($sum<3*$stop)
																{
																	#print qq!2\n!;
																	my $length_=3*$stop-$sum;
																	my $length__=$length-$length_;
																	$lengths.="$length_,";
																	$offsets.="$offset,";
																	$end=$start+$offset+$length_;
																	$seq_line_[$seq_line_count_]=substr $seq_line[$k],0,length($seq_line[$k])-$length__;
																	if ($length__>0)
																	{
																		$seq_line_leftover=substr $seq_line[$k],0,length($seq_line[$k])-$length;
																		$seq_line_leftover.=substr $seq_line[$k],length($seq_line[$k])-$length__,$length__;
																		$seq_line_leftover=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t+1\t$3\t$4\t$length__\t$6\t$7\t$8/;
																	}
																	$seq_line_[$seq_line_count_]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t$2\t$3\t$4\t$length_\t$6\t$7\t$8/;
																	$seq_line_count_++;
																}
															}
															$sum+=$length;
														}
														$k++;
													}
												} 
												else
												{
													$seq_line_after=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t$2\t$3\t$4\t$5\t$6\t$7\t$8/;
													my $k = ($lengths_=~tr/,/,/)-1;
													while($lengths_=~s/([^\,]+)\,$//)
													{
														my $length=$1;
														if($offsets_=~s/([^\,]+)\,$//)
														{
															my $offset=$1;
															#print qq!$offset $length $sum $stop\n!;
															if ($sum+$length<3*$stop)
															{
																#print qq!1.$seq_line_count_.$k\n!;
																$lengths="$length,$lengths";
																$offsets="$offset,$offsets";
																$seq_line_[$seq_line_count_]=$seq_line[$k];
																$seq_line_[$seq_line_count_]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t$2\t$3\t$4\t$5\t$6\t$7\t$8/;
																$seq_line_count_++;
															}
															else
															{
																if ($sum<3*$stop)
																{
																	#print qq!2.$seq_line_count_.$k\n!;
																	my $length_=3*$stop-$sum;
																	my $length__=$length-$length_;
																	my $offset_=$offset+$length__;
																	$lengths="$length_,$lengths";
																	$offsets="$offset_,$offsets";
																	my $offset_abs=$offset_+$start;
																	$seq_line_[$seq_line_count_]=substr $seq_line[$k],0,length($seq_line[$k])-$length;
																	$seq_line_[$seq_line_count_].=substr $seq_line[$k],length($seq_line[$k])-$length_,$length_;
																	$seq_line_[$seq_line_count_]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t$2\t$3\t$offset_abs\t$length_\t$6\t$7\t$8/;
																	if ($length__>0)
																	{
																		$offset_abs=$offset_+$start-$length__;
																		$seq_line_leftover=substr $seq_line[$k],0,length($seq_line[$k])-$length_;
																		$seq_line_leftover=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$name__\t-1\t$3\t$offset_abs\t$length__\t$6\t$7\t$8/;
																	}
																	$seq_line_count_++;
																}
															}
															$sum+=$length;
														}
														$k--;
													}
													$offsets_=$offsets;
													$offsets="0,";
													if ($offsets_=~s/^([^\,]+)\,//)
													{
														my $offset0=$1;
														$start+=$offset0;
														while ($offsets_=~s/^([^\,]+)\,//)
														{
															my $offset=$1-$offset0;
															$offsets.="$offset,";
														}
													}
												}
												$lengths=~s/\,+$//;
												$offsets=~s/\,+$//;
												$map="(MAP:$chr:$start$strand $lengths $offsets$map_variants)";
												if ($bed=~s/([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$1\t$start\t$end\t$name__\t$5\t$strand\t$start\t$end\t$9\t$seq_line_count_\t$lengths\t$offsets/)
												{
													print OUT_BED_DNA qq!$bed\n!;
													print OUT_BED_DNA qq!>$name__ $description$map!;
													if ($var_protein=~/\w/) { print OUT_BED_DNA " (VAR:$var_protein)"; }
													print OUT_BED_DNA qq!\n!;
													if ($strand=~/\+/)
													{
														if ($seq_line_before=~/\w/) { print OUT_BED_DNA qq!$seq_line_before\n!; }
														my $k=0;
														for(;$k<$seq_line_count_;$k++)
														{
															print OUT_BED_DNA qq!$seq_line_[$k]\n!;
														}
														if ($seq_line_leftover!~/\w/)
														{
															if ($seq_line_[$seq_line_count_-1]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/)
															{
																my $start_=$4;
																my $length_=$5;
																$start_+=$length_;
																$seq_line_leftover="$1\t+1\t$3\t$start_\t0\t$6\t0\t";
															}
														}
														if ($seq_line_leftover=~/\w/)
														{
															print OUT_BED_DNA qq!$seq_line_leftover\n!;
														}
													}
													else
													{
														if ($seq_line_leftover!~/\w/)
														{
															if ($seq_line_[0]=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/)
															{
																$seq_line_leftover="$1\t-1\t$3\t$5\t0\t$6\t0\t";
															}
														}
														if ($seq_line_leftover=~/\w/)
														{
															print OUT_BED_DNA qq!$seq_line_leftover\n!;
														}
														my $k=$seq_line_count_-1;
														for(;$k>=0;$k--)
														{
															my $k_=$seq_line_count_-$k-1;
															$seq_line_[$k]=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)$/$1\t$k_\t$3\t$4\t$5\t$6\t$7\t$8/;
															print OUT_BED_DNA qq!$seq_line_[$k]\n!;
														}
														if ($seq_line_after=~/\w/) { print OUT_BED_DNA qq!$seq_line_after\n!; }
													}
												} else { print qq!BED error: $bed\n!; }
											} else { print qq!MAP error: $map\n!; }
										}
									}
									if ($extra!~/\w/)
									{
										print OUT_BED_DNA qq!$bed\n!;
										print OUT_BED_DNA qq!>$name__ $description$map!;
										if ($var_protein=~/\w/) { print OUT_BED_DNA " (VAR:$var_protein)"; }
										print OUT_BED_DNA qq!\n!;
										if ($seq_line_before=~/\w/) { print OUT_BED_DNA qq!$seq_line_before\n!; }
										for(my $k=0;$k<$seq_line_count;$k++)
										{
											print OUT_BED_DNA qq!$seq_line[$k]\n!;
										}
										if ($seq_line_after=~/\w/) { print OUT_BED_DNA qq!$seq_line_after\n!; }
									}
									if (length($protein)>6)
									{
										print OUT qq!>$name__ $description$map!;
										if ($var_protein=~/\w/) { print OUT " (VAR:$var_protein)"; }
										print OUT qq!\n$protein\n!;
									}
									print OUT_BED qq!$bed\n!;
								} else { print qq!translate.pl:$line_number:$filename:Sequence Error: $name\n!; }
							} 
							else { print qq!translate.pl:$line_number:$filename:Error parsing: "$line"\n!; }
						}
					}
				}
			}
		}
		if ($block!=0) { print qq!translate.pl:$line_number:$filename: Warning: end of block not found: $name\n!;}
		close(IN);
	}	
	close(OUT);
	close(OUT_BED);
	close(OUT_BED_DNA);
}
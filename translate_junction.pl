#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $filename="";

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }

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
	if (open (IN,"$filename"))
	{
		open (OUT,">$filename.fasta");
		open (OUT_BED,">$filename.fasta.bed");
		open (OUT_BED_DNA,">$filename.fasta.bed.dna");
		my $proteins_count=0;
		my $name="";
		my $name_old="";
		my $qual="";
		my $bed="";
		my $map="";
		my $description="";
		my $seq="";
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
		my $from="";
		my $to="";
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
				$qual=$5;
				$bed=$line;
				$description="";
				$seq="";
				$seq_before="";
				$seq_after="";
				$seq_line_leftover="";
				$seq_line_before="";
				$seq_line_after="";
				if ($block!=0) { print qq!translate_junction.pl:$line_number:$filename: Warning: end of block not found: $name_old\n!;}
				$block=1;
			} 
			else 
			{
				if ($line=~/^>/)
				{
				} 
				else 
				{
					if ($line=~/^([^\t]*)\t(\-1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
					{
						$seq_before="\U$8";
					} 
					else 
					{
						if ($line=~/^([^\t]*)\t(0)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
						{
							$seq_before.="\U$8";
							$from=$4+$5;
						} 
						else 
						{ 
							if ($line=~/^([^\t]*)\t(1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
							{
								$seq_after="\U$8";
								$to=$4;
							} 
							else 
							{ 
								if ($line=~/^([^\t]*)\t(\+1)\t(chr[^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
								{
									$block=0;
									$seq_after.="\U$8";
									my $before_length=0;
									foreach my $strand ('+','-')
									{
										if ($strand=~/\+/) { $before_length=int(length($seq_before)/3); } else { $before_length=int(length($seq_after)/3); }
										for(my $k=0;$k<3;$k++)
										{
											my $protein="";
											$seq="$seq_before$seq_after";
											if ($strand=~/\-/)
											{
												$seq = reverse $seq;
												$seq=~tr/ATCG/TAGC/;
											}
											for(my $n=$k;$n<length($seq);$n=$n+3)
											{
												my $triplet = substr($seq, $n, 3);
												if (length($triplet)==3)
												{
													if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
													$protein.=$mapping{$triplet}; 
												}
											}
											my $before = substr $protein,0,$before_length;
											my $after = substr $protein,$before_length;
											$before=~s/^.*\*([^\*]+)$/$1/;
											$after=~s/^([^\*]+)\*.*$/$1/;
											my $cleavage=0;
											my $before_="";
											my $before__=$before;
											while($before__=~s/(.)(.)$//)
											{
												my $aa1=$1;
												my $aa2=$2;
												if ($cleavage<2)
												{
													$before_="$aa2$before_";
													$before__.="$aa1";
													if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
													{
														$cleavage++;
													}
												}
											}
											if ($cleavage<2) { $before_="$before__$before_"; }
											if ($cleavage==0) { $before_=""; }
											$cleavage=0;
											my $after_="";
											my $after__=$after;
											while($after__=~s/^(.)(.)//)
											{
												my $aa1=$1;
												my $aa2=$2;
												if ($cleavage<2)
												{
													$after_.="$aa1";
													$after__="$aa2$after__";
													if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
													{
														$cleavage++;
													}
												}
											}
											if ($cleavage<2) { $after_.=$after__; }
											if ($cleavage==0) { $after_=""; }
											if (length($before_)>0 and length($after_)>0 and length($before_)+length($after_)>6)
											{
												my $name_=$name;
												$name_=~s/^([^\-]+)\-.*$/$1/;
												if ($description!~/\w/) { $description="$descriptions{$name_}"; }
												if ($description=~/\w/ and $description!~/\s$/) { $description.=" "; }
												my $name__="$name-$qual-$k$strand-bridge";
												print OUT qq!>$name__ $description\n$before_$after_\n!;
												#print OUT qq!-$before $after\n--$before_ $after_\n!;
											}
										}
									} 
								} 
								else { print qq!translate_junction.pl:$line_number:$filename:Error parsing: "$line"\n!; }
							}
						}
					}
				}
			}
		}
		if ($block!=0) { print qq!translate_junction.pl:$line_number:$filename: Warning: end of block not found: $name\n!;}
		close(IN);
	}	
	close(OUT);
	close(OUT_BED);
	close(OUT_BED_DNA);
}
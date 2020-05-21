#!/usr/bin/perl
#
use strict;
use File::Copy;
#use warnings; 

my $error=0;
my $somatic="";
my $germline=""; 
my $junction="";
my $fusion="";
my $result_dir="";
my $ref_version=""; 
#my $script_root="quilts";
my $script_root="/ifs/data/proteomics/tcga/scripts/quilts/v2.0";
#my $db_root="quilts";
my $db_root="/ifs/data/proteomics/tcga/databases";
my $genome="";
my $db="";
my $db_original="";
my $var=0;
my $g_final="";
my $s_final=""; 
my $threshold_A=2;
my $threshold_AN=3;
my $threshold_N=3;
my $star=0;
if ($ARGV[0]=~/\w/) { $ref_version=$ARGV[0];} else { $ref_version="ensembl_human_37.70"; }
if ($ARGV[1]=~/\w/) { $somatic=$ARGV[1]; } 
if ($ARGV[2]=~/\w/) { $germline=$ARGV[2]; } 
if ($ARGV[3]=~/\w/) { $junction=$ARGV[3]; } 
if ($ARGV[4]=~/\w/) { $fusion=$ARGV[4]; }
if ($ARGV[5]=~/\w/) { $result_dir=$ARGV[5]; } else { $error=1; }
if ($ARGV[6]=~/\w/) { $genome=$ARGV[6]; } else { $genome="$db_root/genome_human"; }
if ($ARGV[7]=~/\w/) { $threshold_A=$ARGV[7];} else { $threshold_A=2 } 
if ($ARGV[8]=~/\w/) { $threshold_AN=$ARGV[8];} else { $threshold_AN=3; } 
if ($ARGV[9]=~/\w/) { $threshold_N=$ARGV[9];} else { $threshold_N=3; } 
if ($ARGV[10]=~/\w/) { $star=$ARGV[10];} else { $star=0; } 
$db="$db_root/$ref_version";
$db_original=$db;
$db_original=~s/\_([^\_]+)$/_original_$1/;
open (LOG, ">$result_dir/cgi_log.txt");
if ($somatic!~/\w/ and $germline!~/\w/ and $junction!~/\w/ and $fusion!~/\w/)
{
	if (open (STATUS, ">$result_dir/status.txt")) { print STATUS "No input\n"; close(STATUS); }
	$error=1;
}
if ($error==0)
{
	system(qq!date!);
	my $step=0;
	my $steps_total=0;
	if ($somatic=~/\w/ or $germline=~/\w/) { $steps_total++; }
	if ($junction=~/\w/) { $steps_total+=3; }
	if (open (STATUS, ">>$result_dir/status.txt")) { print STATUS "Started\n"; close(STATUS); }
	print LOG "QUILTS Version 1.0\nReference DB used: $ref_version\n"; 
	system(qq!mkdir "$result_dir"!);
	system(qq!rm -f "$result_dir/*.txt"!);
	system(qq!rm -Rf "$result_dir/log"!);
	system(qq!rm -Rf "$result_dir/fasta"!); 
	system(qq!mkdir "$result_dir/log"!);
	system(qq!mkdir "$result_dir/fasta"!); 
	system(qq!mkdir "$result_dir/fasta/parts"!); 
	system(qq!cp "$db_original/proteome.fasta" "$result_dir/fasta/parts/$ref_version.fasta"!);
	if (open(TEST,"$result_dir/log/frameshift.bed.dna")) { close(TEST);	system(qq!rm "$result_dir/log/frameshift.bed.dna"!); }
	if ($somatic=~/\w/)
	{
		system(qq!perl $script_root/merge_vcf.pl $somatic!);
	}
	if ($germline=~/\w/)
	{
		system(qq!perl $script_root/merge_vcf.pl $germline!);
	}
	if ($somatic=~/\w/ or $germline=~/\w/)
	{
		system(qq!cp "$db/proteome.bed" "$result_dir/log/proteome.bed"!);
		system(qq!$script_root/read_chr_bed "$result_dir/log/proteome.bed"!);
		if ($somatic=~/\w/)
		{
			system(qq!perl $script_root/get_variants.pl "$result_dir/log/proteome.bed" $somatic/merged/merged.vcf S!);
			system(qq!mv $result_dir/log/proteome.bed.var $result_dir/log/proteome.bed.S.var!);
		}
		if ($germline=~/\w/)
		{
			system(qq!perl $script_root/get_variants.pl "$result_dir/log/proteome.bed" $germline/merged/merged.vcf G!);
			system(qq!mv $result_dir/log/proteome.bed.var $result_dir/log/proteome.bed.G.var!);
		}	
		if ($somatic=~/\w/) { system(qq!cat $result_dir/log/proteome.bed.S.var > "$result_dir/log/proteome.bed.var"!); }
		if ($germline=~/\w/)
		{
			if ($somatic=~/\w/) { system(qq!cat $result_dir/log/proteome.bed.G.var >> "$result_dir/log/proteome.bed.var"!); }
			else { system(qq!cat $result_dir/log/proteome.bed.G.var > "$result_dir/log/proteome.bed.var"!); }
		}
		system(qq!perl $script_root/sort_variants.pl "$result_dir/log/proteome.bed.dna" "$result_dir/log/proteome.bed.var"!);
		system(qq!cp "$db/proteome-descriptions.txt" "$result_dir/log/proteome.aa.var-descriptions.txt"!);
		system(qq!cp "$db/proteome-genes.txt" "$result_dir/log/proteome.aa.var-genes.txt"!);
		system(qq!perl $script_root/translate.pl "$result_dir/log/proteome.aa.var.bed.dna"!);
		system(qq!cp "$result_dir/log/proteome.aa.var.bed.dna.fasta" "$result_dir/fasta/parts/variants.fasta"!);
		system(qq!cp "$result_dir/log/proteome.aa.var.bed.dna.fasta.bed" "$result_dir/fasta/parts/variants.bed"!);
		system(qq!cat "$result_dir/fasta/parts/variants.fasta" >> "$result_dir/fasta/parts/proteome.fasta"!);
		system(qq!cat "$result_dir/fasta/parts/variants.bed" >> "$result_dir/fasta/parts/proteome.bed"!);
		system(qq!perl $script_root/translate.pl "$result_dir/log/proteome.stop.var.bed.dna"!);
		system(qq!cat "$result_dir/log/proteome.stop.var.bed.dna.fasta.bed.dna" > "$result_dir/log/frameshift2.bed.dna"!);
		system(qq!cat "$result_dir/log/proteome.stop.var.bed.dna.fasta.bed" > "$result_dir/log/frameshift2.bed"!);
		
		system(qq!perl $script_root/translate.pl "$result_dir/log/proteome.stop-removed.var.bed.dna" 1!);
		system(qq!cat "$result_dir/log/proteome.stop-removed.var.bed.dna.fasta.bed.dna" >> "$result_dir/log/frameshift2.bed.dna"!);
		system(qq!cat "$result_dir/log/proteome.stop-removed.var.bed.dna.fasta.bed" >> "$result_dir/log/frameshift2.bed"!);
		
		if ($somatic=~/\w/)
		{
			system(qq!perl $script_root/get_variants.pl "$result_dir/log/frameshift2.bed" $somatic/merged/merged.vcf S!);
			system(qq!mv $result_dir/log/frameshift2.bed.var $result_dir/log/frameshift2.bed.S.var!);
		}
		if ($germline=~/\w/)
		{
			system(qq!perl $script_root/get_variants.pl "$result_dir/log/frameshift2.bed" $germline/merged/merged.vcf G!);
			system(qq!mv $result_dir/log/frameshift2.bed.var $result_dir/log/frameshift2.bed.G.var!);
		}
		if ($somatic=~/\w/) { system(qq!cat $result_dir/log/frameshift2.bed.S.var > "$result_dir/log/frameshift2.bed.var"!); }
		if ($germline=~/\w/)
		{
			if ($somatic=~/\w/) { system(qq!cat $result_dir/log/frameshift2.bed.G.var >> "$result_dir/log/frameshift2.bed.var"!); }
			else { system(qq!cat $result_dir/log/frameshift2.bed.G.var > "$result_dir/log/frameshift2.bed.var"!); }
		}
		system(qq!perl $script_root/sort_variants.pl "$result_dir/log/frameshift2.bed.dna" "$result_dir/log/frameshift2.bed.var"!);
		system(qq!cat $result_dir/log/frameshift2.bed.dna >> "$result_dir/log/frameshift.bed.dna"!);
		system(qq!cat $result_dir/log/frameshift2.aa.var.bed.dna >> "$result_dir/log/frameshift.bed.dna"!); 
	}
	
	if ($junction=~/\w/)
	{
		if($star==1)
		{
			print qq!STAR:$junction\n!;
			print qq!perl $script_root/merge_junctions_bed_dir_star.pl $junction\n!;
			system(qq!perl $script_root/merge_junctions_bed_dir_star.pl $junction!);
			$junction.="/merged";
			print qq!STAR:modified:$junction\n!;
		}
		system(qq!perl $script_root/merge_junctions_bed_dir.pl "$junction" "$result_dir/log/"!);
		system(qq!perl $script_root/filter_known_transcripts.pl "$result_dir/log/merged-junctions.bed" "$db/transcriptome.bed"!);
		system(qq!perl $script_root/filter_A.pl "$result_dir/log/merged-junctions.filter.bed" "$db/proteome.bed" 0 $threshold_A $threshold_AN $threshold_N!);
		system(qq!$script_root/read_chr_bed "$result_dir/log/merged-junctions.filter.bed.A.bed"!);
		system(qq!cp "$db/proteome-descriptions.txt" "$result_dir/log/alternative-descriptions.txt"!);
		system(qq!cp "$db/proteome-genes.txt" "$result_dir/log/alternative-genes.txt"!);
		system(qq!cp $result_dir/log/merged-junctions.filter.bed.A.bed.dna "$result_dir/log/alternative.bed.dna"!);
		if ($somatic=~/\w/ or $germline=~/\w/)
		{
			if ($somatic=~/\w/)
			{
				system(qq!perl $script_root/get_variants.pl "$result_dir/log/merged-junctions.filter.bed.A.bed" $somatic/merged/merged.vcf S!);
				system(qq!mv $result_dir/log/merged-junctions.filter.bed.A.bed.var $result_dir/log/merged-junctions.filter.bed.A.bed.S.var!);
			}
			if ($germline=~/\w/)
			{
				system(qq!perl $script_root/get_variants.pl "$result_dir/log/merged-junctions.filter.bed.A.bed" $germline/merged/merged.vcf G!);
				system(qq!mv $result_dir/log/merged-junctions.filter.bed.A.bed.var $result_dir/log/merged-junctions.filter.bed.A.bed.G.var!);
			}
			if ($somatic=~/\w/) { system(qq!cat $result_dir/log/merged-junctions.filter.bed.A.bed.S.var > "$result_dir/log/merged-junctions.filter.bed.A.bed.var"!); }
			if ($germline=~/\w/)
			{
				if ($somatic=~/\w/) { system(qq!cat $result_dir/log/merged-junctions.filter.bed.A.bed.G.var >> "$result_dir/log/merged-junctions.filter.bed.A.bed.var"!); }
				else { system(qq!cat $result_dir/log/merged-junctions.filter.bed.A.bed.G.var > "$result_dir/log/merged-junctions.filter.bed.A.bed.var"!); }
			}
			system(qq!perl $script_root/sort_variants.pl "$result_dir/log/merged-junctions.filter.bed.A.bed.dna" "$result_dir/log/merged-junctions.filter.bed.A.bed.var"!);
			system(qq!cat $result_dir/log/merged-junctions.filter.bed.A.aa.var.bed.dna >> "$result_dir/log/alternative.bed.dna"!);
		}
		system(qq!perl $script_root/translate.pl "$result_dir/log/alternative.bed.dna"!);
		system(qq!cp "$result_dir/log/alternative.bed.dna.fasta" "$result_dir/fasta/parts/alternative.fasta"!);
		system(qq!cp "$result_dir/log/alternative.bed.dna.fasta.bed" "$result_dir/fasta/parts/alternative.bed"!);
		system(qq!cat "$result_dir/fasta/parts/alternative.fasta" >> "$result_dir/fasta/parts/proteome.fasta"!);
		system(qq!cat "$result_dir/fasta/parts/alternative.bed" >> "$result_dir/fasta/parts/proteome.bed"!);

		system(qq!cp "$result_dir/log/merged-junctions.filter.bed.A_.bed" "$result_dir/log/frameshift1.bed"!);
		system(qq!cat "$result_dir/log/merged-junctions.filter.bed.AN.bed" >> "$result_dir/log/frameshift1.bed"!);
		system(qq!$script_root/read_chr_bed "$result_dir/log/frameshift1.bed"!);
		system(qq!perl $script_root/translate.pl "$result_dir/log/frameshift1.bed.dna" 1!);
		system(qq!cp "$result_dir/log/frameshift1.bed.dna.fasta.bed" "$result_dir/log/frameshift3.bed"!);
		system(qq!cp "$result_dir/log/frameshift1.bed.dna.fasta.bed.dna" "$result_dir/log/frameshift3.bed.dna"!);
					
		if ($somatic=~/\w/ or $germline=~/\w/)
		{
			if ($somatic=~/\w/)
			{
				system(qq!perl $script_root/get_variants.pl "$result_dir/log/frameshift3.bed" $somatic/merged/merged.vcf S!);
				system(qq!mv $result_dir/log/frameshift3.bed.var $result_dir/log/frameshift3.bed.S.var!);
			}
			if ($germline=~/\w/)
			{
				system(qq!perl $script_root/get_variants.pl "$result_dir/log/frameshift3.bed" $germline/merged/merged.vcf G!);
				system(qq!mv $result_dir/log/frameshift3.bed.var $result_dir/log/frameshift3.bed.G.var!);
			}
			if ($somatic=~/\w/) { system(qq!cat $result_dir/log/frameshift3.bed.S.var > "$result_dir/log/frameshift3.bed.var"!); }
			if ($germline=~/\w/)
			{
				if ($somatic=~/\w/) { system(qq!cat $result_dir/log/frameshift3.bed.G.var >> "$result_dir/log/frameshift3.bed.var"!); }
				else { system(qq!cat $result_dir/log/frameshift3.bed.G.var > "$result_dir/log/frameshift3.bed.var"!); }
			}
			system(qq!perl $script_root/sort_variants.pl "$result_dir/log/frameshift3.bed.dna" "$result_dir/log/frameshift3.bed.var"!);
		}
		system(qq!cat $result_dir/log/frameshift3.bed.dna >> "$result_dir/log/frameshift.bed.dna"!);
		if ($somatic=~/\w/ or $germline=~/\w/) 
		{ 
			system(qq!cat $result_dir/log/frameshift3.aa.var.bed.dna >> "$result_dir/log/frameshift.bed.dna"!); 
		}
		
		system(qq!date!);
		system(qq!$script_root/read_chr_bed $result_dir/log/merged-junctions.filter.bed.notA.bed!);
		system(qq!perl $script_root/translate_junction.pl $result_dir/log/merged-junctions.filter.bed.notA.bed.dna!);
		system(qq!cp "$result_dir/log/merged-junctions.filter.bed.notA.bed.dna.fasta" "$result_dir/fasta/parts/other.fasta"!);
		system(qq!cat "$result_dir/fasta/parts/other.fasta" >> "$result_dir/fasta/parts/proteome.fasta"!);
		system(qq!date!);
	}
	else
	{
		if ($somatic=~/\w/ or $germline=~/\w/)
		{

		}
	}
	if ($somatic=~/\w/ or $germline=~/\w/ or $junction=~/\w/)
	{
		system(qq!cp "$db/proteome-descriptions.txt" "$result_dir/log/frameshift-descriptions.txt"!);
		system(qq!cp "$db/proteome-genes.txt" "$result_dir/log/frameshift-genes.txt"!);
		system(qq!perl $script_root/translate.pl "$result_dir/log/frameshift.bed.dna" 1!);
		system(qq!cat "$result_dir/log/frameshift.bed.dna.fasta" >> "$result_dir/fasta/parts/frameshift.fasta"!);
		system(qq!cat "$result_dir/log/frameshift.bed.dna.fasta.bed" >> "$result_dir/fasta/parts/frameshift.bed"!);
		system(qq!cat "$result_dir/fasta/parts/frameshift.fasta" >> "$result_dir/fasta/parts/proteome.fasta"!);
		system(qq!cat "$result_dir/fasta/parts/frameshift.bed" >> "$result_dir/fasta/parts/proteome.bed"!);
	}
	#if ($fusion!~/\w/)
	#{
	#	system (qq!perl $script_root/fusion.pl "$fusion" "$result_dir/log"!);
	#	system (qq!cat "$result_dir/log/fusion.fasta">> "$result_dir/fasta/parts/proteome.fasta"!);
	#}
	system(qq!perl $script_root/check_and_clean_fasta.pl "$result_dir/fasta/parts" Y!);
	system(qq!perl $script_root/check_and_clean_fasta.pl "$result_dir/fasta/parts" N!);
	system(qq!cp "$result_dir/fasta/parts/proteome.fasta.cleaned.fasta" "$result_dir/fasta/proteome.fasta"!);
	system(qq!cp "$result_dir/fasta/parts/$ref_version.fasta.cleaned.fasta" "$result_dir/fasta/$ref_version.fasta"!);
	system(qq!perl $script_root/check_and_clean_fasta.pl "$result_dir/fasta" N!);
	system(qq!date!);
	if (open (STATUS, ">>$result_dir/status.txt")) { print STATUS "Done\n"; close(STATUS); }
}
close LOG; 

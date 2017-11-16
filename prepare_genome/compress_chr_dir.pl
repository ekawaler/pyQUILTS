#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $dir="";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }

if ($error==0)
{
	opendir(DIR,"$dir");
	my @files = readdir DIR;
	closedir DIR;
	foreach my $filename (@files)
	{
		if ($filename=~/\.fa$/i)
		{
			system(qq!/ifs/data/proteomics/tcga/scripts/quilts/v2.0/compress_chr $dir/$filename!)
		}
	}
}

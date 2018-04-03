#!c:/perl/bin/perl.exe
#
use strict;
use Cwd;

my $error=0;
my $dir="";
my $curdir = getcwd;

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir=$curdir; }

if ($error==0)
{
	opendir(DIR,"$dir");
	my @files = readdir DIR;
	closedir DIR;
	foreach my $filename (@files)
	{
		if ($filename=~/\.fa$/i)
		{
			system(qq!'$curdir/compress_chr' '$dir/$filename'!)
		}
	}
}

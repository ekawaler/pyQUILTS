#!c:/perl/bin/perl.exe
#
# Unzips every .gz file in a directory.

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
		if ($filename=~/\.gz$/i)
		{
			system(qq!gunzip $dir/$filename!)
		}
	}
}

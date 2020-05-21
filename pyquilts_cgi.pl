#!/usr/bin/perl

#local $ENV{PATH}="$ENV{PATH}:/data/python/2.7.15/bin";
#local $ENV{CPATH}="$ENV{CPATH}:/data/python/2.7.15/include/python2.7";
#local $ENV{MANPATH}="$ENV{MANPATH}:/data/python/2.7.15/share/man";
#local $ENV{PYTHONPATH}="$ENV{PYTHONPATH}:/data/python/2.7.15/lib/python2.7/site-packages";

use CGI;
print "Content-type: text/html\n\n";
my $query = new CGI;
my $status=$query->param("status");
my $quilts_id = $query->param("quilts_id");
my $result_dir = $quilts_id;
my $junc_type = $query->param("junc_type");
$result_dir=~s/^quilts\-//;
$result_dir=~s/\-[^\-]+$//;
$result_dir=~s/\-/\./m;
$result_dir="results_$result_dir";
my $reference= $query->param("reference");
my $ok=0;
if ($status!~/\w/)
{
	print qq!
		<table cellpadding="0" cellspacing="0" border="0"><tr>
		<td width="120" valign="top"><img alt="" src="/quilts/quilts.jpg" width="90" height="90"></td>
		<td width="5"></td>
		<td width="500" valign="middle">
		<b>Welcome to QUILTS</b>, a tool for creating sample specific protein sequence databases. It uses genomic and transcriptomic information to create comprehensive sample specific protein database that supports the identification of novel proteins, resulting from single nucleotide variants, splice variants and fusion genes.   
		</td>
		</tr></table><p>
		<FORM ACTION="pyquilts_cgi.pl" METHOD="post" ENCTYPE="multipart/form-data">
		<INPUT TYPE="hidden" NAME="status" VALUE="QUILTS">
		To create your sample specific protein database, please <b>choose a reference database</b>: <br>
		<INPUT TYPE="radio" NAME="reference" VALUE="ensembl_human_37.70" checked>Ensembl (GRCh37/hg19)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
		<INPUT TYPE="radio" NAME="reference" VALUE="refseq_hg19_20170922">RefSeq (GRCh37/hg19)<br><br>
		If you're using junctions, please <b>choose your junction software</b> (if not, ignore this): <br>
		<INPUT TYPE="radio" NAME="junc_type" VALUE="mapsplice" checked>MapSplice&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
		<INPUT TYPE="radio" NAME="junc_type" VALUE="tophat">TopHat&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
		<INPUT TYPE="radio" NAME="junc_type" VALUE="star">STAR<br><br>
		And upload at least one sample-specific file: <br>
		<b>Germline variants (VCF)</b> (<a href="/quilts/1000G_chr22.vcf">Example</a>, only first 6 columns necessary) : <INPUT TYPE="file" NAME="data_file1"><br>
		<b>Somatic variants (VCF)</b> (<a href="/quilts/somatic_chr22.vcf">Example</a>, only first 6 columns necessary) : <INPUT TYPE="file" NAME="data_file2"><br>
		<b>Junction (TXT for MapSplice, BED for TopHat, TAB for STAR)</b> (<a href="/quilts/junctions_chr22.bed">Example</a>) : <INPUT TYPE="file" NAME="data_file3"><br>
		<b>Fusion (TXT)</b> (<a href="/quilts/fusions.txt">Example</a>, fifth column can be any quality measure) : <INPUT TYPE="file" NAME="data_file4"><br><br>
		<b>Variant quality threshold:</b> <INPUT TYPE="text" NAME="variant_quality_threshold" VALUE="0.0"><br><br>
		<b>Threshold of supporting reads for novel splice junctions</b><br>
		Both boundaries annotated: <INPUT TYPE="text" NAME="threshA" VALUE="2"><br>
		Left boundary annotated: <INPUT TYPE="text" NAME="threshAN" VALUE="3"><br>
		No boundaries annotated: <INPUT TYPE="text" NAME="threshN" VALUE="3"><br><br>
		<INPUT TYPE="checkbox" NAME="no_missed_cleavage" VALUE="no_missed_cleavage"> <b>No missed cleavage</b><br><br>
	!;
						
	print qq!
		Once all data is uploaded click here: 
		<INPUT TYPE="submit" NAME="Submit" VALUE="QUILTS">
		</FORM><p>
		<br><br>
	!;
}
else
{		
	if ($status=~/^QUILTS$/)
	{
		my $reference = $query->param("reference");
		my $junc_type = $query->param("junc_type");
		my $variant_quality_threshold = $query->param("variant_quality_threshold");
		my $threshA = $query->param("threshA");
		my $threshAN = $query->param("threshAN");
		my $threshN = $query->param("threshN");
		my $no_missed_cleavage = $query->param("no_missed_cleavage");
		$id = GetDateTime();
		$id .= "-" . int(123456789*rand());
		$quilts_id="quilts-$id";
		mkdir("../html/temp/$quilts_id");
		mkdir("../html/temp/$quilts_id/germline");
		mkdir("../html/temp/$quilts_id/somatic");
		mkdir("../html/temp/$quilts_id/junctions");
		mkdir("../html/temp/$quilts_id/fusion");
		
		my $upload_germline = $query->upload("data_file1");
		my $germline="";
		if (open UPLOADFILE, ">../html/temp/$quilts_id/germline/germline.vcf")
		{
			while ( $line=<$upload_germline> ) 
			{
				chomp($line);
				$line=~s/\r$//;
				$line=~s/\r([^\n])/\n$1/g;
				#if ($line=~/^chr[0-9XYM]+\t[0-9]+/)
				#{
					print UPLOADFILE "$line\n";
					$ok++;
				#}
			}
			close UPLOADFILE;
		}
		if ($upload_germline=~/\.vcf$/){$germline="/var/www/html/temp/$quilts_id/germline";}
		else { $germline=""; }
		
		my $upload_somatic = $query->upload("data_file2");
		my $somatic="";
		if (open UPLOADFILE, ">../html/temp/$quilts_id/somatic/somatic.vcf")
		{
			while ( $line=<$upload_somatic> ) 
			{
				chomp($line);
				$line=~s/\r$//;
				$line=~s/\r([^\n])/\n$1/g;
				#if ($line=~/^chr[0-9XYM]+\t[0-9]+/)
				#{
					print UPLOADFILE "$line\n";
					$ok++;
				#}
			}
			close UPLOADFILE;
		}
		if ($upload_somatic=~/\.vcf$/){$somatic="/var/www/html/temp/$quilts_id/somatic";}
		else { $somatic=""; }
		
		my $upload_junction = $query->upload("data_file3");
		my $junction="";
		if ($junc_type=='tophat'){
			if (open UPLOADFILE, ">../html/temp/$quilts_id/junctions/junctions.bed")
			{
				while ( $line=<$upload_junction> ) 
				{
					chomp($line);
					$line=~s/\r$//;
					$line=~s/\r([^\n])/\n$1/g;
					#if ($line=~/^chr[0-9XYM]+\t[0-9]+/)
					#{
						print UPLOADFILE "$line\n";
						$ok++;
					#}
				}
				close UPLOADFILE;
			}
		}
		elsif ($junc_type=='star'){
			if (open UPLOADFILE, ">../html/temp/$quilts_id/junctions/junctions.tab")
			{
				while ( $line=<$upload_junction> ) 
				{
					chomp($line);
					$line=~s/\r$//;
					$line=~s/\r([^\n])/\n$1/g;
					#if ($line=~/^chr[0-9XYM]+\t[0-9]+/)
					#{
						print UPLOADFILE "$line\n";
						$ok++;
					#}
				}
				close UPLOADFILE;
			}
		}
		elsif ($junc_type=='mapsplice'){
			if (open UPLOADFILE, ">../html/temp/$quilts_id/junctions/junctions.txt")
			{
				while ( $line=<$upload_junction> ) 
				{
					chomp($line);
					$line=~s/\r$//;
					$line=~s/\r([^\n])/\n$1/g;
					#if ($line=~/^chr[0-9XYM]+\t[0-9]+/)
					#{
						print UPLOADFILE "$line\n";
						$ok++;
					#}
				}
				close UPLOADFILE;
			}
		}
		if ($upload_junction=~/\.bed$/ || $upload_junction=~/\.txt$/ || $upload_junction=~/\.tab$/){$junction="/var/www//html/temp/$quilts_id/junctions";}
		else { $junction="";}
		
		$fusion="";
		my $upload_fusion = $query->upload("data_file4");
		if (open UPLOADFILE, ">../html/temp/$quilts_id/fusion/fusion.txt")
		{
			while ( $line=<$upload_fusion> ) 
			{
				chomp($line);
				$line=~s/\r$//;
				$line=~s/\r([^\n])/\n$1/g;
				print UPLOADFILE "$line\n";
				$ok++;
			}
			close UPLOADFILE;
		}
		if ($upload_fusion=~/\.txt$/){$fusion="../html/temp/$quilts_id/fusion";}
		else { $fusion=""; }
		
		print qq!The processing has started and the results will be available <a href="pyquilts_cgi.pl?status=QuiltsResult\&quilts_id=$quilts_id\&reference=$reference">here</a>.<p>\n!;
		if ($ok>0)
		{
			my $count_input_files=0;
			$command_line=qq!python pyquilts/quilts.py --output_dir ../html/temp/$quilts_id --proteome pyquilts/$reference --genome pyquilts/hg19!;
			if ($germline=~/\w/)
			{
				$command_line.=qq! --germline $germline!;
				$count_input_files++;
			}
			if ($somatic=~/\w/)
			{
				$command_line.=qq! --somatic $somatic!;
				$count_input_files++;
			}
			if ($junction=~/\w/)
			{
				$command_line.=qq! --junction $junction --junction_file_type $junc_type!;
				$count_input_files++;
			}
			if ($fusion=~/\w/)
			{
				$command_line.=qq! --fusion $fusion!;
				$count_input_files++;
			}
									
			if ($variant_quality_threshold=~/\w/)
			{
				$command_line.=qq! --variant_quality_threshold $variant_quality_threshold!;
			}
			if ($threshA=~/\w/)
			{
				$command_line.=qq! --threshA $threshA!;
			}
			if ($threshAN=~/\w/)
			{
				$command_line.=qq! --threshAN $threshAN!;
			}
			if ($threshN=~/\w/)
			{
				$command_line.=qq! --threshN $threshN!;
			}
			if ($no_missed_cleavage=~/\w/)
			{
				$command_line.=qq! --no_missed_cleavage!;
			}
			$command_line.=qq! > ../html/temp/$quilts_id/log.txt 2>\&1!;
			if ($count_input_files>0)
			{
				print qq!$command_line\n!;
				if (open UPLOADFILE, ">../html/temp/$quilts_id/pyquilts.sh")
				{
					print UPLOADFILE qq!
export PATH=/data/python/2.7.15/bin:\$PATH 
export CPATH=/data/python/2.7.15/include/python2.7:\$CPATH 
export MANPATH=/ldata/python/2.7.15/share/man:\$MANPATH 
export PYTHONPATH=/data/python/2.7.15/lib/python2.7/site-packages:\$PYTHONPATH 
cd /var/www/cgi-bin/
$command_line
										!;
					close UPLOADFILE;
				}
				system("chmod u=rwx,g=rwx,o=rwx ../html/temp/$quilts_id/pyquilts.sh");
				system("chgrp slice ../html/temp/$quilts_id/pyquilts.sh"); 
				system("chcon -t httpd_sys_script_exec_t ../html/temp/$quilts_id/pyquilts.sh"); 
				#system("../html/temp/$quilts_id/pyquilts.sh >../html/temp/$quilts_id/test.log 2>&1");
				system("$command_line");
				system("chmod -R u=rwx,g=rwx,o=rwx ../html/temp/$quilts_id");
			}
			else
			{
				print qq!No input files were submitted\n!;
			}
		}
	}		
	if ($status=~/^QuiltsResult$/)
	{
		if (open(IN,"../html/temp/$quilts_id/$result_dir/status.txt"))
		{
			my $message=qq!Processing in progress...<p>\n!;
			while($line=<IN>)
			{
				if ($line=~/^Started/) { $message=qq!Processing in progress...<br>\n!; }
				if ($line=~/^No input/) { $message=qq!Please submit at least one sample-specific file<br>\n!; }
				if ($line=~/^([0-9\.]+)\%/) { $message.=qq!Completed $1\%<br>\n!; }
				if ($line=~/\-([0-9\.]+)\%/) { $message.=qq!-\n!; }
				if ($line=~/^DONE/) 
				{
					$status="QuiltsResultDone"; 
					$message.=qq!Processing done<br>\n!; 
					#$message="";
				}
			}
			close(IN);
			if ($status!~/QuiltsResultDone/) { print qq!<head>\n<META HTTP-EQUIV="refresh" CONTENT="10">\n</head>\n!; }
			print qq!<body>\n!;
			print qq!$message\n!;
			print qq!</body>\n!;
		}
	}
	if ($status=~/^QuiltsResultDone$/)
	{
		print qq!<br><b>Results:</b><br>!;
		print qq!<a href="/temp/$quilts_id/$result_dir/fasta/reference_proteome.fasta">reference_proteome.fasta</a><br>!;
		print qq!<a href="/temp/$quilts_id/$result_dir/fasta/variant_proteome.fasta">Sample-specific FASTA file</a><br>!;
		print qq!<a href="/quilts/quilts_label_key_v2.txt">Key to the FASTA header</a><br>!;
	}
}

sub GetDateTime
{
	my $sec="";
	my $min="";
	my $hour="";
	my $mday="";
	my $mon="";
	my $year="";

	($sec,$min,$hour,$mday,$mon,$year) = localtime();

	if ($sec<10) { $sec="0$sec"; }
	if ($min<10) { $min="0$min"; }
	if ($hour<10) { $hour="0$hour"; }
	if ($mday<10) { $mday="0$mday"; }
	$mon++;
	if ($mon<10) { $mon="0$mon"; }
	$year+=1900;
	$date="$year$mon$mday-$hour$min$sec";
	
	return $date;
}

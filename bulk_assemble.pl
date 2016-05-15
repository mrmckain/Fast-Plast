#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0];

while(<$file>){
	chomp;
	my @tarray=split/\s+/;
	
	my @files = </home/mmckain/Sequence_Value/Illinois_05032016/$tarray[0]*.fastq>;

	my $file1;
	my $file2;
	for my $tfile (@files){
		if($tfile =~ /R1/){
			$file1=$tfile;
		}
		else{
			$file2=$tfile;
		}
	}
	`mkdir $tarray[1]_$tarray[0]`;
	`cd $tarray[1]_$tarray[0]`;
	`cp ../plastome-assembly_condor_submit`;
	`perl -pi -e 's/REPLACE/$file1 $file2 $tarray[1]/g' plastome-assembly_condor_submit`;
	`condor_submit plastome-assembly_condor_submit`;
	#`rm map* *sam *trimmed*`;
	`cd ../`;
}	
		

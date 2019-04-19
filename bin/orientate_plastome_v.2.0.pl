#!/usr/bin/env perl
use strict;
use warnings;
my $lsc;
my $irb;
my $ssc;

my %sequences;
my $sid;

my $lsc_id;
my $max_sc=0;
my $ssc_id;
my $lsc_end;
my $ir_end;
my $ir_relative_position;
my $ir_id;

open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	if(/>/){
		$sid = substr($_,1);
		if($sid =~ /sc_/){
		$sid =~ /sc_(\d+)-(\d+)/;
		my $length = $2-$1;
		if($length > $max_sc){
			if($lsc_id){
				$ssc_id=$lsc_id;
				$lsc_id=$sid;
				$max_sc=$length;
				$lsc_end = $2;
			}
			else{
				$lsc_id=$sid;
				$max_sc=$length;
				$lsc_end = $2;
			}
		}
		else{
			$ssc_id = $sid;
		}
	}
		if($sid =~ /ir_/){
			$ir_id=$sid;
			$sid =~ /ir_(\d+)-(\d+)/;
			$ir_end =$2;
		}
	}
	else{
		$sequences{$sid}.=$_;
	}
}

if($lsc_end > $ir_end){
	$ir_relative_position = "-";

}
else{
	$ir_relative_position="+";
}
my $ir_pos;
my $lsc_pos;
my $ssc_pos;
if($ir_relative_position eq "+"){
	$ir_pos = 2;
	$lsc_pos = 1;
	$ssc_pos =3;
}
else{
	$ir_pos = 2;
	$lsc_pos = 3;
	$ssc_pos =1;
}

my %gene_orientations;
$gene_orientations{lsc}{"-"}=0;
$gene_orientations{lsc}{"+"}=0;
$gene_orientations{ssc}{"-"}=0;
$gene_orientations{ssc}{"+"}=0;
$gene_orientations{ir}{"-"}=0;
$gene_orientations{ir}{"+"}=0;

my $ir_rrns;
open my $blast, "<", $ARGV[1];
while(<$blast>){
	chomp;
	my @tarray = split/\s+/;
	if($tarray[0] eq $lsc_id){
		if($tarray[1] =~ /rpl/ || $tarray[1] =~ /rps/){
			if($tarray[9]-$tarray[8] > 0) {
				$gene_orientations{lsc}{"-"}++;
			}
			else{
				$gene_orientations{lsc}{"+"}++;
			}

		}
	}
	elsif($tarray[0] eq $ssc_id){
		if($tarray[9]-$tarray[8] > 0) {
				$gene_orientations{ssc}{"-"}++;
			}
			else{
				$gene_orientations{ssc}{"+"}++;
			}
	}
	else{
		if($tarray[1] =~ /rrn/){
		if($tarray[9]-$tarray[8] > 0) {
				$gene_orientations{ir}{"-"}++;
			}
			else{
				$gene_orientations{ir}{"+"}++;
			}
	}
}
}
if($gene_orientations{lsc}{"+"} > $gene_orientations{lsc}{"-"}){
	$lsc = reverse($sequences{$lsc_id});
	$lsc =~ tr/ATGCatgc/TACGtacg/;
}
else{
	$lsc = $sequences{$lsc_id};
}
if($gene_orientations{ssc}{"+"} > $gene_orientations{ssc}{"-"}){
	$ssc = reverse($sequences{$ssc_id});
	$ssc =~ tr/ATGCatgc/TACGtacg/;
}
else{
	$ssc = $sequences{$ssc_id};
}
if($gene_orientations{ir}{"+"} > $gene_orientations{ir}{"-"}){
			$irb = reverse($sequences{$ir_id});
			$irb =~ tr/ATGCatgc/TACGtacg/;
}
else{
	$irb = $sequences{$ir_id};
}

my $ira = reverse($irb);
$ira =~ tr/ATCGatcg/TAGCtagc/;
my $fullcp = $lsc . $irb . $ssc . $ira;
open my $out, ">", $ARGV[2] ."_CP_pieces.fsa";
print $out ">lsc\n$lsc\n>irb\n$irb\n>ssc\n$ssc\n";

open my $fullout, ">", "$ARGV[2]\_FULLCP.fsa";
print $fullout ">$ARGV[2]\n$fullcp\n";


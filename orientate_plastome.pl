#!/usr/bin/perl -w
use strict;
my $lsc;
my $irb;
my $ssc;

my %sequences;
my $sid;
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	if(/>/){
		$sid = substr($_,1);
	}
	else{
		$sequences{$sid}.=$_;
	}
}


open my $blast, "<", $ARGV[1];
while(<$blast>){
	chomp;
	my @tarray = split/\s+/;
	if($tarray[1] eq "psbA"){
		if($tarray[7]-$tarray[6] > 0){
			if($tarray[9]-$tarray[8] > 0) {
				$lsc = reverse($sequences{$tarray[0]});
				$lsc =~ tr/ATGCatgc/TACGtacg/;
			}
			else{
				$lsc = $sequences{$tarray[0]};
			}
		}
		else{
			if($tarray[9]-$tarray[8] < 0) {
				$lsc = reverse($sequences{$tarray[0]});
				$lsc =~ tr/ATGCatgc/TACGtacg/;
			}
			else{
				$lsc = $sequences{$tarray[0]};
			}
		}
	}
	if($tarray[1] eq "rps19"){
		if($tarray[7]-$tarray[6] > 0){
			if($tarray[9]-$tarray[8] > 0 ){
				$irb = reverse($sequences{$tarray[0]});
				$irb =~ tr/ATGCatgc/TACGtacg/;
			}
			else{
				$irb = $sequences{$tarray[0]};
			}
		}
		else{
			if($tarray[9]-$tarray[8] < 0) {
				$irb = reverse($sequences{$tarray[0]});
				$irb =~ tr/ATGCatgc/TACGtacg/;
			}
			else{
				$irb = $sequences{$tarray[0]};
			}
		}
	}
	if($tarray[1] eq "ndhF"){
		if($tarray[7]-$tarray[6] > 0){
			if($tarray[9]-$tarray[8] > 0) {
				$ssc = reverse($sequences{$tarray[0]});
				$ssc =~ tr/ATGCatgc/TACGtacg/;
			}
			else{
				$ssc = $sequences{$tarray[0]};
			}
		}
		else{
			if($tarray[9]-$tarray[8] < 0) {
				$ssc = reverse($sequences{$tarray[0]});
				$ssc =~ tr/ATGCatgc/TACGtacg/;
			}
			else{
				$ssc = $sequences{$tarray[0]};
			}
		}
	}
}
my $ira = reverse($irb);
$ira =~ tr/ATCGatcg/TAGCtagc/;
my $fullcp = $lsc . $irb . $ssc . $ira;
open my $out, ">", $ARGV[2] ."_CP_pieces.txt";
print $out ">lsc\n$lsc\n>irb\n$irb\n>ssc\n$ssc\n";

open my $fullout, ">", "$ARGV[2]\_FULLCP.fsa";
print $fullout ">$ARGV[2]\n$fullcp\n";



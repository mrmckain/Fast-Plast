#!/usr/bin/perl -w
use strict;
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

my %cp_genes;
open my $file, "<", $ARGV[1];
while(<$file>){
	chomp;
	if(/>/){
		$sid = substr($_,1);
	}
	else{
		$sequences{$sid}.=$_;
	}
}
my %gene_orientations;
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
		if($tarray[9]-$tarray[8] > 0) {
				$gene_orientations{ir}{"-"}++;
			}
			else{
				$gene_orientations{ir}{"+"}++;
			}
	}
}

if($gene_orientations{lsc}{"+"} > $gene_orientations{lsc}{"-"}){
	$lsc = reverse($sequences{$lsc_id});
	$lsc =~ tr/ATGCatgc/TACGtacg/;
}
else{
	$lsc = $sequences{lsc_id};
}
if($gene_orientations{ssc}{"+"} > $gene_orientations{ssc}{"-"}){
	$ssc = reverse($sequences{$ssc_id});
	$ssc =~ tr/ATGCatgc/TACGtacg/;
}
else{
	$ssc = $sequences{$ssc_id};
}

#####Think this out on white board.
if($ir_relative_position eq "-"){
	if($gene_orientations{ir}{"+"} > $gene_orientations{ir}{"-"}){
		if($gene_orientations{lsc}{"+"} > $gene_orientations{lsc}{"-"}){
			$irb = reverse($sequences{$ir_id});
			$irb =~ tr/ATGCatgc/TACGtacg/;
		}
		else{

			$irb = reverse($sequences{$ir_id});
			$irb =~ tr/ATGCatgc/TACGtacg/;
		}

	}
	else{
		if($gene_orientations{lsc}{"+"} > $gene_orientations{lsc}{"-"}){
			$irb = reverse($sequences{$ir_id});
			$irb =~ tr/ATGCatgc/TACGtacg/;
		}
		else{

			$irb = reverse($sequences{$ir_id});
			$irb =~ tr/ATGCatgc/TACGtacg/;
		}
					$irb = $sequences{$ir_id};

	}
}

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
	if($tarray[1] eq "rpl23"){
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
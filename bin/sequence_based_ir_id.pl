#!/usr/bin/perl
use strict;
use warnings;
use v5.10; # Minimum perl providing defined-or "//", which allows defaults to be zero (unlike regular "||")

my $DEFAULT_MIN_REGION_LENGTH = 10000;

# Get commandline arguments
my $current_afin        = shift // die 'current_afin (first command-line argument) required';
my $name                = shift // die 'name (second command-line argument) required';
my $sc_region_to_split_on = shift // die 'sc_region_to_split_on (third command-line argument) required'; 
my $min_region_legnth   = shift // $DEFAULT_MIN_REGION_LENGTH;

# Throw error if unknown argument used
shift() && die 'Unexpected fifth argument not supported';

my $sequence;
my $sid;

open my $file, "<", $current_afin;
while(<$file>){
	chomp;
	if(/>/){
		$sid = $_;
	}
	else{
		$sequence.=$_;
	}
}

my $seqlen = length($sequence);
my %regions;
my $current_start;
my $current_end;
my $ir;
my $count=0;

for (my $i=0; $i<$seqlen-151; $i++){
	my $tempseq = substr($sequence, $i, 151);
	my $rcseq = reverse($tempseq);
	$rcseq =~ tr/ATCGNatcgn/TAGCNtagcn/;
	if($sequence =~ /$rcseq/){
		if($ir){
			$current_end = $i;
		}
		else{
			if(defined $current_end){
				$regions{$current_start}{$current_end}="sc";#. "_" . $current_start . "-" . $current_end;
				$count++;
				$current_start=$i;
				$current_end=$i;
				$ir=1;
			}
			else{
				$current_start=$i;
				$current_end=$i;
				$ir=1;
			}
		}
	}
	else{
	
		if($ir){
			$regions{$current_start}{$current_end+150}="ir";#"contig_" . $count. "_" . $current_start . "-" . $current_end;
				$count++;
				$current_start=$i+150;
				$current_end=$i;
				$ir=();
		}
		else{
			if(defined $current_end){
				$current_end=$i;
				$ir=();
			
			}
			else{
				$current_start=$i;
				$current_end=$i;
				$ir=();
			}
		}
	}
}
if($ir){
	$regions{$current_start}{$current_end}="ir";
}
else{
	$regions{$current_start}{$current_end}="sc";
}
for my $start (sort {$a<=>$b} keys %regions){
        for my $end (keys %{$regions{$start}}){
                if ($end-$start < 300){
                        delete $regions{$start};
		
                }
	}
}
my $ssc=0;
for my $start (sort {$a<=>$b} keys %regions){
	for my $end (keys %{$regions{$start}}){
		if($regions{$start}{$end} eq "sc"){
			$ssc++;

			# Remove specified region
			if($ssc == $sc_region_to_split_on) {
				delete $regions{$start};
			}
		}
	}
}

my %final_regions;
my $final_start;
my $final_end;
my $final_regionid;
for my $start (sort {$a<=>$b} keys %regions){
        for my $end (keys %{$regions{$start}}){
                if(defined $final_start){
                        if($final_regionid ne $regions{$start}{$end}){
                                $final_regions{$final_start}{$final_end}=$final_regionid;
                                $final_start=$start;
                                $final_end = $end;
                                $final_regionid=$regions{$start}{$end};
                        }
                        else{
                                $final_end=$end;
                        }
                }
                else{
                        $final_start=$start;
                                $final_end = $end;
                                $final_regionid=$regions{$start}{$end};
                }
        }
}

$final_regions{$final_start}{$final_end}=$final_regionid;

open my $out, ">", $name . "_regions_split" . $sc_region_to_split_on . ".fsa";
for my $start (sort {$a<=>$b} keys %final_regions){
        for my $end (keys %{$final_regions{$start}}){
                if($end-$start < $min_region_legnth){
                        next;
                }
                my $tempseq = substr($sequence,$start,($end-$start+1));
                my $tempid;
		if ($final_regions{$start}{$end} eq "ir"){
			$tempid = "ir_" . $start . "-" . $end;
		}
		else{
			$tempid = "sc_" . $start . "-" . $end;
		}
                print $out ">$tempid\n$tempseq\n";
        }
}
	

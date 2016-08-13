#!/usr/bin/perl -w
use strict;

#USAGE: 1: dumpfile 2: assembly 3: kmersize

my %assemkmer;
my $count;
my %finalcount;
my $seq;
my $id;
$ARGV[0] =~ /(.*?)_\d+dump/;
$id=$1;
open my $assem, "<", $ARGV[1];
while(<$assem>){
        chomp;
        unless(/^>/){
		$seq .= $_;
	}
}
                $seq =~ s/\s//g; 
                my $slen = length($seq);
                for (my $i=0; $i <= ($slen-1); $i++){
                        if($i <= ($slen-$ARGV[2])){
                                my $tseq = substr($seq, $i, $ARGV[2]);
				$tseq =~ tr/atcg/ATCG/;
                                $assemkmer{$tseq}{$i}=0;
			        $finalcount{$i}{$tseq}=0;	
                        }

                        else{
                                my $cleft = ($slen-$i);
                                my $tseq = substr($seq, $i, $cleft);
                                $tseq .= substr($seq, 0, ($ARGV[2]-$cleft));
                                $tseq =~ tr/atcg/ATCG/;
				$assemkmer{$tseq}{$i}=0;
				$finalcount{$i}{$tseq}=0;
                        }
                }
        



open my $dfile, "<", $ARGV[0];
while(<$dfile>){
        chomp;
        if(/^>/){
                $count = substr($_,1);
        }
        else{
		my $revcomp = reverse($_);
		$revcomp =~ tr/ATCGatcg/TAGCtagc/;
                if(exists $assemkmer{$_} && !exists $assemkmer{$revcomp}){ #|| exists $assemkmer{$revcomp}){
                	for my $pos (sort keys %{$assemkmer{$_}}){
                		$finalcount{$pos}{$_}+=$count;
                	}
                }
		elsif(exists $assemkmer{$revcomp} && !exists $assemkmer{$_}){
			for my $pos (sort keys %{$assemkmer{$revcomp}}){
                                $finalcount{$pos}{$revcomp}+=$count;
                        }
                }

		elsif( exists $assemkmer{$_} && exists $assemkmer{$revcomp}){
			for my $pos (sort keys %{$assemkmer{$revcomp}}){
                                $finalcount{$pos}{$revcomp}+=$count;
			}
 			for my $pos (sort keys %{$assemkmer{$_}}){
                                $finalcount{$pos}{$_}+=$count;
                        }
		}


        }
}

open my $out, ">", $id . ".coverage_25kmer.txt";
for my $pos (sort {$a <=> $b} keys %finalcount){
	for my $seq (sort keys %{$finalcount{$pos}}){
		print $out "$seq\t$pos\t$finalcount{$pos}{$seq}\n";
	}
}


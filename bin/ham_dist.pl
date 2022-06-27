#!/usr/bin/perl -w
use strict;
print "$ARGV[0]\n";
my $file = $ARGV[0];
my %align;

open my $tfile, "<", $ARGV[0];
my $id;
while(<$tfile>){
    chomp;
    if(/>/){
           $id = substr($_,1);
        }
        elsif($id){
           $align{$id} .= $_;

        }
}

my %distpairs;
my @seq_ids = keys %align;
for (my $i=0; $i < scalar @seq_ids-1; $i++){
        for (my $k = $i+1; $k <= scalar @seq_ids-1; $k++){
                $distpairs{$seq_ids[$i]}{$seq_ids[$k]}=1;
        }
}


for my $pair1 (keys %distpairs){
	for my $pair2 (keys %{$distpairs{$pair1}}){
                my $hamming = ($align{$pair1} ^ $align{$pair2}) =~ tr/\001-\255//;
                my $norm_ham = $hamming/length($align{$pair1});
                #$all_dist{$pair1}{$pair2}=$norm_ham;
                #$all_dist{$pair2}{$pair1}=$norm_ham;
                #$all_taxa{$pair1}=0;
                #$all_taxa{$pair2}=0;
                print "$hamming/" . length($align{$pair1}) . ": $norm_ham\n";
	}
}


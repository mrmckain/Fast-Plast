#!/usr/bin/perl -w
use strict;

####Split VCF Chloroplast file Ignoring MT-CP regions and only within cp genes

####USAGE: ignore_regions plastid_regions vcf_file

my %plastid_regions;

my %regions_to_keep;

my %ignore;
open my $ignore_file, "<", $ARGV[0];
while(<$ignore_file>){
	chomp;
	my @tarray = split /\s+/;

	for (my $i=$tarray[0]; $i <= $tarray[1]; $i++){
		$ignore{$i}=1;
	}
}
my @follow;
for (my $i = 1; $i <= $ARGV[3]; $i++){
	push (@follow, $i);
}
open my $annotation_file, "<", $ARGV[1];
ANNO: while(<$annotation_file>){
		chomp;
		if(/\~\~\~/){
			next;
		}
		else{
			my @tarray = split /\s+/;
			for (my $i = $tarray[1]; $i <= $tarray[2]; $i++){
				if(exists $ignore{$i}){
					delete $plastid_regions{$tarray[0]};
					next ANNO;
				}
				elsif($i >= $ARGV[3]){
					next ANNO;
				}
				else{
					$plastid_regions{$tarray[0]}{$i}=1;
				}
			}
		}
}
my %gene_sample;
for my $gene (keys %plastid_regions){
	for my $pos (sort {$a<=>$b} keys %{$plastid_regions{$gene}}){
		$regions_to_keep{$pos}=$gene;
		$gene_sample{$gene}{$pos}=0;
	}
}
my %all_genes;
my @labels;
my %bases_sample;

open my $vcf_file, "<", $ARGV[2];

while(<$vcf_file>){
	if(/^\#\#/){
		next;
	}
	elsif(/CHROM/){
		@labels = split /\s+/;
	}

	else{
		
		my @tarray = split /\s+/;
		if($tarray[1] >= $ARGV[3]){
			next;
		}
		
		if(exists $regions_to_keep{$tarray[1]}){
			my @bases;
			push(@bases, $tarray[3]);
			if($tarray[4] =~ /,/){
				push(@bases, $tarray[3]);
				push(@bases, split(",", $tarray[4]));
			}
			else{
				if($tarray[4] !~ /\./){
					push(@bases, $tarray[4]);
				}
			}
			if(exists $bases_sample{$tarray[1]} && $tarray[4] eq "."){
				next;
			}	
			for (my $i = 9; $i <= ((scalar @tarray) - 1); $i++){
				$tarray[$i] =~ /(.):(.)/;
				my $one = $1;
				if($tarray[$i] ne "."){
					if($tarray[4] !~ /\./){
						if(length($bases[$one]) > 1 && length($bases[$one]) > length($tarray[3])){
							my @seqarray = split(//,$bases[$one]);
							for (my $j=0; $j<length($bases[$one]); $j++){
								$all_genes{$regions_to_keep{$tarray[1]+$j}}{$labels[$i]}{$tarray[1]+$j}=$seqarray[$j];
								$bases_sample{$tarray[1]+$j}=$seqarray[$j];
							}
						
						}
						elsif(length($bases[$one]) > 1 && length($tarray[3]) > length($bases[$one])){
							my @seqarray = split(//,$bases[$one]);
							for (my $j=0; $j<length($tarray[3]); $j++){
								if($seqarray[$j]){
                                                                	$all_genes{$regions_to_keep{$tarray[1]+$j}}{$labels[$i]}{$tarray[1]+$j}=$seqarray[$j];
                                                                	$bases_sample{$tarray[1]+$j}=$seqarray[$j];
								}
								else{
									$all_genes{$regions_to_keep{$tarray[1]+$j}}{$labels[$i]}{$tarray[1]+$j}=0;
                                                                        $bases_sample{$tarray[1]+$j}=0;
								}
                                                        }
						}
						else{
							$all_genes{$regions_to_keep{$tarray[1]}}{$labels[$i]}{$tarray[1]}=$bases[$one];
                                                	$bases_sample{$tarray[1]}=$bases[$one];
						}
					}
					else{ 
						$all_genes{$regions_to_keep{$tarray[1]}}{$labels[$i]}{$tarray[1]}=$bases[$one];
						$bases_sample{$tarray[1]}=$bases[$one];
					}
				}
			}
		}
		
	}
	

}
my $counter =1;
for my $geneid (keys %all_genes){
	open my $tout, ">", $ARGV[4] . "_cp_from_vcf.fasta";
	
	for my $sample (keys %{$all_genes{$geneid}}){
		print $tout ">$sample\n";
		for my $pos (sort {$a<=>$b} keys%{$gene_sample{$geneid}}){
			if(exists $all_genes{$geneid}{$sample}{$pos}){
				if($all_genes{$geneid}{$sample}{$pos} eq "0"){
					next;
				}
				else{
					print $tout "$all_genes{$geneid}{$sample}{$pos}";
				}
			}
			else{
				print $tout "N";

			}
		}
	}
	print $tout "\n";
}

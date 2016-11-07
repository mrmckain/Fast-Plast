#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use File::Spec;
use File::Which;
use File::Which qw(which where);
use Cwd;
use Cwd 'abs_path';

BEGIN {

    $ENV{FP_HOME} = "$FindBin::Bin";

}


###directories
my $FPROOT = "$FindBin::RealBin";
my $AFIN_DIR = "$FPROOT/afin";
my $COVERAGE_DIR = "$FPROOT/Coverage_Analysis";
my $FPBIN = "$FPROOT/bin";
my $TRIMMOMATIC;
my $BOWTIE2;
my $SPADES;
my $BLAST;


my $help;
my $paired_end1;
my $paired_end2;
my $single_end;
my $name;
my $bowtie_index=$FPBIN . "/Verdant";
my $posgenes;
my $coverage_check;
my $threads = 4;
my $adapters = $FPBIN . "NEB-PE.fa";

GetOptions('help|?' => \$help, "1=s" => \$paired_end1, "2=s" => \$paired_end2, "single=s" => \$single_end, "bowtie_index=s" => \$bowtie_index, "name=s" => \$name, 'coverage_analysis' => \$coverage_check,'positional_genes' => \$posgenes, "threads=i" => \$threads, "adapters=s" => \$adapters)  or pod2usage( { -message => "ERROR: Invalid parameter." } );


if ($help) {
    pod2usage( { -exitstatus => 0 } );
}

if ( !$paired_end1 && !$single_end ) {
    pod2usage( { -message => "ERROR: Missing reads file(s)." } );
}

if (!$paired_end1 && $paired_end2 || !$paired_end2 && $paired_end1){
	pod2usage( { -message => "ERROR: Missing other paired end file." } );
}

if ( !$name ) {
    pod2usage( { -message => "ERROR: Missing sample name." } );
}

### Get full paths for files.  Glob would work for all of them, but it requires perl 5.6+.  Only using it for the ~ calls, just in case. ####
my $datestring = localtime();
my $start_time = time;
print "$datestring\tStarting Fast-Plast.\n"

my @p1_array;
if($paired_end1){
	my @temp_array = split(",", $paired_end1);
	for my $tfile (@temp_array){
		if($tfile =~ /~/){
			$tfile = glob ($tfile);
		}
		my $abs_path = abs_path($tfile);
		push(@p1_array, $abs_path);
	}
}
my @p2_array;
if($paired_end2){ 
	my @temp_array = split(",", $paired_end2);
	for my $tfile (@temp_array){
		if($tfile =~ /~/){
			$tfile = glob ($tfile);
		}
		my $abs_path = abs_path($tfile);
		push(@p2_array, $abs_path);
	}
}
my @s_array; 
if($single_end){
	my @temp_array = split(",", $single_end);
	for my $tfile (@temp_array){
		if($tfile =~ /~/){
			$tfile = glob ($tfile);
		}
		my $abs_path = abs_path($tfile);
		push(@s_array, $abs_path);
	}
}

my $pe_libs = @p1_array;
my $s_libs = @s_array;

if(!$pe_libs){
	$pe_libs = 0;
}
if(!$s_libs){
	$s_libs = 0;
}

print "Assemblying plastome with $s_libs single end libraries and $pe_libs paired end libraries.\n"


###Get read size###

my $current_runtime = time - $start_time;
print "$current_runtime\tDetermining best kmer sizes.\n";

my $maxsize=0;

if(@p1_array){
	for my $file (@p1_array){
		open my $tfile, "<", $file;
		my $count = 0;
		SIZE: while(<$tfile>){
				my $seq = readline($tfile);
				chomp($seq);
				if(length($seq) > $maxsize){
						$maxsize = length($seq);
				}
				$count++;
				if($count == 100){
					last SIZE;
				}
				readline($tfile);
				readline($tfile);
			}
	}
}


if(@p2_array){
	for my $file (@p2_array){
		open my $tfile, "<", $file;
		my $count = 0;
		SIZE: while(<$tfile>){
				my $seq = readline($tfile);
				chomp($seq);
				if(length($seq) > $maxsize){
						$maxsize = length($seq);
				}
				$count++;
				if($count == 100){
					last SIZE;
				}
				readline($tfile);
				readline($tfile);
			}
	}
}

if(@s_array){
	for my $file (@s_array){
		open my $tfile, "<", $file;
		my $count = 0;
		SIZE: while(<$tfile>){
				my $seq = readline($tfile);
				chomp($seq);
				if(length($seq) > $maxsize){
						$maxsize = length($seq);
				}
				$count++;
				if($count == 100){
					last SIZE;
				}
				readline($tfile);
				readline($tfile);
			}
	}
}
##########

###Set K-mer Size###
my $spades_kmer;
if($maxsize >= 140){
	$spades_kmer = "55,87,121";
}
elsif($maxsize >= 100){
	$spades_kmer = "55,69,87";
}
elsif($maxsize >= 80){
	$spades_kmer = "45,57,69";
}
elsif($maxsize >= 50){
	$spades_kmer = "31,37,43";
}
else{
	$spades_kmer = "23,27,31";
}

print "K-mer sizes for SPAdes set at $spades_kmer.\n"
##########

########## Create Directory ###########

mkdir("$name");
chdir("$name");

########## Start Trimmomatic ##########

$current_runtime = time - $start_time;
print "$current_runtime\tStarting read trimming with Trimmomatic.\nUsing $TRIMMOMATIC.\n";

mkdir("Trimmed_Reads");
chdir("Trimmed_Reads");

if(@p1_array){
	for (my $i=0; $i < $pe_libs; $i++){
		my $trim_exec = "java -classpath " . $TRIMMOMATIC . " org.usadellab.trimmomatic.TrimmomaticPE -threads " . $threads . " " . $p1_array[$i] . " " . $p2_array[$i] . " " . $name."_".$i.".trimmed_P1.fq " . $name."_".$i.".trimmed_U1.fq " . $name."_".$i.".trimmed_P2.fq " . $name."_".$i.".trimmed_U2.fq " . "ILLUMINACLIP:".$adapters.":1:30:10 SLIDINGWINDOW:10:20 MINLEN:40";
		system($trim_exec);
	}
	`cat $name\*trimmed_P1.fq > $name.trimmed_P1.fq`;
	`cat $name\*trimmed_P2.fq > $name.trimmed_P2.fq`;
	`cat $name\*trimmed_U*.fq > $name.trimmed_UP.fq`;
}

if(@s_array){
	for (my $i=0; $i < $s_libs; $i++){
		my $trim_exec = "java -classpath " . $TRIMMOMATIC . " org.usadellab.trimmomatic.TrimmomaticSE -threads " . $threads . " " . $s_array[$i] . " " . $name."_".$i.".trimmed_SE.fq " . "ILLUMINACLIP:".$adapters.":1:30:10 SLIDINGWINDOW:10:20 MINLEN:40";
		system($trim_exec);
	}
	`cat $name\*trimmed_SE.fq >> $name.trimmed_UP.fq`;
}

chdir("../");

########## Start Bowtie2 ##########

$current_runtime = time - $start_time;
print "$current_runtime\tStarting read mapping with bowtie2.\nUsing $BOWTIE2.\n";


mkdir("Bowtie_Mapping");
chdir("Bowtie_Mapping");

my $bowtie2_exec = $BOWTIE2 " --very-sensitive-local --al map_hits.fq --al-conc map_pair_hits.fq -p " . $threads . " -x " . $bowtie_index . " -1 ../Trimmed_Reads/" . $name . ".trimmed_P1.fq -2 ../Trimmed_Reads/" . $name . ".trimmed_P2.fq -U ../Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";
system($bowtie2_exec);

chdir("../");

########## Start SPAdes ##########

$current_runtime = time - $start_time;
print "$current_runtime\tStarting initial assembly with SPAdes.\nUsing $SPADES.\n";

mkdir("Spades_Assembly");
chdir("Spades_Assembly");

my $spades_exec = "python " . $SPADES . " -o spades_iter1 -1 ../Bowtie_Mapping/map_pair_hits.1.fq -2 map_pair_hits.2.fq -s map_hits.fq --only-assembler -k " . $spades_kmer . " -t " . $threads;
system($spades_exec);

chdir("../");

########## Start Afin ##########

$current_runtime = time - $start_time;
print "$current_runtime\tStarting improved assembly with afin.\n";

mkdir("Afin_Assembly");
chdir("Afin_Assembly");

`perl $FPBIN/filter_coverage_assembly.pl ../Spades_Assembly/spades_iter1/contigs.fasta`;
my $afin_exec = $AFIN_DIR . "/afin -c filtered_spades_contigs.fsa -r ../Trimmed_Reads/" . $name .".trimmed* -l 100 -f .1 -d 100 -x " $maxsize*0.75 . " -p 20 -i 2 -o ". $name . "_afin";
system($afin_exec);

my %contig_lengths;
my $afin_contig;
my $max_afin=0;
my $min_afin=100000000;

open my $afin_file, "<", $name . "_afin_iter0.fa";
while(<$afin_file){
	chomp;
	if(/>/){
		/len_(\d+)/;
		my $afinlen=$1;
		$contig_lengths{$_}=$afinlen;
		if($max_afin < $afinlen){
			$max_afin = $afinlen;
		}
		if($min_afin > $afinlen){
			$min_afin = $afinlen;
		}
	}
}

my $total_afin_contigs = keys %contig_lengths;

print "After first attempt at afin, there are $total_afin_contigs with a maximum size of $max_afin and a minimum size of $min_afin.\n"



chdir("../");

########## Start Plastome Finishing ##########

$current_runtime = time - $start_time;
print "$current_runtime\tStarting plastome finishing.\nUsing $posgenes for LSC, SSC, and IR identification.\n";

mkdir("Plastome_Finishing");
chdir("Plastome_Finishing");



########## USAGE BELOW ##########
=pod

=head1 Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes

fast-plast.pl

=head1 USAGE

    fast-plast.pl [-1 <paired_end_file1> -2 <paired_end_file2> || -single <singe_end_file>] -name <sample_name> [options] 
or
    fast-plast.pl -help

Required options:
	
	-1 <filenames>	 		File with forward paired-end reads. Multiple files can be designated with a comma-delimited list. Read files should be in matching order with other paired end files.
	-2 <filenames> 	 		File with reverse paired-end reads. Multiple files can be designated with a comma-delimited list. Read files should be in matching order with other paired end files.
	-s <filenames>	 		File with unpaired reads. Multiple files can be designated with a comma-delimited list.

	PAIRED END AND/OR SINGLE END FILES CAN BE PROVIDED SIMULTAENOUSLY.

	-n <sample_name> 		Name for current assembly. We suggest a species name/accession combination as Fast-Plast will use this name as the FASTA ID in the final assembly.

Advanced options:

	--threads				Number of threads used by Fast-Plast.  [Default = 4]
	--adapters				Files of adapters used in making sequencing library. [Default = NEB-PE]
	--bowtie_index   		User defined bowtie2 index for mapping step. Useful for non-angiosperm or unrespesented taxa in the default index. [Default = Verdant]
	--posgenes  	 		User defined genes for identification of single copy/IR regions and orientation. Useful when major rearrangments are present in user plastomes.
	--coverage_analysis     Flag to run the coverage analysis of a final chloroplast assembly.

=head1 DESCRIPTION

Fast-Plast is a pipeline that leverages existing and novel programs to quickly assemble, orient, and verify whole chloroplast genome sequences. For most datasets with sufficient data, Fast-Plast is able to produce a full-length de novo chloroplast genome assembly in approximately 30 minutes with no user mediation. 

Currently, Fast-Plast is written to accomodate Illumina data, though most data types could be used with a few changes.

Fast-Plast uses a de novo assembly approach by combining the de bruijn graph-based method of SPAdes with an iterative seed-based assembly implemented in afin to close gaps of contigs with low coverage. The pipeline then identifies regions from the quadripartite structure of the chloroplast genome, assigns identity, and orders them according to standard convention. A coverage analysis is then conducted to assess the quality of the final assembly. 

=head1 REQUIREMENTS

Fast-Plast requires Trimmomatic, bowtie2, SPAdes, and BLAST+.

If you use the coverage analysis to verify the assembly, then Jellyfish 2 and R will be needed. We highly recommend the coverage analysis to check the Fast-Plast assembly. 

afin requires a c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most *nix systems. See https://github.com/mrmckain/Fast-Plast for help with installation.

Fast-Plast is coded to use 4 threads during the Trimmomatic, bowtie2, SPAdes, and afin steps. This can simply be changed by the user if this number is not available.

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 50 million reads and running them through Fast-Plast.

=head1 INSTALLATION

All required programs should be in the user's path. 

To install afin:

		cd afin
		make





=cut
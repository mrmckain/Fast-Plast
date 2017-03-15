#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use File::Spec;
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
my $TRIMMOMATIC; #path to trimmomatic executable
my $BOWTIE2; #path to bowtie2 executable
my $SPADES; #path to spades executable
my $BLAST; #path to blast executable
my $SSPACE; #path to sspace exectuable
my $BOWTIE1; #path to bowtie1 executable
my $JELLYFISH; #path to jellyfish2 excecutable

$ENV{'PATH'} = '$BOWTIE1';

my $help;
my $paired_end1;
my $paired_end2;
my $single_end;
my $name="Fast-Plast";
my $bowtie_index = "All";
my $posgenes= $FPBIN . "/Angiosperm_Chloroplast_Genes.fsa";
my $coverage_check;
my $min_coverage = 0;
my $threads = 4;
my $adapters = $FPBIN . "/NEB-PE.fa";
my $version;
my $current_version = "Fast-Plast v.1.2.1";
my $user_bowtie;

GetOptions('help|?' => \$help,'version' => \$version, "1=s" => \$paired_end1, "2=s" => \$paired_end2, "single=s" => \$single_end, "bowtie_index=s" => \$bowtie_index, "user_bowtie=s" => \$user_bowtie, "name=s" => \$name, 'coverage_analysis' => \$coverage_check,'positional_genes' => \$posgenes, "threads=i" => \$threads, "min_coverage=i" => \$min_coverage, "adapters=s" => \$adapters)  or pod2usage( { -message => "ERROR: Invalid parameter." } );

if($version) {
	pod2usage( { -verbose => 99, -sections => "VERSION" } );
}

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

if($user_bowtie){
	if ( !glob($user_bowtie."*")) {
    	pod2usage( { -message => "ERROR: User supplied Bowtie2 indices do not exist. Check path." } );
	}
}

### Get full paths for files.  Glob would work for all of them, but it requires perl 5.6+.  Only using it for the ~ calls, just in case. ####
my $datestring = localtime();
my $start_time = time;
open(STDERR, '>', $name.'_results_error.log') or die "Can't open log.\n";
open(STDOUT, '>', $name.'_results_out.log') or die "Can't open log.\n";
open my $LOGFILE, ">", $name."_Fast-Plast_Progress.log" or die "Can't open log.\n";
print $LOGFILE "$datestring\tStarting $current_version.\n";

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

print $LOGFILE "\t\t\t\tAssemblying plastome with $s_libs single end libraries and $pe_libs paired end libraries.\n";


###Get read size###
my $current_runtime = localtime(); 
print $LOGFILE "$current_runtime\tDetermining best kmer sizes.\n";

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

print $LOGFILE "\t\t\t\tK-mer sizes for SPAdes set at $spades_kmer.\n";
##########

########## Create Directory ###########

mkdir("$name");
chdir("$name");
open my $SUMMARY, ">", $name."_Plastome_Summary.txt";
print $SUMMARY "Sample:\t$name\nFast-Plast Version:\t$current_version\n";

########## Start Trimmomatic ##########
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting read trimming with Trimmomatic.\n\t\t\t\tUsing $TRIMMOMATIC.\n";

mkdir("1_Trimmed_Reads");
chdir("1_Trimmed_Reads");
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


if (-e $name.".trimmed_P1.fq"){
	my $se_size = &count_lines($name.".trimmed_P1.fq");
	chomp($se_size);
	$se_size=($se_size/4)*2;
	print $SUMMARY "Total Cleaned Pair-End Reads:\t$se_size\n";
}
if (-e $name.".trimmed_UP.fq"){
	my $se_size = &count_lines($name.".trimmed_UP.fq");
	chomp($se_size);
	$se_size=$se_size/4;
	print $SUMMARY "Total Cleaned Single End Reads:\t$se_size\n";
}
chdir("../");

########## Start Bowtie2 ##########

$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting read mapping with bowtie2.\n\t\t\t\tUsing $BOWTIE2.\n";


mkdir("2_Bowtie_Mapping");
chdir("2_Bowtie_Mapping");

if($user_bowtie){

	$bowtie_index = $user_bowtie;
}
else{
	$bowtie_index= &build_bowtie2_indices($bowtie_index);
}
my $bowtie2_exec;
if(@p1_array){
	$bowtie2_exec = $BOWTIE2 . " --very-sensitive-local --al map_hits.fq --al-conc map_pair_hits.fq -p " . $threads . " -x " . $bowtie_index . " -1 ../1_Trimmed_Reads/" . $name . ".trimmed_P1.fq -2 ../1_Trimmed_Reads/" . $name . ".trimmed_P2.fq -U ../1_Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";
}
else{
	$bowtie2_exec = $BOWTIE2 . " --very-sensitive-local --al map_hits.fq -p " . $threads . " -x " . $bowtie_index . " -U ../1_Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";

}

system($bowtie2_exec);


if (-e "map_pair_hits.1.fq"){
	my $se_size = &count_lines("map_pair_hits.1.fq");
	chomp($se_size);
	$se_size=($se_size/4)*2;
	print $SUMMARY "Total Concordantly Mapped Reads:\t$se_size\n";
}
if (-e "map_hits.fq"){
	my $se_size = &count_lines("map_hits.fq");
	chomp($se_size);
	$se_size=$se_size/4;
	print $SUMMARY "Total Non-concordantly Mapped Reads:\t$se_size\n";
}
chdir("../");

########## Start SPAdes ##########

$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting initial assembly with SPAdes.\nUsing $SPADES.\n";
mkdir("3_Spades_Assembly");
chdir("3_Spades_Assembly");

my $spades_exec;
if(@p1_array){
	$spades_exec = "python " . $SPADES . " -o spades_iter1 -1 ../2_Bowtie_Mapping/map_pair_hits.1.fq -2 ../2_Bowtie_Mapping/map_pair_hits.2.fq -s ../2_Bowtie_Mapping/map_hits.fq --only-assembler -k " . $spades_kmer . " -t " . $threads;
}
else{
	$spades_exec = "python " . $SPADES . " -o spades_iter1 -s ../2_Bowtie_Mapping/map_hits.fq --only-assembler -k " . $spades_kmer . " -t " . $threads;
}
system($spades_exec);
chdir("../");
########## Start Afin ##########
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting improved assembly with afin.\n";

mkdir("4_Afin_Assembly");
chdir("4_Afin_Assembly");

`perl $FPBIN/filter_spades_contigs_weigthed.pl ../3_Spades_Assembly/spades_iter1/contigs.fasta`;

my %temp_filtered;
my $temp_filter_id;
open my $file_filter, "<", "filtered_spades_contigs.fsa";
open my $out_filter, ">", "temp_filtered_spades_contigs.fsa";
while(<$file_filter>){
	chomp;
	if(/>/){
		$temp_filter_id=$_;
	}
	else{
		$temp_filtered{$temp_filter_id}.=$_;
	}
}
my @bases = ("A", "C", "T", "G");

my %repeats;
for my $base (@bases){
	for my $base2 (@bases){
		$repeats{$base.$base2}=0;
		$repeats{$base2.$base}=0;
	}
}

for my $temp_seq (keys %temp_filtered){
	my %temprepeats=%repeats;
	for my $re (keys %repeats){
	
		for (my $i =0; $i <= length($temp_filtered{$temp_seq})-1; $i++){
			my $result = index($temp_filtered{$temp_seq}, $re, $i);
			if($result >= 0){
				$temprepeats{$re}++;
				$i = $result;
			}	
		}
		
	}
	   my $max_repeat=0;
                for my $re_val (keys %temprepeats){
                        if($temprepeats{$re_val} >$max_repeat){
                                $max_repeat = $temprepeats{$re_val};
                        }
                }
	if($max_repeat/(length($temp_filtered{$temp_seq})/2) > 0.5){
		next;
	}
	else{
		print $out_filter "$temp_seq\n$temp_filtered{$temp_seq}\n";
	}
}
rename("temp_filtered_spades_contigs.fsa", "filtered_spades_contigs.fsa");


my $current_afin;
my $extension = $maxsize*0.75;
my ($total_afin_contigs, $max_afin, $min_afin) = &run_afin("150,50,50",100,"20,15,10","2,1,1","filtered_spades_contigs.fsa",$extension);
print $LOGFILE "\t\t\t\tAfter afin, there are $total_afin_contigs contigs with a maximum size of $max_afin and a minimum size of $min_afin.\n";
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tRemoving nested contigs.\n";
my $gotofinish;
if( $total_afin_contigs > 1){
	$current_afin = $name . "_afin_iter2.fa";

	`$BLAST/makeblastdb -in $current_afin -dbtype nucl`;
	my $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db " . $current_afin . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $current_afin . ".blastn";
	`$blast_afin_exec`;

	&remove_nested($current_afin, $current_afin.".blastn");
	
	$total_afin_contigs = &count_contigs($current_afin);

	
	my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
	my %contigs_db_genes = %$contigs_db_genes;
	if ($total_afin_contigs > 1){
		&remove_contamination($current_afin, \%contigs_db_genes);
	}

	$total_afin_contigs = &count_contigs($current_afin);
	($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
	%contigs_db_genes = %$contigs_db_genes;

	if ($total_afin_contigs > 1){
		($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
		%contigs_db_genes = %$contigs_db_genes;
		$percent_recovered_genes=$percent_recovered_genes*100;
		print $LOGFILE "\t\t\t\tChecking coverage of afin output with $total_afin_contigs contigs after contamination removal.\n";
		print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";

	if(@p1_array){
		$current_afin = &scaffolding($current_afin,$name);
		rename($current_afin, $name.".final.scaffolds.fasta");
		$current_afin=$name.".final.scaffolds.fasta";
		$total_afin_contigs = &count_contigs($current_afin);
		if($total_afin_contigs > 1){
				open my $cpcomposition, ">", "Chloroplast_gene_composition_of_final_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
			close ($cpcomposition);
				mkdir("../Final_Assembly");
				rename($current_afin, "../Final_Assembly/".$current_afin);
				rename("Chloroplast_gene_composition_of_final_contigs.txt", "../Final_Assembly/Chloroplast_gene_composition_of_final_contigs.txt");
				chdir("../Final_Assembly");
				my $temppwd = `pwd`;
                                chomp($temppwd);
                                $temppwd .= "/". $current_afin;	
				print $LOGFILE "\t\t\t\tCannot scaffold contigs into a single piece.  Coverage is too low or poorly distributed across plastome. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die;
		}
		else{
			my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
			my %contigs_db_genes = %$contigs_db_genes;
			$percent_recovered_genes=$percent_recovered_genes*100;
			print $LOGFILE "\t\t\t\tChecking coverage of scaffolded contigs with $total_afin_contigs.\n";
			print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
			open my $cpcomposition, ">", "Chloroplast_gene_composition_of_afin_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
	close ($cpcomposition);
		}
	}
	else{
		open my $cpcomposition, ">", "Chloroplast_gene_composition_of_final_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
			close ($cpcomposition);
				mkdir("../Final_Assembly");
				rename($current_afin, "../Final_Assembly/".$current_afin);
				rename("Chloroplast_gene_composition_of_final_contigs.txt", "../Final_Assembly/Chloroplast_gene_composition_of_final_contigs.txt");
				chdir("../Final_Assembly");
				my $temppwd = `pwd`;
                                chomp($temppwd);
                                $temppwd .= "/". $current_afin;	
				print $LOGFILE "\t\t\t\tCannot scaffold contigs into a single piece.  Coverage is too low or poorly distributed across plastome. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die;
	}
	}	

	else{

		my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
		my %contigs_db_genes = %$contigs_db_genes;
		$percent_recovered_genes=$percent_recovered_genes*100;
		print $LOGFILE "\t\t\t\tChecking coverage of afin output with $total_afin_contigs contigs.\n";
		print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
		open my $cpcomposition, ">", "Chloroplast_gene_composition_of_afin_contigs_nested_removed.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
		close ($cpcomposition);

	
	}	




}
    


else{
	$current_afin = $name . "_afin_iter2.fa";
	
	my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
	my %contigs_db_genes = %$contigs_db_genes;
	$percent_recovered_genes=$percent_recovered_genes*100;
	print $LOGFILE "\t\t\t\tChecking coverage of afin output with $total_afin_contigs contigs after contamination removal.\n";
	print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
	open my $cpcomposition, ">", "Chloroplast_gene_composition_of_afin_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
	close ($cpcomposition);
}

chdir("../");

########## Start Plastome Finishing ##########

$current_runtime = localtime(); 
print $LOGFILE "$current_runtime\tStarting plastome finishing.\n\t\t\t\tUsing $posgenes for LSC, SSC, and IR identification.\n";

mkdir("5_Plastome_Finishing");
chdir("5_Plastome_Finishing");

&orientate_plastome("../4_Afin_Assembly/".$current_afin, $name, "../Final_Assembly/"); 
chdir("../");

if(!-d "Final_Assembly"){
	mkdir("Final_Assembly");
        rename("4_Afin_Assembly/".$current_afin, "Final_Assembly/".$current_afin);
        my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery("Final_Assembly/".$current_afin);
		my %contigs_db_genes = %$contigs_db_genes;
		$percent_recovered_genes=$percent_recovered_genes*100;
		print $LOGFILE "\t\t\t\tChecking coverage of final assembly. Final assembly is the last afin iteration.\n";
		print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
		open my $cpcomposition, ">", "Final_Assembly/Chloroplast_gene_composition_of_afin_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
	close ($cpcomposition);
	if( -e "4_Afin_Assembly/Chloroplast_gene_composition_of_afin_contigs.txt"){
		my $temppwd = `pwd`;
                chomp($temppwd);
                $temppwd .= "/Final_Assembly/". $current_afin;
                print $LOGFILE "\t\t\t\tCould not properly orientate the plastome. Either your plastome does not have an IR or there was an issue with the assembly. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die;
	}
	if( -e "4_Afin_Assembly/Chloroplast_gene_composition_of_afin_contigs_nested_removed.txt"){
                my $temppwd = `pwd`;
                chomp($temppwd);
                $temppwd .= "/Final_Assembly/". $current_afin;
                print $LOGFILE "\t\t\t\tCould not properly orientate the plastome. Either your plastome does not have an IR or there was an issue with the assembly. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die;
	}
}
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tAssembly finished. Check Final_Assembly directory for chloroplast assembly and accessory files.\n";

########## Start Coverage Analysis ##########
if($coverage_check){
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting coverage analyses.\n";
my $check_finish = "Final_Assembly/".$name."_FULLCP.fsa";
unless(-e $check_finish){
	print $LOGFILE "\t\t\t\tCannot complete coverage analysis. Full chloroplast genome not complete.";
	die;
}
mkdir("Coverage_Analysis");
chdir("Coverage_Analysis");
my $build_bowtie2_exec = $BOWTIE2 . "-build ../Final_Assembly/" . $name . "_FULLCP.fsa " . $name . "_bowtie";
system($build_bowtie2_exec);

my $cov_bowtie2_exec;
if(@p1_array){
	$cov_bowtie2_exec = $BOWTIE2 . " --very-sensitive-local --quiet --al map_hits.fq --al-conc map_pair_hits.fq -p " . $threads . " -x " . $name ."_bowtie" . " -1 ../1_Trimmed_Reads/" . $name . ".trimmed_P1.fq -2 ../1_Trimmed_Reads/" . $name . ".trimmed_P2.fq -U ../1_Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";
}
else{
	$cov_bowtie2_exec = $BOWTIE2 . " --very-sensitive-local --quiet --al map_hits.fq -p " . $threads . " -x " . $name ."_bowtie" . "  -U ../1_Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";
}

system($cov_bowtie2_exec);

my $jellyfish_count_exec = $JELLYFISH . " count -m 25 -t ". $threads . " -C -s 1G map_*";
system($jellyfish_count_exec);

my $jellyfish_dump_exec = $JELLYFISH . " dump mer_counts.jf > " . $name . "_25dump";
system($jellyfish_dump_exec);

my $window_cov_exec = "perl " . $COVERAGE_DIR . "/new_window_coverage.pl " . $name . "_25dump ../Final_Assembly/" . $name . "_FULLCP.fsa 25";
system($window_cov_exec);

my $rscript_exec = "Rscript " . $COVERAGE_DIR . "/plot_coverage.r " . $name . ".coverage_25kmer.txt ". $name;
system($rscript_exec);

my $check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25 " . $min_coverage;
system($check_cov_exec);

chdir("../");
$current_runtime = localtime(); 
print $LOGFILE "$current_runtime\tCoverage analysis finished.\n";
if(-z "Coverage_Analysis/".$name."_problem_regions_plastid_assembly.txt"){
	print $LOGFILE "\t\t\t\tNo issues with assembly coverage were identified.\n";

	coverage_summary("Final_Assembly/", "Coverage_Analysis/");

}
else{

	print $LOGFILE "\t\t\t\tProblem areas identified for assembly coverage. Attempting to repair assembly.\n";
	print $LOGFILE "********************STARTING REASSEMBLY********************\n";
	print $SUMMARY "\nVALUES BELOW FROM REASSEMBLED PLASTOME\n";
	mkdir("4.5_Reassemble_Low_Coverage");
	chdir("4.5_Reassemble_Low_Coverage");

	my $lc_remove_contigs = &reassemble_low_coverage("../Final_Assembly/" . $name . "_FULLCP.fsa", "../Coverage_Analysis/".$name."_problem_regions_plastid_assembly.txt");
	$current_afin=&afin_wrap($lc_remove_contigs);

	&orientate_plastome($current_afin, $name, "../Final_Assembly_Fixed_Low_Coverage/"); 
	mkdir("../Coverage_Analysis_Reassembly");
	chdir("../Coverage_Analysis_Reassembly");	
	my $build_bowtie2_exec = $BOWTIE2 . "-build ../Final_Assembly_Fixed_Low_Coverage/" . $name . "_FULLCP.fsa " . $name . "_bowtie";
	system($build_bowtie2_exec);

	my $cov_bowtie2_exec;
	if(@p1_array){
		$cov_bowtie2_exec = $BOWTIE2 . " --very-sensitive-local --quiet --al map_hits.fq --al-conc map_pair_hits.fq -p " . $threads . " -x " . $name ."_bowtie" . " -1 ../1_Trimmed_Reads/" . $name . ".trimmed_P1.fq -2 ../1_Trimmed_Reads/" . $name . ".trimmed_P2.fq -U ../1_Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";
	}
	else{
		$cov_bowtie2_exec = $BOWTIE2 . " --very-sensitive-local --quiet --al map_hits.fq -p " . $threads . " -x " . $name ."_bowtie" . "  -U ../1_Trimmed_Reads/" . $name . ".trimmed_UP.fq -S " . $name . ".sam";
	}

	system($cov_bowtie2_exec);

	my $jellyfish_count_exec = $JELLYFISH . " count -m 25 -t ". $threads . " -C -s 1G map_*";
	system($jellyfish_count_exec);

	my $jellyfish_dump_exec = $JELLYFISH . " dump mer_counts.jf > " . $name . "_25dump";
	system($jellyfish_dump_exec);

	my $window_cov_exec = "perl " . $COVERAGE_DIR . "/new_window_coverage.pl " . $name . "_25dump ../Final_Assembly_Fixed_Low_Coverage/" . $name . "_FULLCP.fsa 25";
	system($window_cov_exec);

	my $rscript_exec = "Rscript " . $COVERAGE_DIR . "/plot_coverage.r " . $name . ".coverage_25kmer.txt ". $name;
	system($rscript_exec);

	my $check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25 " . $min_coverage;
	system($check_cov_exec);

	if(-z $name."_problem_regions_plastid_assembly.txt"){
		print $LOGFILE "\t\t\t\tAssembly was successful! No issues with assembly coverage were identified.\n";
		print $LOGFILE "\t\t\t\tNew assembly can be found in Final_Assembly_Fixed_Low_Coverage.\n";

		coverage_summary("../Final_Assembly_Fixed_Low_Coverage/", "../Coverage_Analysis_Reassembly/");

	}	
	else{




	print $LOGFILE "\t\t\t\tProblem areas identified for assembly coverage.  Check $name\_problem_regions_plastid_assembly.txt in the 4.5_Reassemble_Low_Coverage directory.\n";
 }
}
}
else{
	print $LOGFILE "\t\t\t\tCoverage analysis was not selected.  We highly recommend this option to verify the assembly.\n";
}


	
##########
sub scaffolding {
        $current_runtime = localtime();
        print $LOGFILE "$current_runtime\tStarting scaffolding with SSPACE.\n";
        mkdir ("Scaffolding");
        chdir("Scaffolding");
        my @p1_temp = <../../1_Trimmed_Reads/*P1*>;
        open my $lib_out, ">", $name . "_lib.txt";
        my $lib_counter=1;
        for my $p1file (@p1_temp){
                        my $libc ="lib" . $lib_counter;
                        my $short_p1=$p1file;
                        my $short_p2 = $short_p1;
                        $short_p2 =~ s/P1/P2/;
                        print $lib_out "$libc\t$short_p1\t$short_p2\t300\t0.75\tFR\n";
        }
        close $lib_out;
        my $sspace_build_exec = "perl " . $SSPACE ." -l " . $name."_lib.txt -s ../" .$_[0] . " -k 1 -g 1 -b" . $_[1];
        system($sspace_build_exec);

        my $scaffolded_assembly = "Scaffolding/".$name . ".final.scaffolds.fasta";

        chdir ("../");
        return($scaffolded_assembly);

}


##########

sub build_bowtie2_indices {
	
	my $bowtie_index = $_[0];
	open my $bt2_seq, ">", $bowtie_index.".fsa";
	open my $default_gb, "<", $FPBIN."/GenBank_Plastomes";
	my $gbid;
	if($bowtie_index !~ /^all$/i || $bowtie_index !~ /genbank/i ){


		while(<$default_gb>){
			chomp;
			if(/>/){
				if(/$bowtie_index/){
					$gbid=$_;
				}
				else{
					$gbid=();
				}
			}
			elsif($gbid){
				print $bt2_seq "$gbid\n$_\n";
			}
		}
		close $default_gb;
	}
	if(-z $bowtie_index || $bowtie_index =~ /^all$/i || $bowtie_index =~ /genbank/i ){
		print $LOGFILE "\t\t\t\tSamples for $bowtie_index are not in current GenBank plastomes. Using one representative from each order to make bowtie2 indices.\n";

		open my $bt2_seq, ">", $bowtie_index;
		open my $default_gb, "<", $FPBIN."/GenBank_Plastomes";
		$gbid=();
		my %used;
		while(<$default_gb>){
			chomp;
			if(/>/){
				/_(.*?ales)_/;
				my $temp_order = $1;
				if(!exists $used{$temp_order}){
					$gbid=$_;
					$used{$temp_order}=1;
				}
				else{
					$gbid=();
				}
			}
			elsif($gbid){
				print $bt2_seq "$gbid\n$_\n";
			}
		}
		close $default_gb;
	
	}
	else{
		print $LOGFILE "\t\t\t\tSamples for $bowtie_index used to make bowtie2 indices.\n";
	}

	my $build_bowtie2_exec = $BOWTIE2 . "-build " . $bowtie_index.".fsa" . " " . $name . "_bowtie";
	system($build_bowtie2_exec);
	$bowtie_index=$name . "_bowtie";
	return($bowtie_index);
}

##########

sub count_lines {

	my $total_lines;
	open my $lfile, "<", $_[0];
	while(<$lfile>){
		chomp;
		$total_lines++;
	}
	return($total_lines);
}




##########

sub count_contigs {
	my %contig_lengths;
	my $afin_contig;
	my $lmax_afin=0;
	my $lmin_afin=100000000;

	open my $afin_file, "<", $_[0];
	while(<$afin_file>){
		chomp;
		if(/>/){
			if(/len_(\d+)/){
			my $afinlen=$1;
			if($lmax_afin < $afinlen){
				$lmax_afin = $afinlen;
			}
			if($lmin_afin > $afinlen){
				$lmin_afin = $afinlen;
			}
			$contig_lengths{$_}=$afinlen;
			}
			else{
			$contig_lengths{$_}=1;
			}
		}
	}
	my $total_afin_contigs = keys %contig_lengths;
	return ($total_afin_contigs);
}
#########
sub afin_wrap {
	my $current_afin=$_[0];
my $extension = $maxsize*0.75;
my ($total_afin_contigs, $max_afin, $min_afin) = &run_afin("150,50,50",100,"20,15,10","2,1,1",$_[0],$extension);
print $LOGFILE "\t\t\t\tAfter afin, there are $total_afin_contigs contigs with a maximum size of $max_afin and a minimum size of $min_afin.\n";
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tRemoving nested contigs.\n";
my $gotofinish;
if( $total_afin_contigs > 1){
	$current_afin = $name . "_afin_iter2.fa";

	`$BLAST/makeblastdb -in $current_afin -dbtype nucl`;
	my $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db " . $current_afin . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $current_afin . ".blastn";
	`$blast_afin_exec`;

	&remove_nested($current_afin, $current_afin.".blastn");
	
	$total_afin_contigs = &count_contigs($current_afin);

	
	my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
	my %contigs_db_genes = %$contigs_db_genes;
	if ($total_afin_contigs > 1){
		&remove_contamination($current_afin, \%contigs_db_genes);
	}

	$total_afin_contigs = &count_contigs($current_afin);
	($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
	%contigs_db_genes = %$contigs_db_genes;

	if ($total_afin_contigs > 1){
		($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
		%contigs_db_genes = %$contigs_db_genes;
		$percent_recovered_genes=$percent_recovered_genes*100;
		print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";

	if(@p1_array){
		$current_afin = &scaffolding($current_afin,$name);
		rename($current_afin, $name.".final.scaffolds.fasta");
		$current_afin=$name.".final.scaffolds.fasta";
		$total_afin_contigs = &count_contigs($current_afin);
		if($total_afin_contigs > 1){
				open my $cpcomposition, ">", "Chloroplast_gene_composition_of_final_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
			close ($cpcomposition);
				mkdir("../Final_Assembly");
				rename($current_afin, "../Final_Assembly/".$current_afin);
				rename("Chloroplast_gene_composition_of_final_contigs.txt", "../Final_Assembly/Chloroplast_gene_composition_of_final_contigs.txt");
				chdir("../Final_Assembly");
				my $temppwd = `pwd`;
                                chomp($temppwd);
                                $temppwd .= "/". $current_afin;	
				print $LOGFILE "\t\t\t\tCannot scaffold contigs into a single piece.  Coverage is too low or poorly distributed across plastome. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die;
		}
	}
	else{
		open my $cpcomposition, ">", "Chloroplast_gene_composition_of_final_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
			close ($cpcomposition);
				mkdir("../Final_Assembly");
				rename($current_afin, "../Final_Assembly/".$current_afin);
				rename("Chloroplast_gene_composition_of_final_contigs.txt", "../Final_Assembly/Chloroplast_gene_composition_of_final_contigs.txt");
				chdir("../Final_Assembly");
				my $temppwd = `pwd`;
                                chomp($temppwd);
                                $temppwd .= "/". $current_afin;	
				print $LOGFILE "\t\t\t\tCannot scaffold contigs into a single piece.  Coverage is too low or poorly distributed across plastome. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die;
	}
	}	

	else{

		my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
		my %contigs_db_genes = %$contigs_db_genes;
		$percent_recovered_genes=$percent_recovered_genes*100;
		print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
		open my $cpcomposition, ">", "Chloroplast_gene_composition_of_afin_contigs_nested_removed.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
		close ($cpcomposition);
	
	}	




}
    


else{
	$current_afin = $name . "_afin_iter2.fa";
	
	my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
	my %contigs_db_genes = %$contigs_db_genes;
	$percent_recovered_genes=$percent_recovered_genes*100;
	print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
	open my $cpcomposition, ">", "Chloroplast_gene_composition_of_afin_contigs.txt";
			for my $contig_name (sort keys %contigs_db_genes){
				for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
						print $cpcomposition "$contig_name\t$gene_name\n";
				}
			}
	close ($cpcomposition);
}
return($current_afin);
}
#########
sub run_afin {
	###Sub must be given number of iterations, trim length, and number of reads needed to fuse contigs.
	my $extension = $_[5];
	my $afin_exec = $AFIN_DIR . "/afin -c " . $_[4] . " -r ../1_Trimmed_Reads/" . $name .".trimmed* -l " . $_[0] . " -f .1 -d " . $_[1] . " -x " . $extension . " -p " . $_[2] . " -i " . $_[3] ." -o ". $name . "_afin";
	print $LOGFILE "\t\t\t\tUsing command $afin_exec.\n";
	system($afin_exec);

	my %contig_lengths;
	my $afin_contig;
	my $max_afin=0;
	my $min_afin=100000000;

	open my $afin_file, "<", $name . "_afin_iter2.fa";
	while(<$afin_file>){
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
	return ($total_afin_contigs, $max_afin, $min_afin);
}

##########

sub cpgene_recovery {
	$current_runtime = localtime(); 
	
	print $LOGFILE "$current_runtime\tChecking chloroplast gene recovery in contigs.\n";
	my $current_afin = $_[0];
	my %chloroplast_db_genes;
	open my $cpdbgenes, "<", $FPBIN . "/Angiosperm_Chloroplast_Genes.fsa";
	while(<$cpdbgenes>){
		chomp;
		if(/>/){
			/>(.*?)\_/;
			$chloroplast_db_genes{$1}=1;
		}
	}
	`$BLAST/makeblastdb -in $FPBIN/Angiosperm_Chloroplast_Genes.fsa -dbtype nucl`;
	my $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db " . $posgenes . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $current_afin."_positional_genes" . ".blastn";
	`$blast_afin_exec`;
	my $total_chloroplast_db_genes= scalar keys %chloroplast_db_genes;

	my %hit_chloroplast_db_genes;
	my %contigs_db_genes;
	open my $hitcpdbgenes, "<", $current_afin . "_positional_genes" .".blastn";
	while(<$hitcpdbgenes>){
		chomp;
		my @tarray = split /\s+/;
		$tarray[1] =~ /(.*?)\_/;
		$hit_chloroplast_db_genes{$1}=1;
		$contigs_db_genes{$tarray[0]}{$1}=1;
	}
	my $total_hit_chlorplast_db_genes= scalar keys %hit_chloroplast_db_genes;
	my $percent_recovered_genes = $total_hit_chlorplast_db_genes/$total_chloroplast_db_genes;

	
	return ($percent_recovered_genes, \%contigs_db_genes);

	
}

##########

sub remove_nested {
        
        my $temp_contigsfile = $_[0]; #fasta file of contigs that were blasted against themselves
        my $temp_blastfile = $_[1]; #blast file of self blast
        my %delete_contigs;
        my %blast_scores;
          my %temp_seqs;

        my $tesid;
        open my $tsfile, "<", $temp_contigsfile;
        while(<$tsfile>){
                chomp;
                if(/>/){
                        $tesid = substr($_, 1);
                }
                else{
                        $temp_seqs{$tesid}.=$_;
                }
        }

        open my $checkblast, "<", $temp_blastfile;
        while(<$checkblast>){
                chomp;
                my @tarray = split /\s+/;
                if($tarray[0] eq $tarray[1]){
                        next;
                }
                my $len1 = length($temp_seqs{$tarray[0]});
                my $len2 = length($temp_seqs{$tarray[1]});

                my $max = 0;
                my $min = 0;

                if ($len1 > $len2){
                        $max = $len1;
                        $min = $len2;
                }
                else{
                        $max = $len2;
                       $min =  $len1;
                }

                if($max == $len1){
                        next;
                }

                my $qstart;
                my $qstop;
                my $qhit = $tarray[7]-$tarray[6];
                
                if($qhit < 0){
                        $qstart=$tarray[7];
                        $qstop=$tarray[6];
                        $qhit=abs($qhit);
                }
                else{
                        $qstart=$tarray[6];
                        $qstop=$tarray[7];
                }

                my $shit = $tarray[9]-$tarray[8];
                if(exists $blast_scores{$tarray[0]}{$tarray[1]}){
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$blast_scores{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
                elsif(exists $blast_scores{$tarray[1]}{$tarray[0]}){
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$blast_scores{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
                else{
                        for (my $j=0; $j<$min; $j++){
                                @{$blast_scores{$tarray[0]}{$tarray[1]}}[$j]=0;
                        }
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$blast_scores{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
        }
        
        for my $blastid1 (keys %blast_scores){
                for my $blastid2 (keys %{$blast_scores{$blastid1}}){
                        my $overlap=0;
                        for my $i (@{$blast_scores{$blastid1}{$blastid2}}){
                                if($blast_scores{$blastid1}{$blastid2}[$i] > 0){
                                        $overlap ++;
                                }
                        }
                        

                        my $len_min = length($temp_seqs{$blastid1});
                        if ($overlap/$len_min >= 0.9){
                                $delete_contigs{">".$blastid1}=0;
                        }

                }

        }

        open my $afinout, ">", $temp_contigsfile . "_fixed";
        open my $oldafin, "<", $temp_contigsfile;

        my $tempsid;
        while(<$oldafin>){
                chomp;
                if(/>/){
                        if(exists $delete_contigs{$_}){;
                                $tempsid=();
                                next;
                        }
                        else{
                                $tempsid=$_;
                                print $afinout  "$tempsid\n";
                        }
                }
                elsif($tempsid){
                        print $afinout "$_\n";
                }
        }

         my $ftemp = $temp_contigsfile . "_fixed";
        rename($ftemp, $temp_contigsfile);

}

##########
sub remove_contamination{
	
	my $current_seq = $_[0];
	my $contigs_cp_genes = $_[1];
	my %contigs_cp_genes = %$contigs_cp_genes;

	open my $temp_seq, "<", $current_seq;
	my $tsid;
	my %seqlens;
	while(<$temp_seq>){
		chomp;
		if(/>/){
			$tsid = substr($_,1);
		}
		else{
			$seqlens{$tsid}.=$_;
		}
	}

	my %genes_by_contig;
	my %count_per_contig;
	my %delete_contigs;
	my %all_contigs;
	for my $contig_name (sort keys %contigs_cp_genes){
		for my $gene_name (sort keys %{$contigs_cp_genes{$contig_name}}){
			$count_per_contig{$contig_name}++;
			$genes_by_contig{$gene_name}{$contig_name}=1;
			$all_contigs{$contig_name}=1;
		}
	}

	for my $gene_name (keys %genes_by_contig){
		if(keys %{$genes_by_contig{$gene_name}}>1){
			my $max = 0;
			my $maxid;
			for my $contig_name (keys %{$genes_by_contig{$gene_name}}){
				if($count_per_contig{$contig_name} > $max){
					$max = $count_per_contig{$contig_name};
					$maxid = $contig_name;
				}
			}
			for my $contig_name (keys %{$genes_by_contig{$gene_name}}){
				if($contig_name ne $maxid){
					$delete_contigs{">".$contig_name}=0;
				}
			}

		}
	}


	 open my $afinout, ">", $current_seq . "_fixed";
        open my $oldafin, "<", $current_seq;

        my $tempsid;
        while(<$oldafin>){
                chomp;
                if(/>/){
                        if(exists $delete_contigs{$_} || !exists $all_contigs{substr($_,1)}){ #removing contigs that do not have cp genes
                                $tempsid=();
                                next;
                        }
                        else{
                                $tempsid=$_;
                                print $afinout  "$tempsid\n";
                        }
                }
                elsif($tempsid){
                        print $afinout "$_\n";
                }
        }

        `mv $current_seq\_fixed $current_seq`;
}

##########
sub orientate_plastome{
        
        my $current_afin = $_[0];
        my $name = $_[1];
        my $path_to_final = $_[2];

        `$BLAST/makeblastdb -in $posgenes -dbtype nucl`;

        for (my $i = 0; $i <=3; $i++){

        	`perl $FPBIN/sequence_based_ir_id.pl $current_afin $name $i`;
        	my $split_fullname= $name ."_regions_split".$i.".fsa";

        	`$BLAST/makeblastdb -in $split_fullname -dbtype nucl`;
        	my $blast_afin_exec = $BLAST . "blastn -query " . $split_fullname . " -db " . $split_fullname . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $split_fullname . ".blastn";
			system($blast_afin_exec);

			&remove_nested($split_fullname,$split_fullname.".blastn");

			my $chloroplast_pieces = &count_contigs($split_fullname);

			if($chloroplast_pieces < 2){
				next;
			}
			my $exists_ssc=0;
			my $exists_ir=0;

			my %cp_piece_pos;
			open my $exists_file, "<", $split_fullname;
			while(<$exists_file>){
				chomp;
				if(/ir/){
					$exists_ir++;
					/ir_(\d+).(\d+)/;
					$cp_piece_pos{"ir"}{"len"}=$2-$1;
					
				}
				if(/sc/){
					$exists_ssc++;
					/sc_(.+)/;
					$cp_piece_pos{"sc"}{$1}=1;
					
				}
			}

			if($exists_ssc != 2 || $exists_ir != 1){
				next;
			}
			
			my @range;
			
			for my $sc_range (keys %{$cp_piece_pos{"sc"}}){		
					$sc_range =~ /(\d+).(\d+)/;
					push(@range, $1);
					push(@range, $2);
			}

			my $temp_range1 = abs($range[3]-$range[0])-2;
			my $temp_range2 = abs($range[2]-$range[1])-2;
			unless($temp_range1 == $cp_piece_pos{"ir"}{"len"} || $temp_range2 == $cp_piece_pos{"ir"}{"len"}){
				next;
			}
					
			$blast_afin_exec = $BLAST . "blastn -query " . $split_fullname . " -db " . $posgenes . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $split_fullname . "_positional_genes" . ".blastn";
			system($blast_afin_exec);
			my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($split_fullname);
			my %contigs_db_genes = %$contigs_db_genes;
			$percent_recovered_genes=$percent_recovered_genes*100;

			if($percent_recovered_genes > 90){
					`perl $FPBIN/orientate_plastome_v.2.0.pl $split_fullname $split_fullname\_positional_genes.blastn $name`;
        			my $final_seq = $name ."_FULLCP.fsa";
        			$blast_afin_exec = $BLAST . "blastn -query " . $final_seq . " -db " . $posgenes . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $current_afin . "_positional_genes" . ".blastn";
        			system($blast_afin_exec);
        			($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($final_seq);
        			%contigs_db_genes = %$contigs_db_genes;
        			$percent_recovered_genes=$percent_recovered_genes*100;
					print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $final_seq.\n";
        			open my $cpcomposition, ">", "Chloroplast_gene_composition_of_final_chloroplast_sequence.txt";
                    for my $contig_name (sort keys %contigs_db_genes){
                        for my $gene_name (sort keys %{$contigs_db_genes{$contig_name}}){
                            print $cpcomposition "$contig_name\t$gene_name\n";
                        }
                    }
                    close ($cpcomposition);
                    my %pieces;
                    my $pid;

                    my $final_assembly_seq;
                    open my $fin, "<", $final_seq;
                    while(<$fin>){
                    	if(/>/){
                    		next;
                    	}
                    	else{
                    		$final_assembly_seq.=$_;
                    	}
                    }
                    close($fin);

                    my $flen=length($final_assembly_seq);
                    print $SUMMARY "Total Chloroplast Genome Length:\t$flen\n";

                    open my $cppieces, "<", $name."_CP_pieces.fsa";
                    while(<$cppieces>){
                    	chomp;
                    	if(/>/){

                    		$pid=substr($_,1);
                    	}
                    	else{
                    		$pieces{$pid}.=$_;
                    	}
                    }
                    if($pieces{"lsc"}){
                    	my $llen=length($pieces{lsc});
                    	print $SUMMARY "Large Single Copy Size:\t$llen\n";
                    }
                    if($pieces{"irb"}){
                    	my $llen=length($pieces{irb});
                    	print $SUMMARY "Inverted Repeat Size:\t$llen\n";
                    }
                    if($pieces{"ssc"}){
                    	my $llen=length($pieces{ssc});
                    	print $SUMMARY "Small Single Copy Size:\t$llen\n";
                    }
                    %pieces=();
                    close($cppieces);

		    		mkdir("$path_to_final");
                    rename($final_seq, $path_to_final.$final_seq);
                    rename("Chloroplast_gene_composition_of_final_chloroplast_sequence.txt", $path_to_final."Chloroplast_gene_composition_of_final_contigs.txt");
                    rename($name."_CP_pieces.fsa", $path_to_final.$name."_CP_pieces.fsa");
		    return; #bad form to have multiple exits but i haven't worked out the "correct" way yet
			}
		}


}
##########
sub coverage_summary{


	my %cplens;
	my $cpsid;
	open my $cppieces, "<", $_[0].$name."_CP_pieces.fsa";
	while(<$cppieces>){
		chomp;
		if(/>/){
			$cpsid=substr($_,1);
		}
		else{
			$cplens{$cpsid}.=$_;
		}
	}

	my $end_lsc = length($cplens{lsc})-1;
	my $end_ir = length($cplens{irb})+$end_lsc;
	my $end_ssc = length($cplens{ssc})+$end_ir;

	my $count_ir;
	my $count_ssc;
	my $count_lsc;

	open my $covin, "<", $_[1].$name.".coverage_25kmer.txt";
	while(<$covin>){
			chomp;
			my @tarray = split/\s+/;
			if($tarray[1] <= $end_lsc){
				$count_lsc+=$tarray[2];
			}
			if($tarray[1] >$end_lsc && $tarray[1]<=$end_ir){
				$count_ir+=$tarray[2];
			}
			if($tarray[1] >$end_ir && $tarray[1]<=$end_ssc){
				$count_ssc+=$tarray[2];
			}
	}
	my $avg_lsc = $count_lsc/length($cplens{lsc});
	my $avg_ssc = $count_ssc/length($cplens{ssc});
	my $avg_ir = $count_ir/length($cplens{irb});

	print $SUMMARY "Average Large Single Copy Coverage:\t$avg_lsc\nAverage Inverted Repeat Coverage:\t$avg_ir\nAverage Small Single Copy Coverage:\t$avg_ssc\n";

}
##########
sub reassemble_low_coverage{
	my $final_seq = $_[0];
	my $coverage_file = $_[1];
	my %break_points;
	open my $tcov, "<", $coverage_file;
	while(<$tcov>){
		chomp;
		my @tarray = split/\s+/;
		$break_points{$tarray[0]}=$tarray[1];

	}

	my %new_substrings;
	my $cass_seq;
	open my $cassembly, "<", $final_seq;
	while(<$cassembly>){
		chomp;
		if(/>/){
			next;
		}
		else{
			$cass_seq = $_;
		}
	}
	my $c_start;
	my $c_stop;
	for my $starts (sort {$a <=> $b} keys %break_points){
		if(!$c_stop){
			if ($starts != "0"){
				my $temp_break = substr($cass_seq, 0, $starts);
				$new_substrings{"0"}{$starts-1}=$temp_break;
				$c_start = $break_points{$starts}+1;
			}
		}
		else{
			my $temp_break = substr($cass_seq, $c_start, $starts-$c_start);
				$new_substrings{$c_start}{$starts-$c_start}=$temp_break;
		}
	}
		my $temp_break = substr($cass_seq, $c_start);
				$new_substrings{$c_start}{length($cass_seq)-1}=$temp_break;
	
	close($cassembly);
	close($tcov);

	open my $new_sub, ">", $name . "_removed_lowcoverage_contigs.fsa";
	for my $new_pos (sort {$a<=>$b} keys %new_substrings){
		for my $new_end (keys %{$new_substrings{$new_pos}}){
			print $new_sub ">$new_pos\-$new_end\n$new_substrings{$new_pos}{$new_end}\n";
		}
	}

	my $new_seqfile = $name . "_removed_lowcoverage_contigs.fsa";
	return($new_seqfile)





}
        

########## USAGE BELOW ##########
=pod

=head1 Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes

fast-plast.pl

=head1 USAGE

    fast-plast.pl [-1 <paired_end_file1> -2 <paired_end_file2> || -single <singe_end_file>] -name <sample_name> [options] 
or
    fast-plast.pl -help
	-1 <filenames>		File with forward paired-end reads. Multiple files can be designated with a comma-delimited list. 
				Read files should be in matching order with other paired end files.
	-2 <filenames>		File with reverse paired-end reads. Multiple files can be designated with a comma-delimited list. 
				Read files should be in matching order with other paired end files.
	-s <filenames>		File with unpaired reads. Multiple files can be designated with a comma-delimited list.

	PAIRED END AND SINGLE END FILES CAN BE PROVIDED SIMULTAENOUSLY.

	-n <sample_name>	Name for current assembly. We suggest a species name/accession combination as Fast-Plast will use 
				this name as the FASTA ID in the final assembly.

Advanced options:

	--threads		Number of threads used by Fast-Plast.  [Default = 4]
	--adapters		Files of adapters used in making sequencing library. [Default = NEB-PE]
	--bowtie_index		Order for sample to draw references for mapping. If order exists, then all available samples for that order will be used. 
				If order does not exist in default set or the terms "all" or "GenBank" are given, one exemplar from each available order is used 
				to build the Bowtie2 indicies. [default="All"]
	--user_bowtie		User supplied bowtie2 indices. If this option is used, bowtie_index is ignored.
	--posgenes		User defined genes for identification of single copy/IR regions and orientation. Useful when major rearrangments are present in user plastomes.
	--coverage_analysis 	Flag to run the coverage analysis of a final chloroplast assembly.

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

=head1 VERSION

Fast-Plast v.1.2.0





=cut

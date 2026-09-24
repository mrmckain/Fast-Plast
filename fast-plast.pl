#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use File::Spec;
use File::Path qw(remove_tree);
use IO::Handle;
use Cwd;
use Cwd 'abs_path';
use Env qw (PATH);

BEGIN {

    $ENV{FP_HOME} = "$FindBin::Bin";

}

###directories

my $FPROOT = "$FindBin::RealBin";

#-----------------------------------------------------------------------------
# Tool and data resolution
#
# Fast-Plast runs in two layouts:
#   (1) from a git clone  - afin, helper scripts, and reference data live under
#       the repository root; dependencies must be on PATH.
#   (2) from a conda install - every external tool is on PATH (conda guarantees
#       this) and Fast-Plast's own data lives under $PREFIX/share/fast-plast.
#
# External tools are located on PATH at runtime. Any tool can be overridden
# without editing this file by exporting an environment variable, e.g.
#   FP_BLAST=/opt/blast/bin/blastn   FP_SPADES=/opt/spades/bin/spades.py
#-----------------------------------------------------------------------------

# Locate an executable on PATH, honoring an override env var (file or dir value).
sub find_exe {
    my ($name, $envvar) = @_;
    my $override = ($envvar && defined $ENV{$envvar}) ? $ENV{$envvar} : undef;
    if (defined $override && length $override) {
        return $override if -x $override && ! -d $override;
        my $cand = File::Spec->catfile($override, $name);
        return $cand if -x $cand && ! -d $cand;
    }
    for my $dir (File::Spec->path) {
        my $cand = File::Spec->catfile($dir, $name);
        return $cand if -x $cand && ! -d $cand;
    }
    return undef;
}

# Resolve a REQUIRED tool or die with an actionable message.
sub require_exe {
    my ($name, $envvar) = @_;
    my $path = find_exe($name, $envvar);
    return $path if $path;
    die "ERROR: required dependency '$name' was not found on your PATH.\n"
      . "       Install the package that provides it (e.g. via bioconda),\n"
      . "       or set an override: $envvar=/full/path/to/$name, then re-run.\n";
}

# Directory portion of a resolved executable path (no trailing slash).
sub exe_dir {
    my ($exe) = @_;
    my (undef, $dir, undef) = File::Spec->splitpath($exe);
    $dir =~ s{/+$}{};
    return $dir;
}

# Reference data + bundled helper scripts: explicit override, then conda
# layout ($PREFIX/share/fast-plast, i.e. ../share/fast-plast from the driver
# in $PREFIX/bin), then the classic clone layout.
my $FP_SHARE = $ENV{'FASTPLAST_SHARE'};
$FP_SHARE  ||= "$FPROOT/../share/fast-plast" if -d "$FPROOT/../share/fast-plast/bin";
$FP_SHARE  ||= $FPROOT                       if -d "$FPROOT/bin";
$FP_SHARE or die "ERROR: cannot locate Fast-Plast data directory. "
               . "Set FASTPLAST_SHARE=/path/to/fast-plast/share.\n";
my $FPBIN        = "$FP_SHARE/bin";
my $COVERAGE_DIR = $ENV{'FASTPLAST_COVERAGE'} || "$FP_SHARE/Coverage_Analysis";

# External executables are resolved in resolve_tools(), which runs after
# option parsing so that --version and --help work on a machine without the
# dependencies installed. Values preserve the ORIGINAL call-site syntax:
#   $BLAST    BLAST+ bin *directory* WITH a trailing slash (call sites use both
#             "$BLAST/makeblastdb" and $BLAST . "blastn");
#   $AFIN_DIR directory of the afin binary (call sites use "$AFIN_DIR/afin");
#   the others are full executable paths. jellyfish (coverage analysis only)
#   and ragtag (paired-end scaffolding only) are resolved only when needed.
my ($AFIN_DIR, $BLAST, $BOWTIE2, $SPADES, $FASTP, $PIGZ, $JELLYFISH, $RAGTAG);

sub resolve_tools {
	# afin (compiled C++ core). On PATH under conda; compiled in-tree from a clone.
	my $afin_exe = find_exe('afin', 'FP_AFIN') || "$FPROOT/afin/afin";
	-e $afin_exe or die "ERROR: afin binary not found. Compile it (cd afin && make) "
	                  . "or install Fast-Plast via bioconda.\n";
	$AFIN_DIR = exe_dir($afin_exe);

	$BLAST   = exe_dir( require_exe('blastn',  'FP_BLAST') ) . "/";
	$BOWTIE2 = require_exe('bowtie2',   'FP_BOWTIE2');
	$SPADES  = require_exe('spades.py', 'FP_SPADES');
	$FASTP   = require_exe('fastp',     'FP_FASTP');

	# Decompressor for gzipped reads: prefer pigz (faster and multi-member safe),
	# fall back to gzip. Streamed via "<tool> -dc"; pigz also honors -p <threads>.
	$PIGZ = find_exe('pigz', 'FP_PIGZ') || find_exe('gzip', 'FP_GZIP')
	     or die "ERROR: neither pigz nor gzip was found on PATH.\n"
	          . "       Install one, or set FP_PIGZ / FP_GZIP.\n";
}

my $help;
my $paired_end1;
my $paired_end2;
my $single_end;
my $name="Fast-Plast";
my $bowtie_index = "All";
my $posgenes= $FPBIN . "/Angiosperm_Chloroplast_Genes.fsa";
my $coverage_check;
my $min_coverage;
my $threads = 4;
my $adapters = $FPBIN . "/adapters/NEB-PE.fa";
my $version;
# ---------------------------------------------------------------------------
# Single source of truth for the Fast-Plast version. Update this ONE line on
# release; it propagates to the startup log, the Plastome Summary file, and the
# --version flag. (The POD VERSION section defers to --version, so there is no
# second copy to keep in sync.)
# ---------------------------------------------------------------------------
my $FP_VERSION = "1.3.1";
my $current_version = "Fast-Plast v.$FP_VERSION";
my $user_bowtie;
my $clean;
my $subsample;
my $cov_only;
my $min_region_length = 10000;
# Minimum read length after trimming. Left undefined here; once the read length
# has been sampled it defaults to 140 for reads >= 150 bp and to 90% of the
# read length otherwise (a flat 140 silently discarded every read from 75/100 bp
# libraries).
my $min_length_trim;
my $skip;
my $min_filter_spades;
# SPAdes read error correction (BayesHammer) is on by default. It costs little
# on the reduced, plastid-only read set and closes contig breaks caused by
# indel errors in homopolymer runs, which otherwise leave dead ends in the
# k=121 graph. --spades_only_assembler restores the pre-1.3.1 behaviour.
my $spades_only_assembler;
# Optional user-supplied reference plastome for RagTag scaffolding. Overrides
# the automatic best-match selection from the bundled GenBank plastomes.
my $scaffold_reference = $ENV{'FP_REFERENCE'};
GetOptions('help|?' => \$help,'version' => \$version, "1=s" => \$paired_end1, "2=s" => \$paired_end2, "single=s" => \$single_end, "bowtie_index=s" => \$bowtie_index, "user_bowtie=s" => \$user_bowtie, "name=s" => \$name, "clean=s" => \$clean, 'coverage_analysis' => \$coverage_check, 'skip=s' => \$skip, 'posgenes|positional_genes=s' => \$posgenes, "threads=i" => \$threads, "min_coverage=i" => \$min_coverage, "adapters=s" => \$adapters, "subsample=i" => \$subsample, "only_coverage=s" => \$cov_only, "min_region_length=i" => \$min_region_length, "min_length_trim=i" => \$min_length_trim, "min_filter_spades=i" => \$min_filter_spades, "spades_only_assembler" => \$spades_only_assembler, "scaffold_reference=s" => \$scaffold_reference)  or pod2usage( { -message => "ERROR: Invalid parameter." } );
# Resolve the scaffold reference to an absolute path now, before any chdir,
# so the scaffolding step (which runs several directories deep) can find it.
$scaffold_reference = File::Spec->rel2abs($scaffold_reference) if $scaffold_reference;
# Same for a user-supplied adapter file (keywords are resolved later).
$adapters = File::Spec->rel2abs($adapters) unless $adapters =~ /^(?:nextera|truseq|neb)$/i;
# And the assembly given to --only_coverage: it is used from inside
# <name>/Coverage_Analysis, so a relative path would silently break.
if($cov_only){
	$cov_only = File::Spec->rel2abs($cov_only);
	-s $cov_only or pod2usage( { -message => "ERROR: --only_coverage file '$cov_only' not found or empty." } );
}

if($version) {
	print "$current_version\n";
	exit;
}

if ($help) {
    pod2usage( { -exitstatus => 0 } );
}

resolve_tools();

if ( !$paired_end1 && !$single_end ) {
    pod2usage( { -message => "ERROR: Missing reads file(s)." } );
}

if (!$paired_end1 && $paired_end2 || !$paired_end2 && $paired_end1){
	pod2usage( { -message => "ERROR: Missing other paired end file." } );
}

if ( !$name ) {
    pod2usage( { -message => "ERROR: Missing sample name." } );
}
# The name is spliced unquoted into shell commands and file globs throughout.
if ( $name !~ /^[A-Za-z0-9._-]+$/ ) {
    pod2usage( { -message => "ERROR: --name may only contain letters, digits, '.', '_' and '-'." } );
}

if($user_bowtie){
	if ( !glob($user_bowtie."*")) {
    	pod2usage( { -message => "ERROR: User supplied Bowtie2 indices do not exist. Check path." } );
	}
}

if($posgenes ne $FPBIN . "/Angiosperm_Chloroplast_Genes.fsa"){
	$posgenes = File::Spec->rel2abs($posgenes);
	-s $posgenes or pod2usage( { -message => "ERROR: --posgenes file '$posgenes' not found or empty." } );
}

if($skip && $skip ne "trim"){
	pod2usage( { -message => "ERROR: --skip only accepts 'trim'." } );
}
if($clean && $clean ne "light" && $clean ne "deep"){
	pod2usage( { -message => "ERROR: --clean only accepts 'light' or 'deep'." } );
}

# Fail before any work starts if this run needs something that is missing.
# (These checks print to the terminal; STDERR is redirected to a log below.)
if($coverage_check || $cov_only){
	$JELLYFISH = require_exe('jellyfish', 'FP_JELLYFISH');
	find_exe('Rscript', 'FP_RSCRIPT')
	  or warn "WARNING: Rscript not found on PATH; the coverage plot will be skipped.\n";
}
if($paired_end1 && !$cov_only){
	$RAGTAG = require_exe('ragtag.py', 'FP_RAGTAG');
}
if(!$user_bowtie && !$cov_only){
	-s "$FPBIN/GenBank_Plastomes"
	  or die "ERROR: reference plastome database not found at $FPBIN/GenBank_Plastomes.\n"
	       . "       Download it with:  bash $FPBIN/fetch_plastome_db.sh\n"
	       . "       or supply your own bowtie2 index with --user_bowtie.\n";
}

### Get full paths for files.  Glob would work for all of them, but it requires perl 5.6+.  Only using it for the ~ calls, just in case. ####
my $datestring = localtime();
my $start_time = time;
open(STDERR, '>', $name.'_results_error.log') or die "Can't open log.\n";
open(STDOUT, '>', $name.'_results_out.log') or die "Can't open log.\n";
open my $LOGFILE, ">", $name."_Fast-Plast_Progress.log" or die "Can't open log.\n";
$LOGFILE->autoflush(1);   # so `tail -f` on the progress log actually shows progress
STDOUT->autoflush(1);
STDERR->autoflush(1);
print $LOGFILE "$datestring\tStarting $current_version.\n";

# Run an external command, logging it and dying with the exit status on
# failure. Every pipeline step goes through here so a failed tool stops the run
# at the step that failed instead of several steps later on a missing file.
sub run_cmd {
	my ($cmd, $what) = @_;
	$what ||= "command";
	print $LOGFILE "\t\t\t\tRunning: $cmd\n";
	my $rc = system($cmd);
	return if $rc == 0;
	my $why = ($rc == -1) ? "could not be started ($!)"
	        : ($rc & 127)  ? "was killed by signal " . ($rc & 127)
	        :                "exited with status " . ($rc >> 8);
	print $LOGFILE "\t\t******************ERROR: $what $why.\n\t\t\t\tCommand: $cmd\n";
	die "Fast-Plast: $what $why. See $name\_Fast-Plast_Progress.log and $name\_results_error.log.\n";
}

# Like run_cmd, but a failure is logged as a warning and the run continues.
# For steps that have a fallback (RagTag) or are cosmetic (the coverage plot).
sub run_optional {
	my ($cmd, $what) = @_;
	$what ||= "command";
	print $LOGFILE "\t\t\t\tRunning: $cmd\n";
	my $rc = system($cmd);
	return 1 if $rc == 0;
	print $LOGFILE "\t\t\t\tWARNING: $what failed (status " . ($rc >> 8) . "); continuing.\n";
	return 0;
}

# Run a command and return its STDOUT, dying if it fails.
sub capture_cmd {
	my ($cmd, $what) = @_;
	$what ||= "command";
	print $LOGFILE "\t\t\t\tRunning: $cmd\n";
	my $out = `$cmd`;
	if($? != 0){
		print $LOGFILE "\t\t******************ERROR: $what exited with status " . ($? >> 8) . ".\n";
		die "Fast-Plast: $what failed. See $name\_Fast-Plast_Progress.log and $name\_results_error.log.\n";
	}
	return $out;
}

# Normal termination for paths that stop early on purpose (e.g. --only_coverage).
# Previously these used die(), which exited non-zero and confused workflow managers.
sub finish {
	my ($msg) = @_;
	print $LOGFILE "$msg\n" if $msg;
	close $LOGFILE;
	exit 0;
}

my @p1_array;
if($paired_end1){
	my @temp_array = split(",", $paired_end1);
	for my $tfile (@temp_array){
		if($tfile =~ /~/){
			$tfile = glob ($tfile);
		}
		my $abs_path = abs_path($tfile);
		unless(-e $tfile){
			die "$tfile does not exist.";
		}
		unless(-r $tfile){
			die "$tfile is not readable.";
		}
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
		unless(-e $tfile){
			die "$tfile does not exist.";
		}
		unless(-r $tfile){
			die "$tfile is not readable.";
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
		unless(-e $tfile){
			die "$tfile does not exist.";
		}
		unless(-r $tfile){
			die "$tfile is not readable.";
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

# Longest sequence among the first 100 records of each input file. Reads are
# streamed through open_read_stream so gzipped input is decompressed; opening
# a .gz directly (as older versions did) measured compressed bytes and chose
# SPAdes k-mers and the afin extension length from garbage.
sub max_read_length {
	my ($file, $nrecords) = @_;
	my $in = open_read_stream($file);
	my $longest = 0;
	my $count = 0;
	while(my $h = <$in>){
		my $seq = <$in>;
		last unless defined $seq;
		<$in>; <$in>;                       # '+' line and quality line
		chomp($seq);
		$seq =~ s/\r$//;
		$longest = length($seq) if length($seq) > $longest;
		last if ++$count >= $nrecords;
	}
	close $in;
	return $longest;
}

for my $file (@p1_array, @p2_array, @s_array){
	my $len = max_read_length($file, 100);
	$maxsize = $len if $len > $maxsize;
}
if($maxsize == 0){
	print $LOGFILE "\t\t******************ERROR: could not read any sequences from the input files.\n";
	die "Fast-Plast: could not read any sequences from the input read files.\n";
}
print $LOGFILE "\t\t\t\tMaximum read length sampled from input: $maxsize.\n";
##########



###K-mer sizes###
# SPAdes k-mers are chosen AFTER trimming, from the distribution of trimmed
# read lengths (see below), because quality trimming can leave most reads
# well short of the longest raw read: a 151 bp library whose reads mostly
# trim to 100-120 bp contributes almost nothing at k=121. The table maps a
# representative read length to a k-mer set; it is applied to the 25th
# percentile of trimmed lengths so that at least three quarters of the reads
# are longer than the largest k by a useful margin.
my $spades_kmer;
sub kmers_for_length {
	my ($len) = @_;
	return "55,87,121" if $len >= 140;
	return "55,69,87"  if $len >= 100;
	return "45,57,69"  if $len >= 80;
	return "31,37,43"  if $len >= 50;
	return "23,27,31";
}

# Coverage-based choice, used once the reads have been mapped. Each read of
# length L contributes L-k+1 k-mers, so expected k-mer coverage falls with k
# and falls faster when reads have been trimmed short; a fixed length table
# cannot see that. The top k is the largest odd k <= 127 whose expected k-mer
# coverage over the mapped (plastid) reads stays above $floor, which leaves
# headroom for coverage dips and sequencing errors. Below it the ladder uses
# SPAdes' own spacing, 55,77,99,127, keeping 55 as the smallest k so that
# very deep libraries do not tangle at small k on mitochondrial and nuclear
# plastid-like reads. Plastome size is taken as a nominal 150 kb.
# Returns ($kmer_string, $depth, $top_k, $nreads, \%kcov) or () when there
# are too few mapped reads to decide.
sub kmers_from_coverage {
	my ($floor, @files) = @_;
	my $genome = 150000;
	my %kc; my ($n, $bases) = (0, 0);
	for my $file (@files){
		next unless -s $file;
		my $in = open_read_stream($file);
		while(my $h = <$in>){
			my $seq = <$in>;
			last unless defined $seq;
			<$in>; <$in>;
			chomp($seq);
			my $L = length($seq);
			$n++; $bases += $L;
			for(my $k = 33; $k <= 127; $k += 2){
				my $c = $L - $k + 1;
				last if $c <= 0;
				$kc{$k} += $c;
			}
		}
		close $in;
	}
	return () if $n < 1000;
	$kc{$_} = ($kc{$_} || 0) / $genome for keys %kc;
	my $top = 0;
	for(my $k = 33; $k <= 127; $k += 2){
		$top = $k if ($kc{$k} || 0) >= $floor;
	}
	return () if $top < 55;
	my @ladder = grep { $_ <= $top } (55, 77, 99, 127);
	push @ladder, $top if $top > $ladder[-1];
	return (join(",", @ladder), $bases / $genome, $top, $n, \%kc);
}

# Length statistics over the first $nrecords records of each file given:
# returns (n, max, median, 25th percentile).
sub read_length_stats {
	my ($nrecords, @files) = @_;
	my @len;
	for my $file (@files){
		next unless -s $file;
		my $in = open_read_stream($file);
		my $count = 0;
		while(my $h = <$in>){
			my $seq = <$in>;
			last unless defined $seq;
			<$in>; <$in>;
			chomp($seq);
			$seq =~ s/\r$//;
			push @len, length($seq);
			last if ++$count >= $nrecords;
		}
		close $in;
	}
	return (0, 0, 0, 0) unless @len;
	@len = sort { $a <=> $b } @len;
	return (scalar @len, $len[-1], $len[int(@len / 2)], $len[int(@len / 4)]);
}

###Set minimum post-trim read length###
if(!defined $min_length_trim){
	$min_length_trim = ($maxsize >= 150) ? 140 : int($maxsize * 0.9);
	print $LOGFILE "\t\t\t\tMinimum read length after trimming set to $min_length_trim (from a sampled read length of $maxsize).\n";
}
elsif($min_length_trim > $maxsize){
	print $LOGFILE "\t\t******************WARNING: --min_length_trim $min_length_trim exceeds the sampled read length of $maxsize; few or no reads will survive trimming.******************\n";
}
##########

########## Create Directory ###########

mkdir("$name");
chdir("$name");
open my $SUMMARY, ">", $name."_Plastome_Summary.txt";
print $SUMMARY "Sample:\t$name\nFast-Plast Version:\t$current_version\n";



########## Start fastp ##########
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting read trimming with fastp.\n\t\t\t\tUsing $FASTP.\n";

mkdir("1_Trimmed_Reads");
chdir("1_Trimmed_Reads");

###Subsample Data###

my $total_input_files;
if(@p1_array){
	$total_input_files+=scalar @p1_array;
}
if(@p2_array){
	$total_input_files+=scalar @p2_array;
}
if(@s_array){
	$total_input_files+=scalar @s_array;
}


if($subsample){
	my $readsperfile = $subsample/$total_input_files;
	$readsperfile = int($readsperfile);
	if(@p1_array){
		open my $sub_p1, ">", "subset_file1.fq";
		for my $p1f (@p1_array){
			my $in = &open_read_stream($p1f);
			my $count = 0;
			while(my $h = <$in>){
				my $s = <$in>; my $plus = <$in>; my $q = <$in>;
				last unless defined $q;          # truncated final record
				print $sub_p1 $h, $s, $plus, $q;
				last if ++$count >= $readsperfile;
			}
			close $in;
		}
	}
	if(@p2_array){
		open my $sub_p2, ">", "subset_file2.fq";
		for my $p2f (@p2_array){
			my $in = &open_read_stream($p2f);
			my $count = 0;
			while(my $h = <$in>){
				my $s = <$in>; my $plus = <$in>; my $q = <$in>;
				last unless defined $q;
				print $sub_p2 $h, $s, $plus, $q;
				last if ++$count >= $readsperfile;
			}
			close $in;
		}
	}
	if(@s_array){
		open my $sub_s, ">", "subset_files.fq";
		for my $psf (@s_array){
			my $in = &open_read_stream($psf);
			my $count = 0;
			while(my $h = <$in>){
				my $s = <$in>; my $plus = <$in>; my $q = <$in>;
				last unless defined $q;
				print $sub_s $h, $s, $plus, $q;
				last if ++$count >= $readsperfile;
			}
			close $in;
		}
	}

	# Replace only the libraries that were actually subsampled, and update the
	# library counts the trimming loops iterate over. Previously all three
	# arrays were overwritten unconditionally, which invented a nonexistent
	# single-end file for paired-end runs and left the counts stale.
	@p1_array = ("subset_file1.fq") if @p1_array;
	@p2_array = ("subset_file2.fq") if @p2_array;
	@s_array  = ("subset_files.fq") if @s_array;
	$pe_libs = scalar @p1_array;
	$s_libs  = scalar @s_array;
	print $LOGFILE "\t\t\t\tSubsampled $readsperfile reads from each of $total_input_files input files.\n";
}
if($skip && $skip eq "trim"){
		if(@p1_array){
			for my $p1array (@p1_array){
				&append_reads($p1array, "$name.trimmed_P1.fq");
			}
		}
		if(@p2_array){
			for my $p2array (@p2_array){
				&append_reads($p2array, "$name.trimmed_P2.fq");
			}
		}
		if(@s_array){
			for my $sarray (@s_array){
				&append_reads($sarray, "$name.trimmed_UP.fq");
			}
		}
		if (-s $name.".trimmed_P1.fq"){
			my $se_size = &count_lines($name.".trimmed_P1.fq");
			chomp($se_size);
			$se_size=($se_size/4)*2;
			print $SUMMARY "Total Cleaned Pair-End Reads:\t$se_size\n";
		}
		if (-s $name.".trimmed_UP.fq"){
			my $se_size = &count_lines($name.".trimmed_UP.fq");
			chomp($se_size);
			$se_size=$se_size/4;
			print $SUMMARY "Total Cleaned Single End Reads:\t$se_size\n";
		}	
}
else{
# Keyword selection of a bundled adapter set. Anchored: a user path that merely
# contains "neb" or "truseq" must not be replaced.
if($adapters =~ /^nextera$/i){
	$adapters=$FPBIN."/adapters/NexteraPE-PE.fa";
}
elsif($adapters =~ /^truseq$/i){
	$adapters=$FPBIN."/adapters/TruSeq3-PE.fa";
}
elsif($adapters =~ /^NEB$/i){
	$adapters=$FPBIN."/adapters/NEB-PE.fa";
}
else{
	unless(-s $adapters){
		print $LOGFILE "\t\t******************ERROR: adapter file $adapters not found or empty.\n";
		die "Fast-Plast: adapter file $adapters not found or empty.\n";
	}
}
if(@p1_array){
	for (my $i=0; $i < $pe_libs; $i++){
		my $trim_exec = $FASTP
		  . " --in1 " . $p1_array[$i] . " --in2 " . $p2_array[$i]
		  . " --out1 " . $name."_".$i.".trimmed_P1.fq"
		  . " --out2 " . $name."_".$i.".trimmed_P2.fq"
		  . " --unpaired1 " . $name."_".$i.".trimmed_U1.fq"
		  . " --unpaired2 " . $name."_".$i.".trimmed_U2.fq"
		  . " --adapter_fasta " . $adapters
		  . " --detect_adapter_for_pe"
		  . " --cut_right --cut_right_window_size 10 --cut_right_mean_quality 20"
		  . " --length_required " . $min_length_trim
		  . " --overrepresentation_analysis"
		  . " --thread " . $threads
		  . " --json " . $name."_".$i.".fastp.json --html " . $name."_".$i.".fastp.html"
		  . " 2> " . $name."_".$i.".fastp.log";
		run_cmd($trim_exec, "fastp");
	}
	`cat $name\*trimmed_P1.fq > $name.trimmed_P1.fq`;
	`cat $name\*trimmed_P2.fq > $name.trimmed_P2.fq`;
	if(glob("$name*trimmed_U*.fq")){
		`cat $name\*trimmed_U*.fq > $name.trimmed_UP.fq`;
	}
	
}


if(@s_array){
	for (my $i=0; $i < $s_libs; $i++){
		my $trim_exec = $FASTP
		  . " --in1 " . $s_array[$i]
		  . " --out1 " . $name."_".$i.".trimmed_SE.fq"
		  . " --adapter_fasta " . $adapters
		  . " --cut_right --cut_right_window_size 10 --cut_right_mean_quality 20"
		  . " --length_required " . $min_length_trim
		  . " --overrepresentation_analysis"
		  . " --thread " . $threads
		  . " --json " . $name."_".$i.".SE.fastp.json --html " . $name."_".$i.".SE.fastp.html"
		  . " 2> " . $name."_".$i.".SE.fastp.log";
		run_cmd($trim_exec, "fastp");
	}
	`cat $name\*trimmed_SE.fq >> $name.trimmed_UP.fq`;
}


if (-s $name.".trimmed_P1.fq"){
	my $se_size = &count_lines($name.".trimmed_P1.fq");
	chomp($se_size);
	$se_size=($se_size/4)*2;
	print $SUMMARY "Total Cleaned Pair-End Reads:\t$se_size\n";
}
if (-s $name.".trimmed_UP.fq"){
	my $se_size = &count_lines($name.".trimmed_UP.fq");
	chomp($se_size);
	$se_size=$se_size/4;
	print $SUMMARY "Total Cleaned Single End Reads:\t$se_size\n";
}
unlink(glob("$name\_*"));
my @tfile_read = glob("$name.trimmed*");
for my $check_tfile (@tfile_read){
	if(-z $check_tfile){
		unlink($check_tfile);
	}
}
}

# Choose SPAdes k-mers from the trimmed reads (see kmers_for_length above).
{
	my ($n, $max, $median, $p25) = read_length_stats(200000,
		map { "$name.trimmed_$_.fq" } qw(P1 P2 UP));
	if($n){
		$spades_kmer = kmers_for_length($p25);
		print $LOGFILE "\t\t\t\tTrimmed read lengths ($n sampled): longest $max, median $median, 25th percentile $p25.\n";
	}
	else{
		$spades_kmer = kmers_for_length($maxsize);
		print $LOGFILE "\t\t\t\tNo trimmed reads to sample; choosing k-mers from the raw read length of $maxsize.\n";
	}
	print $LOGFILE "\t\t\t\tK-mer sizes for SPAdes set at $spades_kmer.\n";
}
chdir("../");
##########
if($cov_only){
	$current_runtime = localtime();
	print $LOGFILE "$current_runtime\tStarting coverage analyses.\n";
	my $check_finish = $cov_only;
	unless(-e $check_finish){
        	print $LOGFILE "\t\t\t\tCannot complete coverage analysis. File empty or not found.";
        die "Fast-Plast: coverage analysis could not complete (see the run log).\n";
	}
	mkdir("Coverage_Analysis");
	chdir("Coverage_Analysis");
	my $build_bowtie2_exec = $BOWTIE2 . "-build $cov_only " . $name . "_bowtie";
run_cmd($build_bowtie2_exec, "bowtie2-build");

my $read_args = trimmed_read_args("../1_Trimmed_Reads");
unless(defined $read_args){
        print $LOGFILE "Could not find reads to complete Coverage Analysis.  Check 1_Trimmed_Reads.\n";
        die "Could not find reads to complete Coverage Analysis.  Check 1_Trimmed_Reads.\n";
}
my $cov_bowtie2_exec = "$BOWTIE2 --very-sensitive-local --quiet -p $threads -x ${name}_bowtie$read_args -S $name.sam";
run_cmd($cov_bowtie2_exec, "bowtie2");

my $jellyfish_count_exec = $JELLYFISH . " count -m 25 -t ". $threads . " -C -s 1G map_*";
run_cmd($jellyfish_count_exec, "jellyfish count");

my $jellyfish_dump_exec = $JELLYFISH . " dump mer_counts.jf > " . $name . "_25dump";
run_cmd($jellyfish_dump_exec, "jellyfish dump");

my $window_cov_exec = "perl " . $COVERAGE_DIR . "/new_window_coverage.pl " . $name . "_25dump $cov_only  25";
run_cmd($window_cov_exec, "new_window_coverage.pl");

my $rscript_exec = "Rscript " . $COVERAGE_DIR . "/plot_coverage.r " . $name . ".coverage_25kmer.txt ". $name;
run_optional($rscript_exec, "Rscript coverage plot");

my $check_cov_exec;
if(defined $min_coverage){
        $check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25 ".$min_coverage;
}
else{
        $check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25";
}
my $coverage_used = capture_cmd($check_cov_exec, "check_plastid_coverage.pl");
chomp($coverage_used);
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tMinimum coverage of $coverage_used for verifying assembly.\n";
print $SUMMARY "Minimum Coverage Used for Verification: $coverage_used\n";

chdir("../");
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tCoverage analysis finished.\n";
if(-z "Coverage_Analysis/".$name."_problem_regions_plastid_assembly.txt"){
        print $LOGFILE "\t\t\t\tNo issues with assembly coverage were identified.\n";

        coverage_summary($cov_only, "Coverage_Analysis/");
        if($clean){
                        if ($clean eq "light"){
                                unlink(glob("5_Plastome_Finishing/*fsa.n*"));
                                unlink(glob("*/*bt2"));
                                unlink(glob("Coverage_Analysis/*25dump"));
                                unlink(glob("Coverage_Analysis/mer_counts.jf"));
                                unlink(glob("*/*.sam"));
                        }
                        if ($clean eq "deep"){
                                remove_tree("1_Trimmed_Reads");
                                remove_tree("2_Bowtie_Mapping");
                                remove_tree("3_Spades_Assembly");
                                remove_tree("4_Afin_Assembly");
                                remove_tree("5_Plastome_Finishing");
                                unlink(glob("Coverage_Analysis/*25dump"));
                                unlink(glob("Coverage_Analysis/mer_counts.jf"));
                                unlink(glob("*/*bt2"));

                        }
        }
	 finish("Fast-Plast finished.");

}
else{
	print $LOGFILE "\t\t\t\tProblem regions identified with coverage analysis.  Check output.\n";
	finish("Fast-Plast finished.");
}
}
	
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
my $read_args = trimmed_read_args("../1_Trimmed_Reads");
unless(defined $read_args){
        print $LOGFILE "\t\t******************ERROR: No trimmed read files were identified to run SPAdes.  Please check 1_Trimmed_Reads.******************\n";
        die ("No trimmed read files were identified to run SPAdes.  Please check 1_Trimmed_Reads.\n");
}
my $bowtie2_exec = "$BOWTIE2 --very-sensitive-local -p $threads -x $bowtie_index$read_args -S $name.sam";
run_cmd($bowtie2_exec, "bowtie2");


if (-s "map_pair_hits.1.fq"){
	my $se_size = &count_lines("map_pair_hits.1.fq");
	chomp($se_size);
	$se_size=($se_size/4)*2;
	print $SUMMARY "Total Concordantly Mapped Reads:\t$se_size\n";
}
if (-s "map_hits.fq"){
	my $se_size = &count_lines("map_hits.fq");
	chomp($se_size);
	$se_size=$se_size/4;
	print $SUMMARY "Total Non-concordantly Mapped Reads:\t$se_size\n";
}
chdir("../");

########## Start SPAdes ##########

$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting initial assembly with SPAdes.\n\t\t\t\tUsing $SPADES.\n";
mkdir("3_Spades_Assembly");
chdir("3_Spades_Assembly");

# Final k-mer choice from the mapped reads (see kmers_from_coverage). The
# length-based choice made after trimming stands only if this cannot decide.
{
	my ($kmers, $depth, $top, $n, $kc) = kmers_from_coverage(50,
		map { "../2_Bowtie_Mapping/$_" } qw(map_pair_hits.1.fq map_pair_hits.2.fq map_hits.fq));
	if($kmers){
		my $prof = join("  ", map { sprintf("k%d=%.0fx", $_, $kc->{$_} || 0) } (55, 77, 99, 111, 121, 127));
		printf $LOGFILE "\t\t\t\tMapped reads: %d, plastid depth ~%.0fx (150 kb nominal). Expected k-mer coverage: %s.\n", $n, $depth, $prof;
		print $LOGFILE "\t\t\t\tLargest k with >= 50x expected k-mer coverage: $top. K-mer sizes for SPAdes set at $kmers (was $spades_kmer from read length).\n";
		$spades_kmer = $kmers;
	}
	else{
		print $LOGFILE "\t\t\t\tToo few mapped reads to choose k-mers from coverage; keeping $spades_kmer from read length.\n";
	}
}

my $spades_mode = $spades_only_assembler ? " --only-assembler" : "";
print $LOGFILE "\t\t\t\tSPAdes read error correction " . ($spades_only_assembler ? "disabled (--spades_only_assembler)" : "enabled") . ".\n";
my $spades_exec;
if(-s "../2_Bowtie_Mapping/map_pair_hits.1.fq" && -s "../2_Bowtie_Mapping/map_pair_hits.2.fq" && -s "../2_Bowtie_Mapping/map_hits.fq"){
	$spades_exec = $SPADES . " -o spades_iter1 -1 ../2_Bowtie_Mapping/map_pair_hits.1.fq -2 ../2_Bowtie_Mapping/map_pair_hits.2.fq -s ../2_Bowtie_Mapping/map_hits.fq" . $spades_mode . " -k " . $spades_kmer . " -t " . $threads;
}
elsif(-s "../2_Bowtie_Mapping/map_pair_hits.1.fq" && -s "../2_Bowtie_Mapping/map_pair_hits.2.fq" && (! -e "../2_Bowtie_Mapping/map_hits.fq" || -z "../2_Bowtie_Mapping/map_hits.fq")){
	$spades_exec = $SPADES . " -o spades_iter1 -1 ../2_Bowtie_Mapping/map_pair_hits.1.fq -2 ../2_Bowtie_Mapping/map_pair_hits.2.fq" . $spades_mode . " -k " . $spades_kmer . " -t " . $threads;
}
elsif (-s "../2_Bowtie_Mapping/map_hits.fq"){
	$spades_exec = $SPADES . " -o spades_iter1 -s ../2_Bowtie_Mapping/map_hits.fq" . $spades_mode . " -k " . $spades_kmer . " -t " . $threads;
}
else{
	print $LOGFILE "\t\t******************ERROR: No mapped reads files were identified to run SPAdes.  Please check 2_Bowtie_Mapping.******************\n";
	die ("No mapped reads files were identified to run SPAdes.  Please check 2_Bowtie_Mapping.\n");
}
run_cmd($spades_exec, "SPAdes");
chdir("../");
########## Start Afin ##########
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tStarting improved assembly with afin.\n";

mkdir("4_Afin_Assembly");
chdir("4_Afin_Assembly");
if($min_filter_spades){
	run_cmd("perl $FPBIN/filter_spades_contigs_weigthed.pl ../3_Spades_Assembly/spades_iter1/contigs.fasta $min_filter_spades", "filter_spades_contigs_weigthed.pl");
}
else{
	run_cmd("perl $FPBIN/filter_spades_contigs_weigthed.pl ../3_Spades_Assembly/spades_iter1/contigs.fasta", "filter_spades_contigs_weigthed.pl");
}

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
my ($total_afin_contigs, $max_afin, $min_afin) = &run_afin(10,100,20,2,"filtered_spades_contigs.fsa",$extension);
print $LOGFILE "\t\t\t\tAfter afin, there are $total_afin_contigs contigs with a maximum size of $max_afin and a minimum size of $min_afin.\n";
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tRemoving nested contigs.\n";
my $gotofinish;
if( $total_afin_contigs > 1){
	$current_afin = $name . "_afin_iter0.fa";

	run_cmd("$BLAST/makeblastdb -in $current_afin -dbtype nucl", "makeblastdb");
	my $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db " . $current_afin . " -evalue 1e-40 -outfmt 6 -max_target_seqs 100000000 > " . $current_afin . ".blastn";
	run_cmd($blast_afin_exec, "blastn");

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
		$current_afin=&afin_wrap($current_afin, "extend");		

		run_cmd("$BLAST/makeblastdb -in $current_afin -dbtype nucl", "makeblastdb");
		 $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db " . $current_afin . " -evalue 1e-40 -outfmt 6 -max_target_seqs 100000000 > " . $current_afin . ".blastn";
		run_cmd($blast_afin_exec, "blastn");

		&remove_nested($current_afin, $current_afin.".blastn");
		 ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
		 %contigs_db_genes = %$contigs_db_genes;
		 &remove_contamination($current_afin, \%contigs_db_genes);
		 print $LOGFILE "\t\t\t\tRemoved mitochondrial-like sequences.\n";
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
				die "Fast-Plast: could not scaffold the contigs into a single plastome (coverage too low or uneven, or the reference was too distant). Best contigs are in Final_Assembly/. See the run log for details.\n";
		}
		else{
			my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($current_afin);
			my %contigs_db_genes = %$contigs_db_genes;
			$percent_recovered_genes=$percent_recovered_genes*100;
			print $LOGFILE "\t\t\t\tChecking coverage of scaffolded contigs with $total_afin_contigs.\n";
			print $LOGFILE "\t\t\t\t$percent_recovered_genes\% of known angiosperm chloroplast genes were recovered in $current_afin.\n";
			open my $cpcomposition, ">", "Chloroplast_gene_composition_of_scaffolded_contigs.txt";
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
				die "Fast-Plast: could not scaffold the contigs into a single plastome (coverage too low or uneven, or the reference was too distant). Best contigs are in Final_Assembly/. See the run log for details.\n";
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
	$current_afin = $name . "_afin_iter0.fa";
	
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
				die "Fast-Plast: could not orient the plastome (no detectable IR, or an assembly issue). Best contigs are in Final_Assembly/. See the run log for details.\n";
	}
	if( -e "4_Afin_Assembly/Chloroplast_gene_composition_of_afin_contigs_nested_removed.txt"){
                my $temppwd = `pwd`;
                chomp($temppwd);
                $temppwd .= "/Final_Assembly/". $current_afin;
                print $LOGFILE "\t\t\t\tCould not properly orientate the plastome. Either your plastome does not have an IR or there was an issue with the assembly. Best contigs are in $temppwd\. A list of genes in each contig can be found in \"Chloroplast_gene_composition_of_final_contigs.txt\"\.\n";
				die "Fast-Plast: could not orient the plastome (no detectable IR, or an assembly issue). Best contigs are in Final_Assembly/. See the run log for details.\n";
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
	die "Fast-Plast: coverage analysis could not complete (see the run log).\n";
}
mkdir("Coverage_Analysis");
chdir("Coverage_Analysis");
my $build_bowtie2_exec = $BOWTIE2 . "-build ../Final_Assembly/" . $name . "_FULLCP.fsa " . $name . "_bowtie";
run_cmd($build_bowtie2_exec, "bowtie2-build");

my $read_args = trimmed_read_args("../1_Trimmed_Reads");
unless(defined $read_args){
        print $LOGFILE "Could not find reads to complete Coverage Analysis.  Check 1_Trimmed_Reads.\n";
        die "Could not find reads to complete Coverage Analysis.  Check 1_Trimmed_Reads.\n";
}
my $cov_bowtie2_exec = "$BOWTIE2 --very-sensitive-local --quiet -p $threads -x ${name}_bowtie$read_args -S $name.sam";
run_cmd($cov_bowtie2_exec, "bowtie2");
my $jellyfish_count_exec = $JELLYFISH . " count -m 25 -t ". $threads . " -C -s 1G " . "map_*";
run_cmd($jellyfish_count_exec, "jellyfish count");

my $jellyfish_dump_exec = $JELLYFISH . " dump mer_counts.jf > " . $name . "_25dump";
run_cmd($jellyfish_dump_exec, "jellyfish dump");

my $window_cov_exec = "perl " . $COVERAGE_DIR . "/new_window_coverage.pl " . $name . "_25dump ../Final_Assembly/" . $name . "_FULLCP.fsa 25";
run_cmd($window_cov_exec, "new_window_coverage.pl");

my $rscript_exec = "Rscript " . $COVERAGE_DIR . "/plot_coverage.r " . $name . ".coverage_25kmer.txt ". $name;
run_optional($rscript_exec, "Rscript coverage plot");

my $check_cov_exec;
if(defined $min_coverage){
	$check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25 ".$min_coverage;
}
else{
	$check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25";
}
my $coverage_used = capture_cmd($check_cov_exec, "check_plastid_coverage.pl");
chomp($coverage_used);
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tMinimum coverage of $coverage_used for verifying assembly.\n";
print $SUMMARY "Minimum Coverage Used for Verification: $coverage_used\n";
chdir("../");
$current_runtime = localtime(); 
print $LOGFILE "$current_runtime\tCoverage analysis finished.\n";
if(-z "Coverage_Analysis/".$name."_problem_regions_plastid_assembly.txt"){
	print $LOGFILE "\t\t\t\tNo issues with assembly coverage were identified.\n";

	coverage_summary("Final_Assembly/${name}_CP_pieces.fsa", "Coverage_Analysis/");
	if($clean){
			if ($clean eq "light"){
				unlink(glob("5_Plastome_Finishing/*fsa.n*"));
				unlink(glob("*/*bt2"));
				unlink(glob("Coverage_Analysis/*25dump"));
				unlink(glob("Coverage_Analysis/mer_counts.jf"));
				unlink(glob("*/*.sam"));
			}
			if ($clean eq "deep"){
				remove_tree("1_Trimmed_Reads");
				remove_tree("2_Bowtie_Mapping");
				remove_tree("3_Spades_Assembly");
				remove_tree("4_Afin_Assembly");
				remove_tree("5_Plastome_Finishing");
				unlink(glob("Coverage_Analysis/*25dump"));
				unlink(glob("Coverage_Analysis/mer_counts.jf"));
				unlink(glob("*/*bt2"));

			}
	}

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
	run_cmd($build_bowtie2_exec, "bowtie2-build");

	my $read_args = trimmed_read_args("../1_Trimmed_Reads");
	unless(defined $read_args){
	        print $LOGFILE "Could not find reads to complete Coverage Analysis.  Check 1_Trimmed_Reads.\n";
	        die "Could not find reads to complete Coverage Analysis.  Check 1_Trimmed_Reads.\n";
	}
	my $cov_bowtie2_exec = "$BOWTIE2 --very-sensitive-local --quiet -p $threads -x ${name}_bowtie$read_args -S $name.sam";
	run_cmd($cov_bowtie2_exec, "bowtie2");
	my $jellyfish_count_exec = $JELLYFISH . " count -m 25 -t ". $threads . " -C -s 1G map_*";
	run_cmd($jellyfish_count_exec, "jellyfish count");

	my $jellyfish_dump_exec = $JELLYFISH . " dump mer_counts.jf > " . $name . "_25dump";
	run_cmd($jellyfish_dump_exec, "jellyfish dump");

	my $window_cov_exec = "perl " . $COVERAGE_DIR . "/new_window_coverage.pl " . $name . "_25dump ../Final_Assembly_Fixed_Low_Coverage/" . $name . "_FULLCP.fsa 25";
	run_cmd($window_cov_exec, "new_window_coverage.pl");

	my $rscript_exec = "Rscript " . $COVERAGE_DIR . "/plot_coverage.r " . $name . ".coverage_25kmer.txt ". $name;
	run_optional($rscript_exec, "Rscript coverage plot");
	if(defined $min_coverage){
        	$check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25 ".$min_coverage;
	}
	else{
        	$check_cov_exec = "perl " . $COVERAGE_DIR . "/check_plastid_coverage.pl " . $name . ".coverage_25kmer.txt 25";
	}
	
	my $coverage_used = capture_cmd($check_cov_exec, "check_plastid_coverage.pl");
	chomp($coverage_used);
	$current_runtime = localtime();
	print $LOGFILE "$current_runtime\tMinimum coverage of $coverage_used for verifying assembly.\n";
	print $SUMMARY "Minimum Coverage Used for Verification: $coverage_used\n";
	if(-z $name."_problem_regions_plastid_assembly.txt"){
		print $LOGFILE "\t\t\t\tAssembly was successful! No issues with assembly coverage were identified.\n";
		print $LOGFILE "\t\t\t\tNew assembly can be found in Final_Assembly_Fixed_Low_Coverage.\n";

		coverage_summary("../Final_Assembly_Fixed_Low_Coverage/${name}_CP_pieces.fsa", "../Coverage_Analysis_Reassembly/");
		chdir("../");
		if($clean){
			if ($clean eq "light"){
				unlink(glob("5_Plastome_Finishing/*fsa.n*"));
				unlink(glob("4.5_Reassemble_Low_Coverage/*fsa.n*"));
				unlink(glob("*/*bt2"));
				unlink(glob("Coverage_Analysis/*25dump"));
				unlink(glob("Coverage_Analysis/mer_counts.jf"));
				unlink(glob("*/*.sam"));
				unlink(glob("Coverage_Analysis_Reassembly/mer_counts.jf"));
				unlink(glob("Coverage_Analysis_Reassembly/*25dump"));

			}
			if ($clean eq "deep"){
				remove_tree("1_Trimmed_Reads");
				remove_tree("2_Bowtie_Mapping");
				remove_tree("3_Spades_Assembly");
				remove_tree("4_Afin_Assembly");
				remove_tree("5_Plastome_Finishing");
				unlink(glob("Coverage_Analysis/*25dump"));
				unlink(glob("Coverage_Analysis/mer_counts.jf"));
				unlink(glob("*/*bt2"));
				unlink(glob("Coverage_Analysis_Reassembly/mer_counts.jf"));
				unlink(glob("Coverage_Analysis_Reassembly/*25dump"));
				remove_tree("4.5_Reassemble_Low_Coverage");

			}
	}

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
        my ($contigs, $name) = @_;
        $current_runtime = localtime();
        print $LOGFILE "$current_runtime\tStarting reference-guided scaffolding with RagTag.\n";
        mkdir ("Scaffolding");
        chdir("Scaffolding");

        # The returned path is relative to the PARENT directory (the caller
        # chdir's back up before using it), matching the original contract.
        my $result;

        # Reference selection: honor a user-supplied reference if given
        # (-scaffold_reference / FP_REFERENCE), otherwise auto-pick the best
        # match from the bundled GenBank plastomes.
        my $ref;
        if($scaffold_reference){
                if(-s $scaffold_reference){
                        run_cmd("cp \"$scaffold_reference\" scaffold_reference.fsa", "copy of scaffold reference");
                        $ref = "scaffold_reference.fsa";
                        print $LOGFILE "$current_runtime\tUsing user-supplied scaffold reference: $scaffold_reference\n";
                }
                else{
                        print $LOGFILE "$current_runtime\tWARNING: -scaffold_reference '$scaffold_reference' not found or empty; falling back to automatic selection.\n";
                }
        }
        $ref ||= &select_reference_plastome("../$contigs", "scaffold_reference.fsa");

        if($ref){
                my $ragtag_exec = "$RAGTAG scaffold $ref ../$contigs -o ragtag_out -w -r -t $threads";
                run_optional($ragtag_exec, "RagTag");
                if(-e "ragtag_out/ragtag.scaffold.fasta"){
                        $result = "Scaffolding/ragtag_out/ragtag.scaffold.fasta";
                }
        }

        # Graceful fallback: no usable reference, or RagTag made no joins. Pass
        # the afin contigs through unscaffolded so the pipeline still completes.
        unless($result){
                $current_runtime = localtime();
                print $LOGFILE "$current_runtime\tScaffolding produced no joins; using afin contigs unscaffolded.\n";
                run_cmd("cp ../$contigs ./unscaffolded.fasta", "copy of unscaffolded contigs");
                $result = "Scaffolding/unscaffolded.fasta";
        }

        chdir ("../");
        return($result);

}

##########

sub select_reference_plastome {
        # $query is a path relative to the current (Scaffolding) directory.
        # Writes the chosen reference to $outref and returns it, or undef.
        my ($query, $outref) = @_;
        my $gb = $FPBIN . "/GenBank_Plastomes";
        return undef unless (-e $gb && -e $query);

        # The BLAST db of the whole reference set is large and slow to build, so
        # build it once beside the reference file when that directory is
        # writable and reuse it across runs; otherwise build in the (writable)
        # CWD. A lock directory guards against two concurrent runs building the
        # shared copy at the same time: the second run just builds locally.
        # NOTE: no -parse_seqids -- the bundled multi-plastome file has headers
        # that make strict seqid parsing abort, and we don't need blastdbcmd
        # lookups because we extract the chosen record by FASTA scan.
        my $dbprefix = "gb_ref_db";
        my $shared   = "$FPBIN/GenBank_Plastomes.blastdb";
        if(-e "$shared.nal" || -e "$shared.nsq"){
                $dbprefix = $shared;
        }
        elsif(-w $FPBIN && mkdir("$shared.lock")){
                my $ok = run_optional("$BLAST/makeblastdb -in $gb -dbtype nucl -out $shared > makeblastdb.log 2>&1", "makeblastdb (reference plastomes)");
                rmdir("$shared.lock");
                $dbprefix = $shared if $ok;
        }
        if($dbprefix eq "gb_ref_db"){
                run_optional("$BLAST/makeblastdb -in $gb -dbtype nucl -out $dbprefix > makeblastdb.log 2>&1", "makeblastdb (reference plastomes)")
                        or return undef;
        }

        my $blout = "reference_selection.blastn";
        run_optional("$BLAST/blastn -query $query -db $dbprefix -outfmt \"6 sseqid bitscore\" -max_target_seqs 50 -evalue 1e-20 > $blout 2>/dev/null", "blastn (reference selection)")
                or return undef;

        # Pick the reference with the highest summed bitscore across all contigs.
        my %score;
        open my $bl, "<", $blout or return undef;
        while(<$bl>){
                chomp;
                my ($sid, $bits) = split /\t/;
                next unless defined $bits;
                $score{$sid} += $bits;
        }
        close $bl;
        return undef unless %score;

        my ($best) = sort { $score{$b} <=> $score{$a} } keys %score;
        return &extract_fasta_record($gb, $best, $outref) ? $outref : undef;
}

# ---------------------------------------------------------------------------
# afin strand-switch check.
#
# afin can fuse two contigs through a short repeat (for example a copy of
# IR-boundary sequence elsewhere in the LSC): at the repeat it follows the
# branch with the most reads, which is the double-coverage IR side, and
# everything downstream of the join ends up inverted. The pipeline then has a
# single, wrong contig and nothing downstream looks at it again. So after
# every afin run the contigs are placed on a reference plastome and any
# contig whose segments map to both strands is flagged.
#
# Two legitimate two-strand cases are excluded: the IR, whose query interval
# maps to both strands by definition, and the SSC, which is found in either
# orientation in a population of plastomes and is reported as a note, not a
# problem, when the reversed segment's reference coordinates fall inside the
# reference SSC (the shorter stretch between the two IR copies).
#
# The reference is --scaffold_reference if given, otherwise the best BLAST
# match among the bundled plastomes (chosen once and reused). Findings go to
# the progress log, the summary file and <contigs>.strand_check.txt. Returns
# the number of flagged contigs; the run continues either way.
# ---------------------------------------------------------------------------
my $strand_check_ref;
sub check_strand_switches {
	my ($contigs) = @_;
	return 0 unless -s $contigs;
	my $ref;
	if($scaffold_reference && -s $scaffold_reference){ $ref = $scaffold_reference; }
	elsif($strand_check_ref && -s $strand_check_ref){ $ref = $strand_check_ref; }
	else{
		my $picked = &select_reference_plastome($contigs, "strand_check_reference.fsa");
		unless($picked){
			print $LOGFILE "\t\t\t\tStrand check skipped: no reference plastome could be selected.\n";
			return 0;
		}
		$strand_check_ref = File::Spec->rel2abs($picked);
		$ref = $strand_check_ref;
	}
	my $db = "strand_check_ref_db";
	run_optional("$BLAST/makeblastdb -in \"$ref\" -dbtype nucl -out $db > /dev/null 2>&1", "makeblastdb (strand check)") or return 0;
	my $bl = "$contigs.strand_check.blastn";
	run_optional("$BLAST/blastn -query $contigs -db $db -outfmt \"6 qseqid qlen qstart qend sstart send length slen\" -evalue 1e-50 > $bl 2>/dev/null", "blastn (strand check)") or return 0;

	my %hits; my $reflen = 0;
	open my $in, "<", $bl or return 0;
	while(<$in>){
		chomp; my @c = split /\t/;
		next if $c[6] < 2000;                      # only substantial alignments
		my ($qs, $qe, $ss, $se) = @c[2..5];
		my $str = ($ss <= $se) ? "+" : "-";
		($ss, $se) = ($se, $ss) if $ss > $se;
		push @{$hits{$c[0]}}, { qs => $qs, qe => $qe, ss => $ss, se => $se, str => $str };
		$reflen = $c[7];
	}
	close $in;

	my (@flags, @notes);
	for my $q (sort keys %hits){
		my @hs = @{$hits{$q}};
		# IR-like HSPs: the same query interval is also hit on the opposite strand
		my @ir_ref;
		for my $x (@hs){
			for my $y (@hs){
				next if $x == $y || $x->{str} eq $y->{str};
				my $lo = ($x->{qs} > $y->{qs}) ? $x->{qs} : $y->{qs};
				my $hi = ($x->{qe} < $y->{qe}) ? $x->{qe} : $y->{qe};
				my $short = (($x->{qe} - $x->{qs}) < ($y->{qe} - $y->{qs})) ? ($x->{qe} - $x->{qs}) : ($y->{qe} - $y->{qs});
				if($short > 0 && ($hi - $lo) / $short >= 0.8){
					$x->{ir} = 1;
					push @ir_ref, [$x->{ss}, $x->{se}];
				}
			}
		}
		# reference SSC: the shorter of the two stretches between the IR copies
		my @ssc;    # list of [lo, hi] intervals on the linear reference
		if(@ir_ref >= 2){
			my @c = sort { $a->[0] <=> $b->[0] } @ir_ref;
			my ($first, $last) = ($c[0], $c[-1]);
			my $between = $last->[0] - $first->[1];
			my $around  = ($reflen - $last->[1]) + $first->[0];
			@ssc = ($between <= $around) ? ([$first->[1], $last->[0]]) : ([$last->[1], $reflen], [1, $first->[0]]);
		}
		my @sc = sort { $a->{qs} <=> $b->{qs} } grep { !$_->{ir} } @hs;
		# bases covered per strand, as the union of intervals (HSPs can overlap)
		my %len;
		for my $str ("+", "-"){
			my ($cov, $lo, $hi) = (0);
			for my $h (sort { $a->{qs} <=> $b->{qs} } grep { $_->{str} eq $str } @sc){
				if(defined $hi && $h->{qs} <= $hi){ $hi = $h->{qe} if $h->{qe} > $hi; }
				else{ $cov += $hi - $lo + 1 if defined $hi; ($lo, $hi) = ($h->{qs}, $h->{qe}); }
			}
			$cov += $hi - $lo + 1 if defined $hi;
			$len{$str} = $cov;
		}
		next unless $len{"+"} && $len{"-"};
		my $minor = ($len{"+"} >= $len{"-"}) ? "-" : "+";
		my @seg = grep { $_->{str} eq $minor } @sc;
		my $in_ssc = 0;
		if(@ssc){
			$in_ssc = 1;
			for my $s (@seg){
				my $ok = grep { $s->{ss} >= $_->[0] - 500 && $s->{se} <= $_->[1] + 500 } @ssc;
				$in_ssc = 0 unless $ok;
			}
		}
		my $desc = join("; ", map { sprintf("contig %d-%d = reference %d-%d", $_->{qs}, $_->{qe}, $_->{ss}, $_->{se}) } @seg);
		if($in_ssc){
			push @notes, "$q: the SSC is in the opposite orientation to the reference ($desc). This is a normal isomer, not an error.";
		}
		else{
			push @flags, "$q: reversed segment(s) $desc, while the rest of the contig maps to the other strand. Likely an afin misjoin at a repeat; check this contig before trusting the assembly.";
		}
	}

	open my $rep, ">", "$contigs.strand_check.txt";
	print $rep "Reference: $ref\n";
	print $rep "Contigs placed: " . scalar(keys %hits) . "\n";
	print $rep "FLAG\t$_\n" for @flags;
	print $rep "NOTE\t$_\n" for @notes;
	print $rep "OK\tno strand switches\n" unless @flags || @notes;
	close $rep;

	$current_runtime = localtime();
	if(@flags){
		print $LOGFILE "$current_runtime\t******************WARNING: possible afin misjoin in $contigs.******************\n";
		print $LOGFILE "\t\t\t\t$_\n" for @flags;
		print $LOGFILE "\t\t\t\tDetails in $contigs.strand_check.txt.\n";
		print $SUMMARY "WARNING - possible afin misjoin (strand switch vs. reference) in $contigs: " . scalar(@flags) . " contig(s). See " . (File::Spec->splitdir(getcwd()))[-1] . "/$contigs.strand_check.txt\n";
	}
	else{
		print $LOGFILE "$current_runtime\tStrand check of $contigs against the reference: no strand switches.\n";
	}
	print $LOGFILE "\t\t\t\t$_\n" for @notes;
	return scalar @flags;
}

# Pull a single record (by first-token header id) out of a FASTA into $out.
sub extract_fasta_record {
        my ($fasta, $id, $out) = @_;
        open my $in, "<", $fasta or return 0;
        open my $o,  ">", $out   or return 0;
        my $printing = 0; my $found = 0;
        while(<$in>){
                if(/^>(\S+)/){
                        if($1 eq $id){ $printing = 1; $found = 1; }
                        elsif($printing){ last; }   # passed our record; stop
                        else{ $printing = 0; }
                }
                print $o $_ if $printing;
        }
        close $in; close $o;
        return $found;
}

##########

# Open a (possibly gzipped) reads file for streaming; returns a filehandle.
# Gzipped input is piped through pigz/gzip -dc (pigz also gets -p threads),
# which is fast and transparently handles multi-member gzip.
sub open_read_stream {
        my ($file) = @_;
        my $fh;
        if($file =~ /\.gz$/){
                my @cmd = ($PIGZ);
                push @cmd, "-p", $threads if $PIGZ =~ /pigz/;
                push @cmd, "-dc", $file;
                open($fh, "-|", @cmd) or die "ERROR: cannot decompress $file with $PIGZ: $!\n";
        }
        else{
                open($fh, "<", $file) or die "ERROR: cannot open $file: $!\n";
        }
        return $fh;
}

# Append a (possibly gzipped) reads file onto $dst, decompressing with
# pigz/gzip when needed. Used by the 'skip trim' path.
sub append_reads {
        my ($src, $dst) = @_;
        if($src =~ /\.gz$/){
                my $p = $PIGZ;
                $p .= " -p $threads" if $PIGZ =~ /pigz/;
                system("$p -dc '$src' >> '$dst'") == 0
                        or die "ERROR: failed to decompress $src into $dst\n";
        }
        else{
                system("cat '$src' >> '$dst'") == 0
                        or die "ERROR: failed to read $src into $dst\n";
        }
}


##########

# bowtie2 read arguments for whichever trimmed files exist under $trim_dir
# (paired, unpaired, or both). Returns undef when there are no reads at all.
# Replaces four hand-written if/elsif chains, one of which (the reassembly
# branch) assumed an unpaired file always existed for paired-end runs.
sub trimmed_read_args {
	my ($trim_dir) = @_;
	my $p1 = "$trim_dir/$name.trimmed_P1.fq";
	my $p2 = "$trim_dir/$name.trimmed_P2.fq";
	my $up = "$trim_dir/$name.trimmed_UP.fq";
	my $have_pe = (-s $p1 && -s $p2) ? 1 : 0;
	my $have_up = (-s $up) ? 1 : 0;
	return undef unless $have_pe || $have_up;
	my $args = " --al map_hits.fq";
	$args .= " --al-conc map_pair_hits.fq -1 $p1 -2 $p2" if $have_pe;
	$args .= " -U $up" if $have_up;
	return $args;
}

sub build_bowtie2_indices {

	my $bowtie_index = $_[0];
	my @terms = grep { length } map { my $t = $_; $t =~ s/^\s+|\s+$//g; $t } split /,/, $bowtie_index;
	$bowtie_index = join("_", @terms) if @terms > 1;

	# Whole-token, case-insensitive match against each record's taxonomy
	# (genus, species, tribe, subfamily, family, order). The old substring
	# match made "Poa" pull every Poaceae, and treated regex metacharacters
	# in the query as syntax.
	my $term_re = @terms ? '(?:^|\s)(?:' . join("|", map { quotemeta(lc $_) } @terms) . ')(?:\s|$)' : undef;
	my $want_all = (!defined $term_re || $bowtie_index =~ /^all$/i || $bowtie_index =~ /^genbank$/i) ? 1 : 0;

	# Taxonomy source. The current database uses bare-accession headers with a
	# sidecar metadata TSV (accession genus species tribe subfamily family order
	# source). Load it if present. If it's absent, fall back to the legacy
	# rich-header behavior so an older database (taxonomy encoded in the defline)
	# still works with this script -- and vice versa.
	my %order_of;    # accession -> order
	my %genus_of;    # accession -> genus (for one-per-genus reduction)
	my %tax_blob;    # accession -> lowercased taxonomy string, for --bowtie_index matching
	my $metafile = $FPBIN."/GenBank_Plastomes.metadata.tsv";
	if(open my $mfh, "<", $metafile){
		my $hdr = <$mfh>;   # column header
		while(<$mfh>){
			chomp;
			next unless length;
			my @c = split /\t/;
			my $acc = $c[0];
			next unless defined $acc && length $acc;
			$order_of{$acc} = (defined $c[6] ? $c[6] : "");
			$genus_of{$acc} = (defined $c[1] ? $c[1] : "");
			$tax_blob{$acc} = lc join(" ", grep { defined } @c[1..$#c]);
		}
		close $mfh;
	}
	my $have_meta = %order_of ? 1 : 0;
	unless($have_meta){
		print $LOGFILE "\t\t\t\tNo metadata table ($metafile) found; using legacy header-encoded taxonomy.\n";
	}

	my $gbfile = $FPBIN."/GenBank_Plastomes";
	my $outfsa = $bowtie_index.".fsa";

	# Per-record taxonomy for a header line: (accession, order, genus, tokens).
	my $taxonomy_of = sub {
		my ($hdr) = @_;
		my ($acc) = $hdr =~ /^>(\S+)/;
		$acc = "" unless defined $acc;
		if($have_meta && exists $order_of{$acc}){
			return ($acc, $order_of{$acc}, lc $genus_of{$acc}, $tax_blob{$acc});
		}
		# legacy rich header: Genus_species_tribe_subfamily_family_order_ACC_source
		(my $blob = $acc) =~ s/_/ /g;
		my ($order) = $hdr =~ /_(.*?ales)_/;
		my ($genus) = $hdr =~ /^>([^_]+)/;
		return ($acc, (defined $order ? $order : ""), (defined $genus ? lc $genus : ""), lc $blob);
	};

	# Stream the database, writing every record for which $keep->() is true.
	# Returns the number of records written.
	my $write_selected = sub {
		my ($keep) = @_;
		open my $in,  "<", $gbfile or die "Fast-Plast: cannot open $gbfile: $!\n";
		open my $out, ">", $outfsa or die "Fast-Plast: cannot write $outfsa: $!\n";
		my $printing = 0; my $n = 0;
		while(<$in>){
			if(/^>/){
				$printing = $keep->($_) ? 1 : 0;
				$n++ if $printing;
			}
			print $out $_ if $printing;
		}
		close $in; close $out;
		return $n;
	};

	my $written = 0;
	unless($want_all){
		# One sequence per unique genus among the matching taxa. Pulling every
		# plastome for a broad taxon (e.g. all Poales) builds a needlessly large
		# index; one representative per genus keeps the diversity that matters
		# for read reduction while staying small. Sequences whose genus is
		# unknown are kept (not collapsed), so no reference is silently dropped.
		my %used_genus;
		$written = $write_selected->(sub {
			my (undef, undef, $genus, $blob) = $taxonomy_of->($_[0]);
			return 0 unless $blob =~ /$term_re/;
			return 1 if $genus eq "" || $genus eq "na";
			return 0 if $used_genus{$genus}++;
			return 1;
		});
		if($written){
			print $LOGFILE "\t\t\t\t$written plastomes matching '$bowtie_index' (one per genus) used to make bowtie2 indices.\n";
		}
		else{
			print $LOGFILE "\t\t\t\tNo plastomes in the database match '$bowtie_index'. Using one representative from each order to make bowtie2 indices.\n";
			$want_all = 1;
		}
	}
	if($want_all){
		my %used_order;
		$written = $write_selected->(sub {
			my (undef, $order) = $taxonomy_of->($_[0]);
			return 0 if $order eq "";
			return 0 if $used_order{$order}++;
			return 1;
		});
		print $LOGFILE "\t\t\t\t$written plastomes (one per order) used to make bowtie2 indices.\n";
	}
	if($written == 0){
		print $LOGFILE "\t\t******************ERROR: no reference plastomes could be selected from $gbfile.\n";
		die "Fast-Plast: no reference plastomes could be selected from $gbfile (is the database intact?).\n";
	}

	my $build_bowtie2_exec = $BOWTIE2 . "-build " . $outfsa . " " . $name . "_bowtie";
	run_cmd($build_bowtie2_exec, "bowtie2-build");
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
my ($total_afin_contigs, $max_afin, $min_afin);

if($_[1]){
	($total_afin_contigs, $max_afin, $min_afin)  = &run_afin(50,100,15,2,$_[0],$extension, "extend");
}
else{
	($total_afin_contigs, $max_afin, $min_afin)  = &run_afin(50,100,15,2,$_[0],$extension);
}
print $LOGFILE "\t\t\t\tAfter afin, there are $total_afin_contigs contigs with a maximum size of $max_afin and a minimum size of $min_afin.\n";
$current_runtime = localtime();
print $LOGFILE "$current_runtime\tRemoving nested contigs.\n";
my $gotofinish;
if( $total_afin_contigs > 1){
	$current_afin = $name . "_afin_iter0.fa";

	run_cmd("$BLAST/makeblastdb -in $current_afin -dbtype nucl", "makeblastdb");
	my $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db " . $current_afin . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $current_afin . ".blastn";
	run_cmd($blast_afin_exec, "blastn");

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
				die "Fast-Plast: could not scaffold the contigs into a single plastome (coverage too low or uneven, or the reference was too distant). Best contigs are in Final_Assembly/. See the run log for details.\n";
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
				die "Fast-Plast: could not scaffold the contigs into a single plastome (coverage too low or uneven, or the reference was too distant). Best contigs are in Final_Assembly/. See the run log for details.\n";
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
	$current_afin = $name . "_afin_iter0.fa";
	
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
	my $afin_exec;
	if($_[6]){
		 $afin_exec = $AFIN_DIR . "/afin -c " . $_[4] . " -r ../2_Bowtie_Mapping/map_* -l " . $_[0] . " -f .1 -d " . $_[1] . " -x " . $extension . " -p " . $_[2] . " -i " . $_[3] . " -t " . $threads . " -o ". $name . "_afin --no_fusion";
	}
	else{	
		 $afin_exec = $AFIN_DIR . "/afin -c " . $_[4] . " -r ../2_Bowtie_Mapping/map_* -l " . $_[0] . " -f .1 -d " . $_[1] . " -x " . $extension . " -p " . $_[2] . " -i " . $_[3] . " -t " . $threads . " -o ". $name . "_afin";
	}
	print $LOGFILE "\t\t\t\tUsing command $afin_exec.\n";
	run_cmd($afin_exec, "afin");

	my %contig_lengths;
	my $afin_contig;
	my $max_afin=0;
	my $min_afin=100000000;
	my $iter_num=0;
	for my $input (@_){
		my $temp_count = ($input =~ tr/,//);
		if ($temp_count > $iter_num){
			$iter_num = $temp_count;
		}	
	}
	open my $afin_file, "<", $name . "_afin_iter".$iter_num.".fa";
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
	# afin can misjoin at repeats; place its output on a reference and flag
	# strand switches before anything downstream trusts it.
	check_strand_switches($name . "_afin_iter".$iter_num.".fa");
	return ($total_afin_contigs, $max_afin, $min_afin);
}

##########

sub cpgene_recovery {
	$current_runtime = localtime(); 
	
	print $LOGFILE "$current_runtime\tChecking chloroplast gene recovery in contigs.\n";
	my $current_afin = $_[0];
	my %chloroplast_db_genes;
	# Gene set: $posgenes (the bundled angiosperm set unless --posgenes was
	# given). Previously this sub hard-coded the bundled file, so a custom
	# gene set was never actually used for orientation.
	open my $cpdbgenes, "<", $posgenes or die "Fast-Plast: cannot open gene file $posgenes: $!\n";
	while(<$cpdbgenes>){
		chomp;
		if(/>/){
			/>(.*?)\_/;
			$chloroplast_db_genes{$1}=1;
		}
	}
	close $cpdbgenes;
	# Build the BLAST db once per working directory instead of on every call.
	unless(-e "Angiosperm_Chloroplast_Genes.nsq"){
		run_cmd("$BLAST/makeblastdb -in $posgenes -dbtype nucl -out Angiosperm_Chloroplast_Genes", "makeblastdb");
	}
	my $blast_afin_exec = $BLAST . "blastn -query " . $current_afin . " -db Angiosperm_Chloroplast_Genes -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $current_afin."_positional_genes" . ".blastn";
	run_cmd($blast_afin_exec, "blastn");
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

                if($max == $len1 && $len1 != $len2){
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
                                $blast_scores{$tarray[0]}{$tarray[1]}{$i}++;
                        }
                }
                elsif(exists $blast_scores{$tarray[1]}{$tarray[0]}){
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                $blast_scores{$tarray[1]}{$tarray[0]}{$i}++;
                        }
                }
                else{
                        for (my $j=0; $j<$min; $j++){
                                $blast_scores{$tarray[0]}{$tarray[1]}{$j}=0;
                        }
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                $blast_scores{$tarray[0]}{$tarray[1]}{$i}++;
                        }
                }
        }
        
        for my $blastid1 (keys %blast_scores){
                for my $blastid2 (keys %{$blast_scores{$blastid1}}){
                        my $overlap=0;
                        for my $i (keys %{$blast_scores{$blastid1}{$blastid2}}){
                                if($blast_scores{$blastid1}{$blastid2}{$i} > 0){
                                        $overlap ++;
                                }
                        }
                        

                        my $len_min1 = length($temp_seqs{$blastid1});
			my $len_min2 = length($temp_seqs{$blastid2});
			my $true_min=1000000000;
			my $remove_seq;
			if($len_min1 < $len_min2){
				$true_min = $len_min1;
				$remove_seq=$blastid1;
			}
			else{
				$true_min = $len_min2;
				$remove_seq=$blastid2;
			}
				
                        if ($overlap/$true_min >= 0.9){
                                $delete_contigs{">".$remove_seq}=0;
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
	for my $contig_name (keys %contigs_cp_genes){
		for my $contig_name2 (keys %contigs_cp_genes){
			if ($contig_name eq $contig_name2){
				next;
			}
			my $overlap=0;
			for my $gene_name (sort keys %{$contigs_cp_genes{$contig_name}}){
				if(exists $contigs_cp_genes{$contig_name2}{$gene_name}){
					$overlap++;
				}
			}
			my $genes1 = scalar keys %{$contigs_cp_genes{$contig_name}};
			my $genes2 = scalar keys %{$contigs_cp_genes{$contig_name2}};
                        if ($overlap/$genes1 >= 0.9 && $overlap/$genes2 <= 0.9){
				$delete_contigs{">".$contig_name}=0;
			}
			if ($overlap/$genes1 <= 0.9 && $overlap/$genes2 >= 0.9){
                                $delete_contigs{">".$contig_name2}=0;
                        }
			if ($overlap/$genes1 >= 0.9 && $overlap/$genes2 >= 0.9){
                                if(length($seqlens{$contig_name}) > length($seqlens{$contig_name2})){
					$delete_contigs{">".$contig_name2}=0;
				}
				else{
					$delete_contigs{">".$contig_name}=0;
				}
                        }
		}
	}

#	for my $gene_name (keys %genes_by_contig){
#		if(keys %{$genes_by_contig{$gene_name}}>1){
#			my $max = 0;
#			my $maxid;
			#for my $contig_name (keys %{$genes_by_contig{$gene_name}}){
			#	if($count_per_contig{$contig_name} > $max){
		#			$max = $count_per_contig{$contig_name};
			#		$maxid = $contig_name;
			#	}
			#}
#			for my $contig_name (keys %{$genes_by_contig{$gene_name}}){
#				if($contig_name ne $maxid){
#					$delete_contigs{">".$contig_name}=0;
#				}
#			}
#
#		}
#	}


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

        rename($current_seq."_fixed", $current_seq) or die "Fast-Plast: cannot replace $current_seq: $!\n";
}

##########
sub orientate_plastome{
        
        my $current_afin = $_[0];
        my $name = $_[1];
        my $path_to_final = $_[2];

        run_cmd("$BLAST/makeblastdb -in $posgenes -dbtype nucl -out Angiosperm_Chloroplast_Genes", "makeblastdb");

        for (my $i = 0; $i <=3; $i++){

        	run_cmd("perl $FPBIN/sequence_based_ir_id.pl $current_afin $name $i $min_region_length", "sequence_based_ir_id.pl");
        	my $split_fullname= $name ."_regions_split".$i.".fsa";

        	run_cmd("$BLAST/makeblastdb -in $split_fullname -dbtype nucl", "makeblastdb");
        	my $blast_afin_exec = $BLAST . "blastn -query " . $split_fullname . " -db " . $split_fullname . " -evalue 1e-40 -outfmt 6 -max_target_seqs 1000000 > " . $split_fullname . ".blastn";
			run_cmd($blast_afin_exec, "blastn");

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

			if($exists_ssc != 2){
				next;
			}
			
			if($exists_ir >1){
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
			unless($temp_range1 && $temp_range2 && exists $cp_piece_pos{"ir"} && ($temp_range1 == $cp_piece_pos{"ir"}{"len"} || $temp_range2 == $cp_piece_pos{"ir"}{"len"})){
				next;
			}
					
			# cpgene_recovery runs the positional-gene BLAST itself (and writes
			# the same <split>_positional_genes.blastn the orientation script reads).
			my ($percent_recovered_genes, $contigs_db_genes) = &cpgene_recovery($split_fullname);
			my %contigs_db_genes = %$contigs_db_genes;
			$percent_recovered_genes=$percent_recovered_genes*100;

			if($percent_recovered_genes > 75){
					run_cmd("perl $FPBIN/orientate_plastome_v.2.0.pl $split_fullname $split_fullname\_positional_genes.blastn $name", "orientate_plastome_v.2.0.pl");
        			my $final_seq = $name ."_FULLCP.fsa";
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
			chomp;
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
	if($_[0] =~ /CP_pieces/){
	open my $cppieces, "<", $_[0];
	while(<$cppieces>){
		chomp;
		if(/>/){
			$cpsid=substr($_,1);
		}
		else{
			$cplens{$cpsid}.=$_;
		}
	}
	}
	else{
		my $regions = identify_cp_regions($_[0]);
		unless($regions){
			print $LOGFILE "\t\t\t\tNo inverted repeat could be identified in $_[0]; per-region coverage averages skipped.\n";
			print $SUMMARY "Per-region coverage: not computed (no inverted repeat identified in the assembly)\n";
			return;
		}
		%cplens = %$regions;
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
###########
sub identify_cp_regions{
	my $temp_cpgenome;
	my $tcpid;
	open my $file, "<", $_[0]; #cp genome full
	while(<$file>){
        	chomp;
        	if(/>/){
                $tcpid = substr($_,1);
                next;
        	}
        else{
                $temp_cpgenome .= $_;
        }
	}	

	my $lsc;
	my $ssc;
	my $irb;
	my $ira;

	# The assembly is expected to end inside IRa, so the reverse complement of
	# its last 21 bases marks the start of IRb. If that 21-mer is not found the
	# assembly has no detectable IR (or does not end in one); report that
	# instead of scanning off the end of the sequence (which looped forever).
	return undef if length($temp_cpgenome) < 42;
	my $boundary1= substr(reverse($temp_cpgenome), 0, 21);
	$boundary1 =~ tr/ATCGatcg/TAGCtagc/;
	my $start_irb = index($temp_cpgenome, $boundary1);
	return undef if $start_irb < 0;

	$lsc = substr($temp_cpgenome, 0, $start_irb);

	my $i = $start_irb;

	my $irb_end = substr($temp_cpgenome, $i, 21);
	$irb_end = reverse($irb_end);
	$irb_end =~ tr/ATCGatcg/TAGCtagc/;

	# walk forward while the reverse complement of the current 21-mer still occurs
	until(length($irb_end) < 21 || index($temp_cpgenome, $irb_end) < 0){
        	$i++;
        	$irb_end = substr($temp_cpgenome, $i, 21);
        	$irb_end = reverse($irb_end);
        	$irb_end =~ tr/ATCGatcg/TAGCtagc/;
	}
	return undef if $i + 21 > length($temp_cpgenome);
	$i--;
	$irb_end = substr($temp_cpgenome, $i, 21);
	$irb = substr($temp_cpgenome, $start_irb, ($i+21-$start_irb));
	my $ira_seq = $irb_end;
	$ira_seq = reverse($ira_seq);
	$ira_seq =~ tr/ATCGatcg/TAGCtagc/;
	my $ira_start = index($temp_cpgenome, $ira_seq);
	$ssc = substr($temp_cpgenome, $i+21, $ira_start-($i+21));
	$ira= substr($temp_cpgenome, $ira_start);
        my %return_cp;
	$return_cp{"lsc"}=$lsc;
	$return_cp{"ssc"}=$ssc;
	$return_cp{"irb"}=$irb;
	return \%return_cp;
}	
###########

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
			$cass_seq .= $_;   # was '=', which kept only the last line of a wrapped FASTA
		}
	}
	# Keep the stretches between low-coverage regions. $c_start is the first
	# base after the previous bad region (0 before any). The old code tested
	# truthiness, so a bad region starting at position 0 was never skipped.
	my $c_start = 0;
	for my $starts (sort {$a <=> $b} keys %break_points){
		if($starts > $c_start){
			$new_substrings{$c_start}{$starts-1} = substr($cass_seq, $c_start, $starts-$c_start);
		}
		$c_start = $break_points{$starts}+1;
	}
	if($c_start < length($cass_seq)){
		$new_substrings{$c_start}{length($cass_seq)-1} = substr($cass_seq, $c_start);
	}

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
	--adapters		Files of adapters used in making sequencing library. Users can select "Nextera" for Nextera adapters, "TruSeq" for TruSeq adapters, leave the default (NEB), or provide their own. [Default = NEB-PE]
	--bowtie_index		Order for sample to draw references for mapping. If order exists, then all available samples for that order will be used. 
				If order does not exist in default set or the terms "all" or "GenBank" are given, one exemplar from each available order is used 
				to build the Bowtie2 indicies. [default="All"]
	--user_bowtie		User supplied bowtie2 indices. If this option is used, bowtie_index is ignored.
	--posgenes		User defined genes for identification of single copy/IR regions and orientation. Useful when major rearrangments are present in user plastomes.
	--coverage_analysis 	Flag to run the coverage analysis of a final chloroplast assembly.
	--min_region_length 	Minimum region length (passed on to sequence_based_ir_id.pl)
	--spades_only_assembler	Run SPAdes without read error correction (the behaviour before 1.3.1). Error correction is on by default.
	--min_length_trim	Minimum acceptable length for reads after trimming. [default = 140 for reads of 150 bp or longer, otherwise 90% of the read length]
	--posgenes		FASTA of genes used for LSC/SSC/IR identification and orientation. [default = bundled angiosperm gene set]
	--scaffold_reference	Reference plastome (FASTA) for RagTag scaffolding, overriding automatic selection.

=head1 DESCRIPTION

Fast-Plast is a pipeline that leverages existing and novel programs to quickly assemble, orient, and verify whole chloroplast genome sequences. For most datasets with sufficient data, Fast-Plast is able to produce a full-length de novo chloroplast genome assembly in approximately 30 minutes with no user mediation. 

Currently, Fast-Plast is written to accomodate Illumina data, though most data types could be used with a few changes.

Fast-Plast uses a de novo assembly approach by combining the de bruijn graph-based method of SPAdes with an iterative seed-based assembly implemented in afin to close gaps of contigs with low coverage. The pipeline then identifies regions from the quadripartite structure of the chloroplast genome, assigns identity, and orders them according to standard convention. A coverage analysis is then conducted to assess the quality of the final assembly. 

=head1 REQUIREMENTS

Fast-Plast requires fastp, bowtie2, SPAdes, and BLAST+.

If you use the coverage analysis to verify the assembly, then Jellyfish 2 and R will be needed. We highly recommend the coverage analysis to check the Fast-Plast assembly. 

afin requires a c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most *nix systems. See https://github.com/mrmckain/Fast-Plast for help with installation.

Fast-Plast is coded to use 4 threads during the fastp, bowtie2, SPAdes, and afin steps. This can simply be changed by the user if this number is not available.

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 50 million reads and running them through Fast-Plast.

=head1 INSTALLATION

All required programs should be in the user's path. 

To install afin:

		cd afin
		make

=head1 VERSION

Run C<fast-plast.pl --version> to print the installed version.





=cut

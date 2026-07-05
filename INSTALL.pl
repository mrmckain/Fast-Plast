#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use File::Spec;
use Getopt::Long;

#-----------------------------------------------------------------------------
# INSTALL.pl  --  set up Fast-Plast for use.
#
# As of v1.3.0 Fast-Plast resolves its external tools from your PATH at runtime
# (see the README), so installation no longer bakes absolute paths into the
# driver. This script therefore does three things:
#
#   1. compiles afin (the C++ seed-and-extend core),
#   2. checks that the required dependencies are visible on your PATH, and
#   3. downloads the reference plastome database from Zenodo.
#
# The easiest way to satisfy the dependencies is the provided conda
# environment:  conda env create -f environment.yml && conda activate fast-plast
#
# Usage:
#   perl INSTALL.pl [--skip-afin] [--skip-deps] [--skip-db] [--help]
#-----------------------------------------------------------------------------

my ($skip_afin, $skip_deps, $skip_db, $help);
GetOptions(
    'skip-afin' => \$skip_afin,
    'skip-deps' => \$skip_deps,
    'skip-db'   => \$skip_db,
    'help|?'    => \$help,
) or die "Invalid option. Try: perl INSTALL.pl --help\n";

if ($help) {
    print <<"USAGE";
Fast-Plast installer

  perl INSTALL.pl [options]

  --skip-afin   do not compile afin
  --skip-deps   do not check for dependencies on PATH
  --skip-db     do not download the reference plastome database
  --help        show this message

Dependencies are easiest to obtain via conda:
  conda env create -f environment.yml
  conda activate fast-plast
USAGE
    exit 0;
}

my $ROOT = $FindBin::RealBin;
my @problems;   # collected required-dependency problems -> nonzero exit

# Locate an executable on PATH (honors an optional FP_* override, mirroring the
# resolution the driver uses, so a check here matches what a run will see).
sub which_exe {
    my ($name, $envvar) = @_;
    if ($envvar && defined $ENV{$envvar} && length $ENV{$envvar}) {
        my $o = $ENV{$envvar};
        return $o if -x $o && ! -d $o;
        my $c = File::Spec->catfile($o, $name);
        return $c if -x $c && ! -d $c;
    }
    for my $dir (File::Spec->path) {
        my $c = File::Spec->catfile($dir, $name);
        return $c if -x $c && ! -d $c;
    }
    return undef;
}

print "Fast-Plast installation\n";
print "  repository: $ROOT\n\n";

#-----------------------------------------------------------------------------
# 1. Compile afin
#-----------------------------------------------------------------------------
unless ($skip_afin) {
    print "== Compiling afin ==\n";
    my $afin_dir = File::Spec->catdir($ROOT, "afin");
    my $cxx = which_exe('g++') || which_exe('c++') || which_exe('clang++');
    if (!$cxx) {
        print "  ! No C++ compiler (g++/c++/clang++) found on PATH.\n";
        print "    Install one, or 'conda install compilers', then re-run.\n";
        push @problems, "afin: no C++ compiler";
    }
    elsif (! -d $afin_dir) {
        print "  ! afin/ directory not found at $afin_dir\n";
        push @problems, "afin: source directory missing";
    }
    else {
        my $rc = system("cd \"$afin_dir\" && make");
        my $afin_bin = File::Spec->catfile($afin_dir, "afin");
        if ($rc == 0 && -x $afin_bin) {
            print "  afin compiled: $afin_bin\n";
        }
        else {
            print "  ! afin failed to compile. See the make output above.\n";
            push @problems, "afin: compilation failed";
        }
    }
    print "\n";
}

#-----------------------------------------------------------------------------
# 2. Dependency check
#-----------------------------------------------------------------------------
unless ($skip_deps) {
    print "== Checking dependencies on PATH ==\n";

    # name => [ FP_* override var, required?, note ]
    my @required = (
        ['fastp',       'FP_FASTP',     'read trimming'],
        ['bowtie2',     'FP_BOWTIE2',   'read mapping'],
        ['bowtie2-build','FP_BOWTIE2',  'bowtie2 index build'],
        ['spades.py',   'FP_SPADES',    'initial assembly'],
        ['blastn',      'FP_BLAST',     'gene/reference checks'],
        ['makeblastdb', 'FP_BLAST',     'BLAST database build'],
        ['ragtag.py',   undef,          'scaffolding'],
        ['minimap2',    undef,          'required by RagTag'],
    );
    my @optional = (
        ['pigz',        'FP_PIGZ',      'fast decompression (gzip used if absent)'],
        ['jellyfish',   'FP_JELLYFISH', 'coverage analysis'],
        ['Rscript',     undef,          'coverage plotting'],
    );

    for my $t (@required) {
        my ($name, $env, $note) = @$t;
        my $p = which_exe($name, $env);
        if ($p) { printf "  [ok]   %-14s %s\n", $name, $p; }
        else {
            printf "  [MISS] %-14s (required: %s)\n", $name, $note;
            push @problems, "missing dependency: $name";
        }
    }
    for my $t (@optional) {
        my ($name, $env, $note) = @$t;
        my $p = which_exe($name, $env);
        if ($p) { printf "  [ok]   %-14s %s\n", $name, $p; }
        else    { printf "  [warn] %-14s (optional: %s)\n", $name, $note; }
    }
    # pigz is optional only if gzip exists
    unless (which_exe('pigz','FP_PIGZ') || which_exe('gzip','FP_GZIP')) {
        print "  [MISS] neither pigz nor gzip found (needed to read .gz inputs)\n";
        push @problems, "missing dependency: pigz/gzip";
    }
    print "\n";
}

#-----------------------------------------------------------------------------
# 3. Reference database
#-----------------------------------------------------------------------------
unless ($skip_db) {
    print "== Reference plastome database ==\n";
    my $fetch = File::Spec->catfile($ROOT, "bin", "fetch_plastome_db.sh");
    my $fasta = File::Spec->catfile($ROOT, "bin", "GenBank_Plastomes");
    my $meta  = File::Spec->catfile($ROOT, "bin", "GenBank_Plastomes.metadata.tsv");

    if (-s $fasta && -s $meta) {
        print "  database already present in bin/; skipping download.\n";
    }
    elsif (! -e $fetch) {
        print "  ! fetch script not found: $fetch\n";
        push @problems, "database: fetch script missing";
    }
    else {
        chmod 0755, $fetch;
        my $rc = system("bash \"$fetch\"");
        if ($rc == 0 && -s $fasta && -s $meta) {
            print "  database installed into bin/.\n";
        }
        else {
            print "  ! database download did not complete. You can retry later with:\n";
            print "      bash bin/fetch_plastome_db.sh\n";
            push @problems, "database: download failed";
        }
    }
    print "\n";
}

#-----------------------------------------------------------------------------
# make bundled bin/ scripts executable
#-----------------------------------------------------------------------------
for my $s (glob(File::Spec->catfile($ROOT, "bin", "*.pl")),
           glob(File::Spec->catfile($ROOT, "bin", "*.sh")),
           glob(File::Spec->catfile($ROOT, "bin", "*.py"))) {
    chmod 0755, $s;
}
chmod 0755, File::Spec->catfile($ROOT, "fast-plast.pl");

#-----------------------------------------------------------------------------
# summary
#-----------------------------------------------------------------------------
print "== Summary ==\n";
if (@problems) {
    print "  The following need attention:\n";
    print "    - $_\n" for @problems;
    print "\n  Resolve these (conda is the easy path) and re-run INSTALL.pl.\n";
    exit 1;
}
else {
    print "  Fast-Plast is ready.\n";
    print "  Try:  perl fast-plast.pl --help\n";
    exit 0;
}

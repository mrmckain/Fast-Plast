#!/usr/bin/env perl
use strict;
use warnings;

# normalize_plastome_db.pl
#
# One-time build step for the Fast-Plast reference plastome database.
#
# Converts the legacy "rich header" GenBank plastome FASTA, e.g.
#   >Arabidopsis_thaliana_Camelineae_NA_Brassicaceae_Brassicales_NC_000932.1_Genbank
# into:
#   (1) a clean FASTA whose headers are bare, unique accessions (>NC_000932.1),
#       which are short, BLAST-safe (no makeblastdb id-length abort, and
#       -parse_seqids works again if ever wanted), and
#   (2) a sidecar metadata TSV mapping accession -> taxonomy, which Fast-Plast
#       loads for per-order sampling and --bowtie_index taxon selection.
#
# The legacy header layout is positional:
#   genus _ species _ tribe _ subfamily _ family _ order _ <accession> _ source
# where <accession> may itself contain an underscore (RefSeq "NC_000932.1").
# Accession detection is therefore done by regex (robust to that underscore),
# while the six taxonomy ranks are read positionally. Anything that fails to
# parse cleanly is still emitted (never silently dropped) under a generated id
# and written to a .rejected report for manual review.
#
# Usage:
#   normalize_plastome_db.pl <input_rich.fsa> <output_clean.fsa> <output_meta.tsv>
#
# Streams input (safe for multi-GB files); sequence lines are passed through
# verbatim under the rewritten header.

@ARGV == 3 or die "Usage: $0 <input_rich.fsa> <output_clean.fsa> <output_meta.tsv>\n";
my ($infile, $outfasta, $outmeta) = @ARGV;

open my $IN,   "<", $infile   or die "ERROR: cannot open input '$infile': $!\n";
open my $FA,   ">", $outfasta or die "ERROR: cannot write '$outfasta': $!\n";
open my $META, ">", $outmeta  or die "ERROR: cannot write '$outmeta': $!\n";
open my $REJ,  ">", $outfasta.".rejected" or die "ERROR: cannot write reject log: $!\n";

print $META join("\t", qw(accession genus species tribe subfamily family order source)), "\n";

# Accession patterns: RefSeq (NC_000932.1) or GenBank (MN123456.1 / X12345.1).
# No \b anchors: underscore is a word char, so \b never fires at "_NC". The
# accession is the only digit-bearing token in these headers, so an unanchored
# match is safe; RefSeq is tried first so it wins over the GenBank alternative.
my $ACC_RE = qr/((?:[A-Z]{2}_\d+\.\d+)|(?:[A-Z]{1,2}\d{5,8}\.\d+))/;

my %seen_acc;          # accession -> count, to enforce uniqueness
my ($n_total, $n_clean, $n_reject, $n_dupe) = (0,0,0,0);
my $fallback = 0;

while (my $line = <$IN>) {
    if ($line =~ /^>/) {
        $n_total++;
        chomp(my $header = $line);
        $header =~ s/^>//;

        # --- accession (regex; tolerant of the underscore in RefSeq ids) ---
        my $acc;
        if ($header =~ $ACC_RE) { $acc = $1; }

        # --- taxonomy (positional split on underscore) ---
        my @f = split /_/, $header;
        my ($genus,$species,$tribe,$subfamily,$family,$order,$source) =
            ("NA","NA","NA","NA","NA","NA","NA");
        my $clean_parse = 0;
        if (@f >= 8) {
            ($genus,$species,$tribe,$subfamily,$family,$order) = @f[0..5];
            $source = $f[-1];
            # sanity: botanical family ends -aceae, order ends -ales
            $clean_parse = ($family =~ /aceae$/i && $order =~ /ales$/i) ? 1 : 0;
        }

        # --- resolve a unique, short id ---
        if (!defined $acc) {
            $fallback++;
            $acc = sprintf("FP_unplaced_%05d", $fallback);
            print $REJ "no-accession\t$header\t-> $acc\n";
            $n_reject++;
        }
        elsif (!$clean_parse) {
            # accession is fine but taxonomy fields look off; keep the record,
            # flag it so the ranks can be corrected upstream if needed.
            print $REJ "taxonomy-parse\t$header\t(family='$family' order='$order')\n";
            $n_reject++;
        }

        if ($seen_acc{$acc}++) {
            my $dup = $acc . "_dup" . $seen_acc{$acc};
            print $REJ "duplicate-accession\t$header\t$acc -> $dup\n";
            $acc = $dup;
            $n_dupe++;
        }

        print $FA ">$acc\n";
        print $META join("\t", $acc,$genus,$species,$tribe,$subfamily,$family,$order,$source), "\n";
        $n_clean++;
    }
    else {
        print $FA $line;   # sequence line, verbatim
    }
}

close $IN; close $FA; close $META; close $REJ;

print STDERR "normalize_plastome_db: $n_total records\n";
print STDERR "  clean FASTA + metadata written: $n_clean\n";
print STDERR "  flagged for review (see $outfasta.rejected): $n_reject\n";
print STDERR "  duplicate accessions disambiguated: $n_dupe\n";
print STDERR "  records needing a generated id: $fallback\n";

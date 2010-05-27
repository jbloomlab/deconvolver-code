#!/usr/local/bin/perl -w
use strict;

# Test: Make sure reads are never assigned to a barcode that it does not match
my $outdir = shift @ARGV || die "A output directory is required to validate the results of grid-deconvolve.";

# Read in fuzznuc results and build a table of sequences-barcodes that is unique
my %seqs;

opendir (OUTDIR, $outdir) || die "Could not open directory $outdir\n";
my @entries = readdir(OUTDIR);
foreach my $fuzznuc_file (@entries) {
	if ($fuzznuc_file =~ /^fuzznuc\./) {
		open (FUZZNUC, "< $outdir/$fuzznuc_file") || die "Could not open file $fuzznuc_file\n";
		#print "Fuzznuc file $outdir/$fuzznuc_file\n";
		while (<FUZZNUC>) {
			chomp;
			next if $_=~/SeqName/;
			my ($seq_id, $start, $end, $length, $strand, $barcode, $num_mismatches) = split /\t/, $_;
			$seqs{$barcode}{$seq_id}++;
		}
		close FUZZNUC;
	}
}
die "No fuzznuc search results were found in $outdir\n" unless scalar keys %seqs;

# Read in results and validate the assignments
my $num_invalid = 0; my $num_valid = 0;
foreach my $barcode_dir (@entries) {
	if (-d "$outdir/$barcode_dir") {
		next if $barcode_dir eq "." || $barcode_dir eq "..";
		opendir (BARCODE_DIR, "$outdir/$barcode_dir") || die "Could not open barcode directory $outdir/$barcode_dir\n";
		my @trim_files = readdir(BARCODE_DIR);
		foreach my $trim_file (@trim_files) {
			if ($trim_file =~ /\.trim/) {
				open (ASSIGNMENTS, "< $outdir/$barcode_dir/$trim_file" ) || die "Could not open trim file $outdir/$barcode_dir/$trim_file\n";
				
				# get barcode from the name of the file
				my $barcode = $trim_file;
				$barcode =~ s/\.trim//g;
				
				while (<ASSIGNMENTS>) {
					chomp;
					my ($seq_id) = split /\s+/, $_;
					if ($seq_id) {
						if (!exists $seqs{$barcode}{$seq_id}) {
							$num_invalid++;
							print "$seq_id $barcode\n";
						} else {
							$num_valid++;
						}
					}
				}
				close ASSIGNMENTS;
			}
		}
		close BARCODE_DIR;
	}
}

my $num_barcodes = scalar keys %seqs;
my $num_seqs = 0;
foreach my $barcode( keys %seqs ) {
	my $num_barcode_seqs = scalar keys %{ $seqs{$barcode} };
	$num_seqs += $num_barcode_seqs;
	print "# barcode $barcode $num_barcode_seqs\n";
}
print "Number of barcodes with assigned reads: $num_barcodes\n";
print "Total sequences: $num_seqs\n";
print "Total valid assignments: $num_valid\n";
print "Total invalid assignments: $num_invalid\n";

close OUTDIR;

#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Barcode::Deconvolver;
use Test::More;

# Get the user
my $user=`whoami`; chomp $user;
$user ||= "deconvolve";

# Set inputs files and output directory
my ($fasta, $pattern, $answers_file, $dir, $key) = @ARGV;
$fasta ||= "$Bin/data/test_GEIVOUP02.sff";
$pattern ||= "$Bin/data/test_barcode_metadata_from_GLK.txt.pat";
$answers_file = "$Bin/data/test_GEIVOUP02_expected_answer_trim_BC019CG.txt.trim";
if (!$dir) {
	$dir = "/usr/local/scratch/$user";
	mkdir $dir unless -e $dir;
	$dir .= "/test_deconvolve";
}
mkdir $dir unless -e $dir;
die "Output directory $dir is not writeable.\n" unless -w $dir;
die "Unable to read barcode pattern file: $pattern.\n" unless -r $pattern;
die "Unable to read fasta file: $fasta.\n" unless -r $fasta;
die "Unable to read answers file: $answers_file.\n" unless -r $answers_file;

# Set the input, tmp and output directories to $dir
my ($outdir, $errdir, $tmpdir);
$outdir = $errdir = $tmpdir = $dir;

# Issue the test plan when we call done_testing($num_tests)
my $num_tests = 0;

# Get our Grid::Tools::Barcode::Deconvolver object
# Test (1a): Using options stringified
my $Deconvolver = new Grid::Tools::Barcode::Deconvolver({
			key		=> $key,
			pattern	=> $pattern,
			infile	=> $fasta,
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			cleanup 	=> 0,
			verbose	=> 1
});

# Run the grid deconvolution pipeline
$Deconvolver->run();

# Test fuzznuc grid searches
my $Fuzznuc = $Deconvolver->fuzznuc;

# Get the trim results made by our deconvolution pipeline
my $trim_table = $Deconvolver->trim_table;

# Test the results match our answer file
my @answers = read_answers_file($answers_file);
foreach my $answer (@answers) {
	my ($seq_id, $barcode, $clear_start, $clear_end) = @$answer;
	my ($assigned_barcode, $assigned_start, $assigned_end) = get_trim_points($trim_table, $seq_id);
	ok($assigned_barcode eq $barcode && $assigned_start == $clear_start && $assigned_end == $clear_end, 
			"$seq_id answer $barcode $clear_start..$clear_end assigned $assigned_barcode $assigned_start..$assigned_end");
	$num_tests++;
}

done_testing($num_tests);


######################## SUB ROUTINES ############################

sub read_answers_file {
	my $answers_file = shift;
	my @answers;
	open (ANSWERS, "< $answers_file" ) || die "Could not open answers file $answers_file\n";
	while (<ANSWERS>) {
		chomp;
		my ($barcode, $seq_id, $clear_start, $clear_end) = split /\s+/, $_;
		push @answers, [ $seq_id, $barcode, $clear_start, $clear_end ];
	}
	return @answers;
}

sub get_trim_points {
	my ($trim_table, $seq_id) = @_;
	my $trim_results = $trim_table->{$seq_id} || print STDERR "Error: No trim results found for sequence $seq_id\n";
	if ($trim_results && ref($trim_results) eq "ARRAY") {
		return @$trim_results;
	} else {
		return ("", "", ""); # [ $clear_start, $clear_end, $reason ]
	}
}



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

# test using relative paths
# ./grid_deconvolve.t data/test_GEIVOUP02.sff data/test_barcode_metadata_from_GLK.txt.pat data/test_GEIVOUP02_expected_answer_trim_BC019CG.txt.trim 

# Set inputs files and output directory
my ($fasta, $pattern, $answers_file, $dir, $key) = @ARGV;
$fasta ||= "$Bin/data/test_GEIVOUP02.sff";
$pattern ||= "$Bin/data/test_barcode_metadata_from_GLK.txt.pat";
$answers_file = "$Bin/data/test_GEIVOUP02_expected_answer_trim_BC019CG.txt.trim";
if (!$dir) {
	$dir = "/usr/local/scratch/$user";
	mkdir $dir unless -e $dir;
	$dir .= "/test_grid_deconvolve";
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

# Get our SGE object based on options
my $grid = new Grid::SGE({
			project 	=> 810001,
			queue	=> "fast.q",
			name		=> "gridDeconvolve",
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			verbose	=> 1,
});

# Get our Grid::Tools::Barcode::Deconvolver object
# Test (1a): Using options stringified
my $Deconvolver = new Grid::Tools::Barcode::Deconvolver({
			key		=> $key,
			grid		=> $grid,
			pattern	=> $pattern,
			infile	=> $fasta,
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			clamplength => 6, 		# Test that barcode-specific clamplength in pattern file is used
			cleanup 	=> 0,
			verbose	=> 1
});

# Run the grid deconvolution pipeline
$Deconvolver->run();

# Test fuzznuc grid searches
my $Fuzznuc = $Deconvolver->fuzznuc;

# Test that all tasks in fuzznuc job array completed successfully
my $num_tasks = $grid->num_tasks;
my $Tasks = $grid->job_tasks;
my ($num_success, $num_failed);
foreach my $T (@$Tasks) {
	if ($T->status eq "complete") {
		$num_success++;
	} else {
		$num_failed++;
	}
}
ok($num_tasks == $num_success, "All fuzznuc tasks are successfully completed");
ok(!$num_failed, "No job failures");
$num_tests += 2;

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

# Report any job failures
print $grid->tasks_report() if $num_failed;

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



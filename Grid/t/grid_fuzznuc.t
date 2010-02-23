#!/usr/local/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Fuzznuc;
use Test::More tests => 26;
use File::Basename;

my $verbose = 1;
my $outdir = "/usr/local/scratch/naxelrod/out";
my $errdir = "/usr/local/scratch/naxelrod/err";
my $tmpdir = "/usr/local/scratch/naxelrod/tmp";
my $fasta = "/usr/local/devel/ANNOTATION/naxelrod/lib/Grid/t/data/test_seqs.fasta";

# Use a pattern file
my $pattern1 = "$Bin/../t/data/barcodes.pat"; 	# /usr/local/projects/VHTNGS/sample_data/HI_6_100/barcodes/barcode_metadata_from_GLK.txt
my $pattern2 = "[CG](5)TG{A}N(1,5)C";

# Get our SGE object based on options
my $grid = new Grid::SGE({
			project 	=> 810001,
			name		=> "gridFuzznuc",
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			verbose	=> $verbose,
			poll_delay => 10,
});

################################################################################
# Test (1a): Using options stringified
################################################################################
my $Fuzznuc = new Grid::Tools::Fuzznuc({
			fasta	=> $fasta,
			verbose	=> $verbose,
			options	=> "-pattern \@$pattern1 -pmismatch 2 -filter -rformat excel -stdout true -complement Yes",
});

$grid->add_command($Fuzznuc->command($grid));

my @Jobs = $grid->submit_and_wait();
foreach my $Job (@Jobs) {
	ok($Job->status eq "complete", "positive test (1a)");
	print $Job->to_string() if $grid->verbose;
}

# Get the output files
my $files = $Fuzznuc->get_outfiles($grid);

# Get an iterator of FileIO::FuzznucHit objects from our output files
my $iter = FileIO::FuzznucUtils->hit_iterator($files);

# Test that we have a Hit result with a sequence id
my $Hit = $iter->();
ok($Hit->seq_id, "test hit iterator");

################################################################################
# Test (1b): using options hash
################################################################################\
# Test: Use an invalid command option that causes fuzznuc to die
# Get our Fuzznuc object
$Fuzznuc = new Grid::Tools::Fuzznuc({
			fasta	=> [ $fasta ],
			options 	=> {
				pattern 		=> $pattern2,
				pmismatch		=> 2,
				filter		=> undef,
				rformat		=> "excel",
				stdout		=> "true",
				complement	=> "Yes"
				
			},
			outdir	=> $outdir,
			tmpdir	=> $tmpdir,
});

$grid->add_command($Fuzznuc->command($grid));

@Jobs = $grid->submit_and_wait();
foreach my $Job (@Jobs) {
	ok($Job->status eq "complete", "positive test (1b)");
	print $Job->to_string() if $grid->verbose;
}

# Test that we have a Hit result with a sequence id
$files = $Fuzznuc->get_outfiles($grid);
$iter = FileIO::FuzznucUtils->hit_iterator($files);
$Hit = $iter->();
ok($Hit->seq_id, "test hit iterator");


################################################################################
# Test (2a): Expected failure due to invalid fuzznuc parameter
################################################################################\
$grid->name("TestBadOpt");
$Fuzznuc = new Grid::Tools::Fuzznuc({
			fasta	=> $fasta,
			options	=> "-nonexistent_option blah -pattern \"[CG](5)TG{A}N(1,5)C\" -pmismatch 2",
});

$grid->add_command($Fuzznuc->command($grid));

@Jobs = $grid->submit_and_wait();
foreach my $Job (@Jobs) {
	ok($Job->status eq "failed", "negative test (2a)");
	print $Job->to_string() if $grid->verbose;
}


################################################################################
# Test (2b): Expected failure due to invalid output directories to Grid object
################################################################################\
my $grid_baddir = new Grid::SGE({
			project 	=> 810001,
			name		=> "TestBadDir",
			verbose	=> $verbose,
			poll_delay => 10
});

$Fuzznuc = new Grid::Tools::Fuzznuc({
			tmpdir	=> "/path/to/nowhere/",
			fasta	=> $fasta,
			options	=> "-pattern \"[CG](5)TG{A}N(1,5)C\" -pmismatch 2",
});

$grid_baddir->add_command($Fuzznuc->command($grid_baddir));
@Jobs = $grid_baddir->submit_and_wait();
foreach my $Job (@Jobs) {
	ok($Job->status eq "failed", "negative test (2a)");
	print $Job->to_string() if $grid->verbose;
}


# not sure how many tests will be run
done_testing();


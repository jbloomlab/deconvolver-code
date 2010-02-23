#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Barcode::Deconvolver;
use Test::More tests => 24;

my $outdir = "/usr/local/scratch/naxelrod/out";
my $errdir = "/usr/local/scratch/naxelrod/err";
my $tmpdir = "/usr/local/scratch/naxelrod/tmp";
my $fasta = "$Bin/../t/s_1_1_sequence.fastq"; 	# /usr/local/devel/ANNOTATION/naxelrod/lib/Grid/t/data/test_seqs.fasta
my $pattern = "$Bin/../t/barcodes.pat"; 	# /usr/local/projects/VHTNGS/sample_data/HI_6_100/barcodes/barcode_metadata_from_GLK.txt

# Get our SGE object based on options
my $grid = new Grid::SGE({
			project 	=> 810001,
			name		=> "gridDeconvolve",
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			verbose	=> 1,
			poll_delay => 30,
});

# Get our Grid::Tools::Barcode::Deconvolver object
# Test (1a): Using options stringified
my $Deconvolver = new Grid::Tools::Barcode::Deconvolver({
			grid		=> $grid,
			pattern	=> $pattern,
			infile	=> $fasta,
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			verbose	=> 1,
			options 	=> "-pmismatch 2 -filter -rformat excel -stdout true -complement Yes"
});

# Validate our input parameters
my $error_message = $Deconvolver->validate_parameters();
die $error_message if $error_message;

# Store the fasta sequences in a hash table
$Deconvolver->fasta_file("/usr/local/scratch/naxelrod/tmp/s_1_1_sequence.fasta");
$Deconvolver->quals_file("/usr/local/scratch/naxelrod/tmp/s_1_1_sequence.quals");
my $fasta_table = $Deconvolver->make_fasta_table;
my $id = "SOLEXA1_1_36_1096_1964#0/1";
my $F = $fasta_table->{$id};
print join(" ", $F->id, $F->seq), "\n";

# Log runtime settings
print $Deconvolver->runtime_settings_string();

# Generate Fasta sequence and quality files (required for fuzznuc searches)
$Deconvolver->write_temp_fasta_files;

# Run fuzznuc searches on the grid
# Note: if no grid object is provided, it runs the searches in serial
my $is_success = $Deconvolver->run_fuzznuc($grid);

# Test that running fuzznuc on the grid was successful
ok($is_success, "grid fuzznuc searches");

# Report a list of failed tasks if searches fail
die $grid->failed_tasks_report() if !$is_success;

# Get an iterator of FileIO::FuzznucHit objects from our output files
my $iter = $Deconvolver->fuzznuc->get_hits_iterator($grid);

# Test that we have a function reference
ok(ref($iter) eq "CODE", "test hit results iterator");

# Read in our barcode file
my $barcode_table = $Deconvolver->read_pattern_file;

# Store the fasta sequences in a hash table
# my $fasta_table = $Deconvolver->make_fasta_table;

# Assign sequences to barcodes, and make our multicode table
my $assignment_table = $Deconvolver->make_assignment_table($iter);

# Log our multicoded sequences
$Deconvolver->write_multicode_report();

# Writes barcode fasta files
$Deconvolver->write_barcode_fasta;

done_testing();



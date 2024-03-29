#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Barcode::Deconvolver;
use Test::More tests => 8;

# Get the user
my $user=`whoami`; chomp $user;
$user ||= "deconvolve";

my ($fasta, $pattern, $dir) = @ARGV;
$fasta ||= "/usr/local/devel/VIRIFX/users/naxelrod/data/FTF2AAH01.sff";
$pattern ||= "/usr/local/devel/VIRIFX/users/naxelrod/data/barcodes.pat";
$dir ||= "/usr/local/scratch/$user";
die "Output directory $dir is not writeable.\n" unless -w $dir;
die "A barcode pattern file is required.\n" unless -r $pattern;
die "An input fasta file is required.\n" unless -r $fasta;

# Set the input files and output directories
my $outdir = "$dir/out";
my $errdir = "$dir/err";
my $tmpdir = "$dir/tmp";

# Get our SGE object based on options
my $grid = new Grid::SGE({
			project 	=> 810001,
			queue	=> "fast.q",
			name		=> "gridDeconvolve",
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			verbose	=> 1,
			poll_delay => 60,
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

# Test that barcode pattern file ($pattern) was read correctly
my $default_readlength = $Deconvolver->readlength();
ok($default_readlength == 50, "Default minimum read length is 50.");
ok($Deconvolver->readlength('BC009CG') == 70, "Minimum read length for barcode BC009CG is 70.");

# Let's set some barcode-specific values (by hash)
$Deconvolver->clamplength({ 
	'BC009CG' => 14, 
	'BC015CG' => 10 
});
ok($Deconvolver->clamplength('BC009CG') == 14, "Test setting clamplength for barcode BC009CG by hash.");

# Or, set as key-value pairs
$Deconvolver->clamplength('BC009CG', 10);
ok($Deconvolver->clamplength('BC009CG') == 10, "Test setting clamplength for barcode BC009CG by argument.");

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

# Test that a "reasonable" number of assignments were made
my $num_seqs = scalar keys %{ $Deconvolver->fasta_table };
my $num_seqs_with_hits = $Deconvolver->num_seqs_with_hits;
my $num_assignments = $Deconvolver->num_assignments;
my $perc_seqs_with_hits = ($num_seqs) ? int(($num_seqs_with_hits/$num_seqs)*100) : 0;
my $perc_assignments = ($num_seqs) ? int(($num_assignments/$num_seqs)*100) : 0;
print join(", ", $num_seqs, $num_seqs_with_hits, $num_assignments, $perc_seqs_with_hits, $perc_assignments), "\n";
cmp_ok($perc_seqs_with_hits, '>=', 90, "Test at least 90\% of sequences have barcode hits ($perc_seqs_with_hits)");
cmp_ok($perc_assignments, '>=', 50, "Test at least 50\% of sequences are assigned to barcodes ($perc_assignments)");

# 601967, 566509, 95, 94, 0
# perc_assign = 95

# Report any job failures
print $grid->tasks_report() if $num_failed;

done_testing();



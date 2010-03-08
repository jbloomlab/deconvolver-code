#!/usr/local/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Blast;
use Test::More;

# blastall -p blastp -d /usr/local/scratch/naxelrod/panda/panda.1433906054260851320.fasta -i /usr/local/scratch/naxelrod/25K-vics-2009-11-10-v23021-r756-FPBMLOC01-560k-without-clear-range-metagene/metagene_mapped_pep.p1.fasta -o some_output_file -v 10 -b 10 -X 15 -e 1e-5 -K 10 -f 11 -Z 25.0 -W 3 -U F -y 7.0 -A 40 -m 7

my $verbose = 0;
my $basedir = "/usr/local/scratch/naxelrod";
my $outdir = "$basedir/out";
my $errdir = "$basedir/err";
my $tmpdir = "$basedir/tmp";
my $fasta = "/usr/local/devel/ANNOTATION/naxelrod/lib/Grid/t/data/test_seqs.fasta";
# my $fasta = "/usr/local/scratch/kkrampis/25K-vics-2009-11-10-v23021-r756-FPBMLOC01-560k-without-clear-range-metagene/metagene_mapped_pep.p1.fasta";
my $database = "/usr/local/scratch/kkrampis/panda/panda.1433906054260851320.fasta";
my $program = "blastp";
my $options = "-v 10 -b 10 -X 15 -e 1e-5 -K 10 -f 11 -Z 25.0 -W 3 -U F -y 7.0 -A 40 -m 7";
my $num_tests = 0;

# TODO: Need to check for failed cases (such as invalid fasta files, etc)
# BUGS
# Need to add a unique key to the temporary files that are created to avoid,
# simultaneous runs from overwriting each other.

# Get our SGE object based on options
my $grid = new Grid::SGE({
			project 	=> 810001,
			name		=> "gridBlastJob",
			tmpdir	=> $tmpdir,
			errdir	=> $errdir,
			outdir	=> $outdir,
			verbose	=> $verbose,
			poll_delay => 300, # every 5 minutes
});

################################################################################
# Test (1a): Using options stringified
################################################################################\
my $Blast = new Grid::Tools::Blast({
			fasta	=> $fasta,
			num_seqs	=> 500,
			options	=> "-p $program -d $database $options",
});

$grid->add_command($Blast->command($grid));

my @Jobs = $grid->submit_and_wait();
foreach my $Job (@Jobs) {
	print $Job->to_string;
	ok($Job->status eq "complete", "positive test (1a)");
}
$num_tests += scalar @Jobs;

done_testing($num_tests);

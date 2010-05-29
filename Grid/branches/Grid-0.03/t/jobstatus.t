#!/usr/local/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Blast;
use Test::More;

# Get our SGE object based on options
my $grid = new Grid::SGE({
			project 	=> 810001,
});

my $job_id = $ARGV[0] || die "Job identifier is required\n";
my @Jobs = $grid->wait_for_jobs($job_id);
foreach my $Job (@Jobs) {
	print $Job->to_string();
	ok($Job->status eq "complete", "positive test (1a)");
}

my $num_tests = scalar @Jobs;
done_testing($num_tests);



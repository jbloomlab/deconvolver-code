#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Cwd;

# Get the user
my $user=`whoami`; chomp $user;
$user ||= "deconvolve";
my $cwd			= cwd();
my $grid_deconvolve = "$Bin/grid_deconvolve.t";
my $sff_file 		= "$Bin/data/test_GEIVOUP02.sff";
my $fastq_file 	= "$Bin/data/test_GEIVOUP02.fastq";
my $pattern_file 	= "$Bin/data/test_barcode_metadata_from_GLK.txt.pat";
my $answers_file 	= "$Bin/data/test_GEIVOUP02_expected_answer_trim_BC019CG.txt.trim";
my $dir = "/usr/local/scratch/$user";
my $key = "tcag";

# (1) Sff file
print STDERR "######################## DECONVOLUTION TEST 1 ############################\n";
my $syscall = "$grid_deconvolve $sff_file $pattern_file $answers_file $dir/test.1.sff > $cwd/log.sff";
system($syscall);

# (2) Sff file with key
print STDERR "######################## DECONVOLUTION TEST 2 ############################\n";
$syscall = "$grid_deconvolve $sff_file $pattern_file $answers_file $dir/test.2.sff $key > $cwd/log.sff.key";
system($syscall);


# (3) Fastq file
print STDERR "\n######################## DECONVOLUTION TEST 3 ############################\n\n";
$syscall = "$grid_deconvolve $fastq_file $pattern_file $answers_file $dir/test.3.fastq > $cwd/log.fastq";
system($syscall);

# (1) Fastq file with key
print STDERR "######################## DECONVOLUTION TEST 4 ############################\n";
$syscall = "$grid_deconvolve $fastq_file $pattern_file $answers_file $dir/test.4.fastq.key $key > $cwd/log.fastq.key";
system($syscall);


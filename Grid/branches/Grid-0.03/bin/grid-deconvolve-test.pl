#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Barcode::Deconvolver;
use Test::More;

=head1 NAME
	
	grid-deconvolve-test.pl -  A command-line interface to run the JCVI
				deconvolution pipeline, and test the results using a 
				manually reviewed answers file.
	
=head1 USAGE
	
	grid-deconvolve-test.pl 
		--project 810001
		--infile seq.fasta
		--pattern t/barcodes.pat
	[	--name myJob
		--mailto naxelrod@jcvi.org
		--mailon abe 
		--notify 0
	     --verbose
	  	--priority 1
	  	--queue fast.q
	  	--trim-points-only
		--o /path/to/out
		--e /path/to/err
		--t /path/to/tmp
		--help 
	 ]

EXAMPLE

	grid-deconvolve-test.pl -P 810001 --infile /usr/local/devel/VIRIFX/Grid/t/data/FTF2AAH01.sff --pattern /usr/local/devel/VIRIFX/Grid/t/data/barcodes.pat -o /usr/local/scratch/naxelrod/out -t /usr/local/scratch/naxelrod/tmp -e /usr/local/scratch/naxelrod/err

=head1 SYNOPSIS
	
	my $grid = new Grid::SGE({
				project => 810001,
				queue	=> "fast.q",
				name	=> "gridDeconvolve",
				tmpdir	=> $tmpdir,
				errdir	=> $errdir,
				outdir	=> $outdir,
				verbose	=> 1,
				poll_delay => 30,
	});
	my $Deconvolver = new Grid::Tools::Barcode::Deconvolver({
				grid	=> $grid,
				pattern	=> $pattern,
				infile	=> $fasta,
				tmpdir	=> $tmpdir,
				errdir	=> $errdir,
				outdir	=> $outdir,
				verbose	=> 1,
				options => "-pmismatch 2 -filter -rformat excel -stdout true -complement Yes"
	});

	my $num_assignments = $Deconvolver->run();

	print $grid->failed_tasks_report() if $grid->failed_tasks;
	
=head1 OPTIONS

B<--project,-p>
	REQUIRED.  Project code id.

B<--name,-n>
	OPTIONAL.  Job name.

B<--queue,-q>
	OPTIONAL.  qsub destination queue.

B<--pattern,-b>
	REQUIRED. Path to barcode pattern file

B<--infile,-f>
	REQUIRED.  Path to read sequence fasta, fastq, or sff file.

B<--outdir,-o>
	OPTIONAL.  Path to output directory.  Directory will be created if not found.

B<--mismatches,-m>
	OPTIONAL.  Number of mismatches to allow in running fuzznuc searches.

B<--mailto>
	OPTIONAL.  Email address to send job status notifications.

B<--mailon>
	OPTIONAL.  Status events to determine when to job notifications.
			 Example: be 
	# a		 mail is sent when the job is aborted by the batch system. 
	# b		 mail is sent when the job begins execution. 
	# e		 mail is sent when the job terminates. 
	
B<--readlength,-l>
	OPTIONAL.  Minimal acceptable read length after barcode trimming is complete.

B<--clamplength,-l>
	OPTIONAL.  The length of the barcode clamp. A/K/A hexamer length.  Default is 6.

B<--key,-k>
	OPTIONAL.  The Roche/454 key sequence.  This is typically a 4bp sequence (eg. TCAG)
			 that is on the 5' end of any read, upstream of the barcode sequence.
			 If the key sequence is provided, we search the key + barcode sequence
			 against the untrimmed read sequences to define the appropriate
			 trim points.  Default is undefined.

B<--trim_points_only>
	OPTIONAL.  Boolean parameter to output only the trim points file for the 
			 sequences for each barcode.

B<--num_seqs>
	OPTIONAL.  Integer parameter to specify the number of sequences per fasta 
			 file for splitting and distributing the Fuzznuc searches.

B<--cleanup>
	OPTIONAL.  Boolean parameter to determine if the temporary files, such as
			 the output of fuzznuc searches, should be removed once the 
			 pipeline is complete.  Default is true.

B<--help,-h>
	Print this message

=head1  DESCRIPTION
	
	This script deconvolves a fasta, fastq or sff file based on the barcode sequences, as 
	determined by running fuzznuc to find the best hits of read to barcode sequences.
	The results are a report of trim points, and a list of barcode fasta files with 
	each entry representing a read sequence that has its unambiguous best hit to 
	the bar code.  The bar code sequences are trimmed by default, unless otherwise 
	specified.

	The rules for trimming reads:
		1) If barcode hits overlap, extend the hit and trim the extended region;
		2) Allow any number of hits, as long as the trimmed read length >= $MIN_READ_LENGTH
		3) Throw out and log any sequences with hits to multiple barcodes

	All results are in space-based coordinates.

=head1  CONTACT
	
	Nelson Axelrod
	naxelrod@jcvi.org

=cut

# Get the user
my $user=`whoami`; chomp $user;
$user ||= "deconvolve";

# Get command-line options
my %opts = ();
my $results = GetOptions( \%opts, 
					'answers|a=s',
					'project|p=s',
					'infile|f=s',
					'pattern=s',
					'name|n=s',
					'queue|q=s',
					'mailto=s',
					'mailon=s',
					'mismatches|m:i',
					'num_seqs:i',
					'poll_delay:i',
					'options|opts=s',
					'tmpdir|t=s',
					'outdir|o=s',
					'errdir|e=s',
					'outformat=s',
					'key|k=s',
					'readlength:i',
					'clamplength:i',
					'cleanup!',
					'notify!',
					'trim_points_only!',
					'verbose!',
					'help|h',
) || &_pod;

# Verify required parameters
&_pod if $opts{help};
&_pod("Error: Input fasta, fastq or sff file(s) are required.") unless defined $opts{'infile'};

my $answers_file = $options{answers};
&_pod("Error: An answers file is required to test deconvolution.") 
	unless defined $answers_file && -r $answers_file;

"$Bin/data/test_GEIVOUP02_expected_answer_trim_BC019CG.corrected.trim";

# Get our SGE object based on options
my $grid = new Grid::SGE( \%opts );

# Get our Fuzznuc object
my $Deconvolver = new Grid::Tools::Barcode::Deconvolver( \%opts );

# Use the grid for deconvolution
$Deconvolver->grid($grid);

# Validate parameters provided
my $error_mssg = $Deconvolver->validate_parameters();
&_pod($error_mssg) if $error_mssg;

# run the deconvolution pipeline
my $num_assignments = $Deconvolver->run();

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

# Issue the test plan when we call done_testing($num_tests)
my $num_tests = 0;

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

sub _pod {
	# display a specific error message if provided
	pod2usage( { -message => shift, -exitval => 0, -verbose => 2, -output => \*STDERR } );
}


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
	my $trim_results = $trim_table->{$seq_id} || warn "Error: No trim results found for sequence $seq_id\n";
	if ($trim_results && ref($trim_results) eq "ARRAY") {
		return @$trim_results;
	} else {
		return ("", "", ""); # [ $clear_start, $clear_end, $reason ]
	}
}


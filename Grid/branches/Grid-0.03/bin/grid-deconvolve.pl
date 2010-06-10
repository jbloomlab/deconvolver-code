#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Barcode::Deconvolver;
use Time::HiRes qw(gettimeofday);

=head1 NAME
	
	grid-deconvolve.pl -  A command-line interface to run the JCVI
				deconvolution pipeline using Sun Grid Engine,
				or optionally without the grid.
	
=head1 USAGE
	
	grid-deconvolve.pl 
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

	grid-deconvolve.pl -P 810001 --infile /usr/local/devel/VIRIFX/Grid/t/data/FTF2AAH01.sff --pattern /usr/local/devel/VIRIFX/Grid/t/data/barcodes.pat -o /usr/local/scratch/naxelrod/out -t /usr/local/scratch/naxelrod/tmp -e /usr/local/scratch/naxelrod/err

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
	
B<--readlength,-r>
	OPTIONAL.  Minimal acceptable read length after barcode trimming is complete.

B<--clamplength,-c>
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

# Get command-line options
my %opts = ();
my $results = GetOptions( \%opts, 
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
					'readlength|r:i',
					'clamplength|c:i',
					'cleanup!',
					'notify!',
					'trim_points_only!',
					'verbose!',
					'help|h',
) || &_pod;

# Verify required parameters
&_pod if $opts{help};
&_pod("Error: Input fasta, fastq or sff file(s) are required.") unless defined $opts{'infile'};

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

# Report any job failures
print $grid->num_failed_tasks, " failed tasks\n";
#print $grid->failed_tasks_report() if $grid->num_failed_tasks;


######################## SUB ROUTINES ############################

sub _pod {
	# display a specific error message if provided
	pod2usage( { -message => shift, -exitval => 0, -verbose => 2, -output => \*STDERR } );
}


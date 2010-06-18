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
	
	deconvolve.pl -  A command-line interface to run the JCVI deconvolution pipeline.
	
=head1 USAGE
	
	deconvolve.pl 
		--infile seq.fasta
		--pattern t/barcodes.pat
	[	--trim-points-only
	     --verbose
		--o /path/to/out
		--e /path/to/err
		--t /path/to/tmp
		--help 
	 ]

EXAMPLE

	deconvolve.pl -P 810001 --infile /usr/local/devel/VIRIFX/Grid/t/data/FTF2AAH01.sff --pattern /usr/local/devel/VIRIFX/Grid/t/data/barcodes.pat -o /usr/local/scratch/naxelrod/out -t /usr/local/scratch/naxelrod/tmp -e /usr/local/scratch/naxelrod/err

=head1 OPTIONS

B<--pattern,-b>
	REQUIRED. Path to barcode pattern file

B<--infile,-f>
	REQUIRED.  Path to read sequence fasta, fastq, or sff file.

B<--outdir,-o>
	OPTIONAL.  Path to output directory.  Directory will be created if not found.

B<--mismatches,-m>
	OPTIONAL.  Number of mismatches to allow in running fuzznuc searches.

B<--readlength,-r>
	OPTIONAL.  Minimal acceptable read length after barcode trimming is complete.

B<--clamplength,-c>
	OPTIONAL.  The length of the barcode clamp. A/K/A hexamer length.  Default is 6.

B<--keylength>
        OPTIONAL.  The length of the key sequence.  Default is 4.

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
					'infile|f=s',
					'pattern=s',
					'mismatches|m:i',
					'options|opts=s',
					'tmpdir|t=s',
					'outdir|o=s',
					'errdir|e=s',
					'outformat=s',
					'key|k=s',
					'keylength:i',
					'readlength|r:i',
					'clamplength|c:i',
					'cleanup!',
					'trim_points_only!',
					'verbose!',
					'help|h',
) || &_pod;

# Verify required parameters
&_pod if $opts{help};
&_pod("Error: Input fasta, fastq or sff file(s) are required.") unless defined $opts{'infile'};

# Get our Fuzznuc object
my $Deconvolver = new Grid::Tools::Barcode::Deconvolver( \%opts );

# Validate parameters provided
my $error_mssg = $Deconvolver->validate_parameters();
&_pod($error_mssg) if $error_mssg;

# run the deconvolution pipeline
$Deconvolver->run();



######################## SUB ROUTINES ############################

sub _pod {
	# display a specific error message if provided
	pod2usage( { -message => shift, -exitval => 0, -verbose => 2, -output => \*STDERR } );
}


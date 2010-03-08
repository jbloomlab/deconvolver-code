#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;
use Grid::Tools::Fuzznuc;


=head1 NAME
	
	grid-fuzznuc.pl -  An API to submit fuzznuc jobs (array)
				    to Sun Grid Engine.
	
=head1 SYNOPSIS
	
	my $Fuzznuc = new Grid::Tools::Fuzznuc({
				fastas 	=> [ "seq1.fa", "seq2.fa" ],
				tmpdir	=> 'tmp',
				outdir	=> 'results',
				options 	=> {
						pmismatch		=> 3,
						complement	=> "Y",
						stdout 		=> "true",
						filter		=> "true",
						sequence		=> "sequences.fasta",
						pattern		=> "[CG](5)TG{A}N(1,5)C",
						rformat		=> "excel"
				}
	});
	
	my $grid = new Grid::SGE({
				project 	=> 810001,
				name		=> "myJob",
				verbose	=> 1,
	});
	
	# Add the Fuzznuc job request
	$grid->add_command($Fuzznuc->command($grid));
	
	# Submit the jobs and wait until they are all complete (success/failure)
	my @Jobs = $grid->submit_and_wait();
	
	# Get the fuzznuc output files
	my $files = $Fuzznuc->get_outfiles($grid);
	
	# Get an iterator of FileIO::FuzznucHit objects
	my $iter = FileIO::FuzznucUtils->hit_iterator($files);
	
	while (my $Hit = $iter->()) {
		print $Hit->to_string();
	}
	
	# Or, you can use the stringified approach to set options
	my $Fuzznuc = new Grid::Tools::Fuzznuc({
				fasta	=> $fasta,
				verbose	=> $verbose,
				options	=> "-pattern @patternfile.pat -pmismatch 2 -filter -rformat excel -stdout true -complement Yes",
	});
	
	
=head1 USAGE
	
	grid-fuzznuc.pl 
		--project 810001
		--fasta seq.fasta
		--pattern t/barcodes.pat
	[	--name myJob
		--mailto naxelrod@jcvi.org
		--mailon abe 
		--notify 0
	     --verbose
	  	--priority 1
	  	--task 1-2:1
		--o /path/to/redirect/stdout
		--e /path/to/redirect/stderr
		--help 
	 ]

EXAMPLE

	grid-fuzznuc.pl 
	
=head1 OPTIONS

B<--project,-p>
	REQUIRED.  Project code id.

B<--name,-n>
	OPTIONAL.  Job name.
	
B<--mailto,-e>
	OPTIONAL.  Email address to send job status notifications.

B<--mailon,-o>
	OPTIONAL.  Status events to determine when to job notifications.
			 Example: be 
	# a		 mail is sent when the job is aborted by the batch system. 
	# b		 mail is sent when the job begins execution. 
	# e		 mail is sent when the job terminates. 
	

B<--help,-h>
	Print this message

=head1  DESCRIPTION

This script makes it easy to submit jobs to Sun Grid Engine.

=head1  CONTACT
	
	Nelson Axelrod
	naxelrod@jcvi.org

=cut

# Get command-line options
my %opts = ();
my $results = GetOptions( \%opts, 
					'fasta|f=s@',
					'project|p=s',
					'name|n=s',
					'mailto|m=s',
					'mailon|on:s',
					'poll_delay:i',
					'options|opts=s',
					'tmpdir|t=s',
					'outdir|o:s',
					'error|e:s',
					'notify!',
					'verbose!',
					'help|h',
) || &_pod;

# Verify required parameters
&_pod if $opts{help};
my $fastas = $opts{'fasta'} || &_pod("Error: Input fasta file(s) are required.");

# Get our SGE object based on options
my $sge = new Grid::SGE( \%opts );

# Get our Fuzznuc object
my $Fuzznuc = new Grid::Tools::Fuzznuc( \%opts );

# Add our fuzznuc job
$sge->add_command($Fuzznuc->command());

# Submit our job the grid and get the return job id
my @Jobs = $sge->submit_and_wait();
foreach my $Job (@Jobs) {
	print $Job->to_string();
}

######################## SUB ROUTINES ############################

sub _pod {
	# display a specific error message if provided
	pod2usage( { -message => shift, -exitval => 0, -verbose => 2, -output => \*STDERR } );
}


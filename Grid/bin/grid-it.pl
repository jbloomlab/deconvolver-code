#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../../";
use Grid::SGE;

=head1 NAME
	
	grid-it.pl - A pure Perl API to submit jobs (or job arrays)
				    to Sun Grid Engine.
	
=head1 SYNOPSIS

USAGE:
	
	grid-it.pl 
		--project 810001
		--command "../t/hello.pl nelson"
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

EXAMPLE:
	
	grid-it.pl -p 810001 -name myJob --verbose -c "hello.pl nelson"
	
	# or to receive emails
	grid-it.pl -p 810001 -name myJob -mailto naxelrod@jcvi.org --mailon abe --verbose -c "hello.pl nelson"

	# Submit a job array
	grid-it.pl -p 810001 -name testarray --verbose -o /scratch/naxelrod/test -c "jobarray.sh <input dir> <output dir>"
	grid-it.pl -p 810001 -name testarray --verbose -t 0-2:1 -c "/home/naxelrod/lib/perl/Grid/t/jobarray.sh /home/naxelrod/lib/perl/Grid/t /home/naxelrod/lib/perl/Grid/t/data/output"
	-m naxelrod@jcvi.org --mailon abe
	
=head1 OPTIONS

B<--project,-p>
	REQUIRED.  Project code id.

B<--command,-c>
	REQUIRED.  Command to execute on the grid

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
my %options = ();
my $results = GetOptions( \%options, 
					'command|c:s',
					'project|P=s',
					'name|n=s',
					'mailto|m=s',
					'mailon|on:s',
					'num_nodes|j:i',
					'tasks|t:s',
					'output|o:s',
					'error|e:s',
					'notify!',
					'verbose!',
					'help|h',
) || &_pod;

# Verify required parameters
&_pod if $options{help};
&_pod("A command statement is required") unless $options{command};

# Get our SGE object based on options
my $sge = new Grid::SGE( \%options );

# Submit our job the grid and get the return job id
my @job_ids = $sge->execute();
foreach my $job_id (@job_ids) {
	print "job_id: $job_id\n";
	print `qstat -j $job_id`;
}

my $delim = "----------------------------------------------------------------------------\n";
print "Checking job status by user...$$sge{user}\n";
print $_->to_string(), "\n" foreach $sge->check_jobs_status_by_user();
print $delim;

#$sge->check_jobs_status($job_id);
#print "Checking job status by job id(s)... $job_id\n";
#print "job ", $_->to_string(), "\n" foreach $sge->check_jobs_status($job_id);

print "Checking grid status...\n";
$sge->check_grid_status();
print "grid jobs status\n";
print "job ", $_->to_string(), "\n" foreach $sge->jobs;
print $delim;
#print "grid nodes status\n";
#print "node ", $_->to_string(), "\n" foreach $sge->nodes;
#print $delim;


######################## SUB ROUTINES ############################

sub _pod {
	# display a specific error message if provided
	pod2usage( { -message => shift, -exitval => 0, -verbose => 2, -output => \*STDERR } );
}


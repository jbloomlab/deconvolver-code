package Grid::Tools::Blast;
use strict;
use Cwd;
use base qw(Grid::Tools::Generic Class::Accessor);
use IO::File;
use File::Basename;

# We can run parallelize based on our input fasta sequences
Grid::Tools::Blast->mk_accessors(qw(
	fasta
	num_seqs
	tmpdir
	outdir
	split_fastas
));

################################################################################
# Package defaults
################################################################################\

my $CWD 			= cwd(); 	
my $EXEC 			= 'blastall';	# path of blast executable
my $NUM_SEQS 		= 100000; 	# number of sequences per fasta file for split

################################################################################

# Initialize new objects in order to set default values
sub new {
	my $self = shift->SUPER::new(@_);
	$self->init();
	return $self;
}

=head2 init()

	Validate and set default values for attributes
=cut

sub init {
	my $self = shift;
	
	# Set the defaults for our globals
	$self->executable($EXEC) unless defined $self->executable;
	$self->num_seqs($NUM_SEQS) unless defined $self->num_seqs;
	
	# Validate our fasta files
	my $fastas = $self->fasta;
	if (ref($fastas) eq "ARRAY") {
		foreach my $fasta (@$fastas) {
			die "Could not open fasta file $fasta for reading.\n"
				unless -r $fasta;
		}
	} else {
		$self->fasta([ $fastas ]);
		die "Could not open fasta file $fastas for reading.\n"
			unless -r $fastas;
		
	}
}

=head2 command()
	
	Prep any data that the blast command requires, and
	return the blast command
	
=cut

sub command {
	my ($self, $grid) = @_;
	
	# The Grid object determines where to store our results
	my $tmpdir = $grid->tmpdir;
	my $outdir = $grid->outdir;
	
	# Give every execution a unique id
	my $guid = $grid->create_guid();
	
	# Split fasta files, and store in tmpdir
	my $split_fastas = $self->split_fasta_files($grid);
	
	# job task for each fasta file
	my $num_tasks = scalar @$split_fastas;
	my $task_option = "1-$num_tasks:1";
	
	# This is used to check that all tasks are complete
	$grid->num_tasks($num_tasks);
	
	# Get our system command
	my $executable_command = $self->executable_command();
	
	# Create our job array shell script that executes our command
	# against the dataset assigned to each task
	my $shellscript = "$tmpdir/grid-blast.$guid.sh";
	sysopen(SCRIPT_FH, $shellscript, O_RDWR|O_CREAT|O_TRUNC, 0755);
	print SCRIPT_FH join("\n", 
		"#!/bin/bash",
		"\n# input and output directories",
		"TMPDIR=\$1",
		"OUTPUT_DIR=\$2",
		"\n# temp directory (on the local grid node)",
		"TMP=\"/tmp/blast.\$USER.\$JOB_ID.\$SGE_TASK_ID\"",
		"\n# get an array of fasta files",
		"files=( `ls \$TMPDIR/split.$guid.*.fasta` )",
		"\n# get the file to work on for this job",
		"task_file=\${files[\$SGE_TASK_ID - 1]}",
		"\n# Execute our command on the file assigned using our SGE_TASK_ID",
		"$executable_command -i \$task_file > \$TMP",
		"FUZZNUC_EXIT_STATUS=\$?",
		"\n# Copy the output from tmp to output directory",
		"mv \$TMP \$OUTPUT_DIR/",
		"exit \$FUZZNUC_EXIT_STATUS "
	);
	close SCRIPT_FH;
	
	# Define our command to execute this shell script
	# shellscript.sh <input_dir> <outdir>
	return "-t $task_option $shellscript $tmpdir $outdir";
}

1;
__END__

=head1 NAME

Grid::Tools::Blast - Makes it easy to run blast jobs on the grid

=head1 SYNOPSIS

	use Grid::Tools::Blast;
	use Grid::SGE;

	my $blast = new Grid::Tools::Blast({
				fastas => @fastas,
				options => {
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
	
	$grid->add_job($blast->job_request());
	
	$grid->execute();


=head1 DESCRIPTION



=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut


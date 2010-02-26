package Grid::Tools::Fuzznuc;
use strict;
use Cwd;
use base qw(Grid::Tools::Generic Class::Accessor);
use IO::File;
use File::Basename;
use FileIO::FuzznucUtils;

# We can run parallelize based on our input fasta sequences
Grid::Tools::Fuzznuc->mk_accessors(qw(
	fasta
	hits
	num_seqs
	tmpdir
	outdir
	outfiles
	split_fastas
));

################################################################################
# Package defaults
################################################################################\

my $CWD 			= cwd(); 	
my $EXEC 			= 'fuzznuc';	# path of fuzznuc executable
my $NUM_SEQS 		= 100000; 	# number of sequences per fasta file for split

# specify any default options
my $OPTIONS = {
	pmismatch		=> 2,
	stdout		=> "true",
	filter		=> "true",
	complement	=> "Yes",
	rformat		=> "excel"
};
################################################################################

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
	
	# Set the defaults for any unspecified options
	my $options = $self->options;
	foreach my $key (keys %$OPTIONS) {
		$self->options($key, $OPTIONS->{$key})
			unless defined $self->{options}{$key};
	}
	
	# These are equired option settings
	$self->options('rformat', 'excel');
	$self->options('filter', 'true');
	$self->options('stdout', 'true');
}

=head2 add_fasta()
	
	Add fasta file to list of files
	Returns an arrayref of file paths
	
=cut
sub add_fasta {
	my ($self, $fasta) = shift;
	
	# Validate our fasta files
	my $fasta_arrayref = $self->fasta || [];
	
	# Verify fasta file is readable
	if (-r $fasta) {
		push @$fasta_arrayref, $fasta;
	} else {
		die "Could not open fasta file $fasta for reading.\n";
	}
	return $fasta_arrayref;
}

=head2 get_options_string()
	
	Same as Grid::Tools::Generic::get_options_string,
	except we need to make sure the pattern is surrounded by quotes
	
=cut

sub get_options_string {
	my $T = shift;
	my $options = $T->options;
	my $command;
	foreach my $key (reverse keys %{ $options }) {
		$command.= "-$key ";
		if (defined $options->{$key}) {
			
			# patterns should be enclosed in quotes
			# pattern "TCACCCCTGA"
			if ($key =~ /pattern/ && $options->{$key} !~ /^['"]/ && $options->{$key} !~ /^\@/) 
			{
				$command.= "\"$$options{$key}\" ";
			} else {
				$command.= "$$options{$key} ";
			}
		}
	}
	return $command;
}

=head2 get_hits_iterator_by_files()
	
	Accepts an arrayref of file paths, and
	returns an iterator of FileIO::FuzznucHit objects
	
=cut
sub get_hits_iterator_by_files {
	my ($self, $files) = @_;
	return FileIO::FuzznucUtils->hit_iterator($files);
}


=head2 get_hits_iterator()
	
	Returns an iterator of FileIO::FuzznucHit objects
	
=cut
sub get_hits_iterator {
	my ($self, $grid) = @_;
	
	# Get all of the output files from our grid fuzznuc searches
	my $outfiles = $self->get_outfiles($grid);
	
	# Return an Iterator of Hit objects
	return $self->get_hits_iterator_by_files($outfiles);
}


=head2 get_hits()
	
	Returns an arrayref of FileIO::FuzznucHit objects
	
=cut

sub get_hits {
	my ($self, $grid) = @_;
	
	# Get all of the output files from our grid fuzznuc searches
	my $outfiles = $self->get_outfiles($grid);
	
	# Build our list of FileIO::FuzznucHit objects
	my @Hits;
	foreach my $file (@$outfiles) {
		print STDERR "Reading $file... ";
		my $HitsByFile = FileIO::FuzznucUtils->parse_fuzznuc_by_csv_file($file);
		print scalar @$HitsByFile, " hits", "\n";
		push @Hits, @$HitsByFile;
	}
	return $self->hits(\@Hits);
}


=head2 get_outfiles()
	
	Returns a list of output files that will be generated
	upon upon successful completion
	
=cut

sub get_outfiles {
	my ($self, $grid) = @_;
	
	# This method can only be run after the grid job has been submitted
	my ($job_id, $user, $num_tasks) = ($grid->job_id, $grid->user, $grid->num_tasks);
	my $outdir = $grid->outdir;
	unless ($job_id) {
		warn "Unable to return location of outfiles without a grid job id.  Have you submitted the job yet?\n";
		return;
	}
	
	my @outfiles;
	for (my $i=1; $i<=$num_tasks; $i++) {
		push @outfiles, "$outdir/fuzznuc.$user.$job_id.$i";
	}
	return $self->outfiles(\@outfiles);
}

=head2 command()
	
	Write a shell script to be executed for each task
	and return the command line string to run it
	
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
	my $shellscript = "$tmpdir/grid-fuzznuc.$guid.sh";
	sysopen(SCRIPT_FH, $shellscript, O_RDWR|O_CREAT|O_TRUNC, 0755);
	
	print SCRIPT_FH join("\n", 
		"#!/bin/bash",
		"\n# input and output directories",
		"TMPDIR=\$1",
		"OUTPUT_DIR=\$2",
		"\n# temp directory (on the local grid node)",
		"TMP=\"/tmp/fuzznuc.\$USER.\$JOB_ID.\$SGE_TASK_ID\"",
		"\n# get an array of fasta files",
		"files=( `ls \$TMPDIR/split.$guid.*.fasta` )",
		"\n# get the file to work on for this job",
		"task_file=\${files[\$SGE_TASK_ID - 1]}",
		"\n# Execute our command on the file assigned using our SGE_TASK_ID",
		"$executable_command -sequence \$task_file > \$TMP",
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

Grid::Tools::Fuzznuc - Makes it easy to run fuzznuc jobs on the grid

=head1 SYNOPSIS

	use Grid::Tools::Fuzznuc;
	use Grid::SGE;

	my $fuzznuc = new Grid::Tools::Fuzznuc({
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
	
	$grid->add_job($fuzznuc->job_request());
	
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


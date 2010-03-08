package Grid::SGE;
use strict;
use base qw(Class::Accessor::Fast);
use Grid::SGE::Status qw/check_grid_status check_jobs_status check_jobs_complete 
	check_jobs_status_by_user _job_by_line _node_by_line/;
use Grid::SGE::Control qw/delete_jobs/;

# Need to find a better way of generating UIDs
#use Data::GUID;
use Time::HiRes qw(gettimeofday);

our $VERSION=0.01;

Grid::SGE->mk_accessors(qw( 
		user
		verbose
		project
		mailto
		mailon
		notify
		tasks
		use_cwd
		status
		commands
		guid
		num_tasks
		poll_delay
		job_id
		job_tasks
));

# Mapping of attributes to command line parameters
my %OPTIONS = (
		'name' 		=> '-N',
		'project'		=> '-P',
		'mailto'		=> '-M',
		'mailon'		=> '-m',
		'tasks'		=> '-t',
		'outdir'		=> '-o',
		'errdir'		=> '-e',
);
my $DEFAULT_DIR = "/usr/local/scratch";

=item new
	Initialize new objects in order to set default values
=cut

sub new {
	my $self = shift;
	$self = $self->SUPER::new(@_);
	$self->init();
	return $self;
}

=item init
	Sets default values for unspecified parameters
=cut

sub init {
	my $self = shift;
	
	# Set the initial status
	$self->status('init');
	
	# Set commands to an empty arrayref
	$self->{commands} = [];
	
	# verbose is off by default
	$self->verbose(0) unless $self->{verbose};
	
	# Set the user unless otherwise provided
	unless (defined $self->user) {
		my $user=`whoami`;chomp $user;
		$self->user($user);
	}
	
	# Set dirs using user or default options
	my $USER_DIR = join("/", $DEFAULT_DIR, $self->user);
	my $tmpdir = ($self->{tmpdir})
				? $self->tmpdir($self->{tmpdir})
				: $self->tmpdir("$USER_DIR/tmp");
	my $errdir = ($self->{errdir})
				? $self->errdir($self->{errdir})
				: $self->errdir("$USER_DIR/err");
	my $outdir = ($self->{outdir})
				? $self->outdir($self->{outdir})
				: $self->outdir("$USER_DIR/out");
	
	# Minimum poll frequency of 10 seconds, default 60 seconds
	$self->poll_delay(60) unless defined $self->poll_delay && $self->poll_delay >= 10;
	
	# Allow users to add a command during construction
	if (defined $self->{command}) {
		$self->add_command($self->{command});
		delete $self->{command};
	}
	
	return $self;
}

=head2 create_guid()
	Quick replaceable hack to create guids
=cut
sub create_guid {
	my $grid = shift;
	my ($secut,$musec) = gettimeofday;
	my ($guid) = sprintf("%010d%06d%05d", $secut, $musec, $$);
	return $grid->guid($guid);
}

=head2 tmpdir()
	Getter and setter for our temp directory
=cut
sub tmpdir {
	my ($self, $dir) = @_;
	# Setter
	if ($dir) {
		$self->{tmpdir} = $dir;
		$self->_valid_write_dir($dir) ||
			die "Unable to write to temp directory $dir\n";
	}
	# Getter
	return $self->{tmpdir};
}

=head2 outdir()
	Getter and setter for our output directory
=cut
sub outdir {
	my ($self, $dir) = @_;
	if ($dir) {
		$self->{outdir} = $dir;
		$self->_valid_write_dir($dir) ||
			die "Unable to write to output directory $dir\n";
	}
	return $self->{outdir};
}
=head2 errdir()
	Getter and setter for our error directory
=cut
sub errdir {
	my ($self, $dir) = @_;
	if ($dir) {
		$self->{errdir} = $dir;
		$self->_valid_write_dir($dir) ||
			die "Unable to write to error directory $dir\n";
	}
	return $self->{errdir};
}

=head2 _valid_write_dir()
	Returns boolean if the directory is writable
=cut
sub _valid_write_dir {
	my ($self, $dir) = @_;
	mkdir $dir unless -e $dir;
	return -w $dir;
}

=head2 add_command()
	
	Add a command string to our commands array
=cut

sub add_command {
	my ($self, @commands) = @_;
	foreach (@commands) {
		push @{ $self->{commands} }, $_ if $_;
	}
}

=head2 add_or_update_job()
	
	Add a Grid::SGE::Job 
	$sge->add_or_update_job( new Grid::SGE::Job ( { job_id => 1341353 } ) );
=cut

sub add_or_update_job {
	my ($self, $J)=@_;
	my $job_key = join(":", $J->job_id, $J->task_id);
	$self->{jobs}{$job_key} = $J
		if $J && $J->isa("Grid::SGE::Job");
}
sub jobs {
	return values %{ shift->{jobs} };
}

=head2 add_or_update_node()
	
	Add a Grid::SGE::Job 
	$sge->add_or_update_node( new Grid::SGE::Node ( { name => "default.q@dell-0-4-8.jcvi.org" } ) );
=cut

sub add_or_update_node {
	my ($self, $N)=@_;
	$self->{nodes}{$N->name} = $N
		if $N && $N->isa("Grid::SGE::Node");
}
sub nodes {
	return values %{ shift->{nodes} };
}

=head2 environment()

	Get and set the environment variables for Grid::SGE. 
	
	my $hashref=$sge->environment(\%vars);

=cut

sub environment {
	my ($self, $hash)=@_;
	my $return;
	foreach my $var (qw[SGE_CELL SGE_EXECD_PORT SGE_QMASTER_PORT SGE_ROOT]) {
		if ($hash->{$var}) {
			$ENV{$var}=$hash->{$var};
		}
		$return->{$var}=$ENV{$var};
	}
	return $return;
}

=head2 perc_complete()

Returns the percentage complete as a floating point number 
between 0 and 100, to one-decimal place.

=cut
sub perc_complete {
	my ($self, $job) = @_;
	my $jobs = $self->jobs;
	my $num_complete = 0;
	foreach my $j (@$jobs) {
		$num_complete++ if $self->status($j) =~ /complete/i;
	}
	
	my $total_jobs = scalar @$jobs;
	my $perc_complete = ($num_complete) ? sprintf("%.1f", (($num_complete/$total_jobs)*100)) : "0.0";
	return $perc_complete;
}


=head2 guess_executable()

Attempt to guess the locations of an executable

=cut
sub guess_executable {
	my ($self, $exec) = @_;
	my $guess_exec = `which $exec`; 
	chomp($guess_exec);
	return $guess_exec;
}

=head2 executable()
	
	Get or set the executables that we will use. This method takes upto two arguments. With no arguments we will try and guess the settings that we need, and if we fail we will die. With a single argument we will return that executable path/program, guess it if we don't know it, and then finally fail. With two arguments we will assume that the second is the location of the executable (incl. path) of the first.
	
	We will also take a reference to a hash as the single argument. In this case, we will use the hash as locations of the executables.
	
	e.g.s:
	
	# using a hash to set all the executables at once (recommended as we don't have to guess anything)
	my $exec={'qsub'=>'/usr/local/bin/qsub', 'qstat'=>'/usr/local/bin/qstat'}
	$sge->exectuable($exec);
	my $pid=$sge->job_id;
	
	
	# guessing all the executables (not recommended)
	$sge->exectuables();
	my $pid=$sge->job_id;
	
	# getting the value for qsub
	my $qsubexec=$sge->executable('qsub');
	
	# setting a single value for qsub only
	my $qsubexec=$sge->executable('qsub', '/usr/local/bin/qsub');


=cut
sub executable {
	my ($self, $exec, $path)=@_;
	
	# first, if it is a reference add it
	# eg. $self->{'execute'}{command} = /path/to/bin
	if (ref($exec) eq "HASH") {
		foreach my $command (keys %$exec) {
			$self->{'execute'}->{$command} = $exec->{$command};
		}
		return $self->{'execute'};
	}
	
	# now if we have both arg/var
	elsif ($exec && $path) {
		$self->{'execute'}->{$exec} = $path;
		return $path;
	}
	
	# now if we only have arg
	elsif ($exec) {
		return $self->{'execute'}->{$exec}
			if $self->{'execute'}->{$exec};
	}
	
	# Otherwise, we need to guess the executable location
	if ($exec && ref($exec) ne "HASH") {
		$self->{'execute'}->{$exec} = $self->guess_executable($exec);
		return $self->{'execute'}->{$exec};
	# or, return the current execute hash
	} else {
		return $self->{'execute'};
	}
}



=head2 name()

Get or set the name of the job used by Grid::SGE. 

=cut

sub name {
	my ($self, $val)=@_;
	if ($val) {
		unless ($val =~ /^[a-zA-Z]/) {
			print STDERR "Name must start with a letter. Name is now a$val\n";
			$val="a".$val;
		}
		if ($val =~ / /) {
			$val =~ s/ /_/g;
			print STDERR "Name can not have spaces in it. Name is now $val\n";
		}
		if (length($val) > 10) { 
			$val=substr($val, 0, 10);
			print STDERR "Name is truncated to 10 letters. Name is now $val\n";
		}
		$self->{'name'}=$val;
	}
	return $self->{'name'};
}

=head2 qsub_command()

	Get the qsub command line string

=cut

sub qsub_command {
	
	my $self = shift;
	
	# Get our qsub executable
	my $qsub_command = $self->executable('qsub') || $self->_dieout('qsub');
	
	# add key-value parameters
	foreach my $tag (keys %OPTIONS) {
		if ($self->{$tag}) {
			$qsub_command .= join(" ", " $OPTIONS{$tag}", $self->{$tag});
		}
	}
	
	# add boolean parameters
	$qsub_command .= " -cwd" if $self->use_cwd;
	$qsub_command .= " -notify" if $self->notify;
	
	return $qsub_command;
}

=head2 submit()

	Execute the qsub command and get the list of jobs

=cut

=head2 submit()
	
	Submit a job and wait for it to complete
	Returns a list of Job objects, where each 
	Job represents either a single Job or a task
	in a Job array.
	
	
=cut

sub submit {
	my ($self, $command) = @_;
	
	# Use the command if provided, or get the list of commands
	my $commands = ($command) ? [ $command ] : $self->commands();
	
	# submit our commands using qsub
	my @Jobs;
	while (scalar @$commands) {
		my $command = join(" ", $self->qsub_command, shift @$commands);
		my $job_id = $self->submit_command($command);
		
		# set our job id
		$self->job_id($job_id);
		
		# Get our Job from doing qstat
		push @Jobs, $self->check_jobs_status($job_id);
	}
	return @Jobs;
}

=head2 submit_and_wait()
	
	Submit a job and wait for it to complete
	Returns a list of Job objects, where each 
	Job represents either a single Job or a task
	in a Job array.
	
	
=cut

sub submit_and_wait {
	my ($self, $command) = @_;
	
	# Use the command if provided, or get the list of commands
	my $commands = ($command) ? [ $command ] : $self->commands();
	
	# submit our commands using qsub
	my @Tasks;
	while (scalar @$commands) {
		my $command = join(" ", $self->qsub_command, shift @$commands);
		my $job_id = $self->submit_command($command);
		
		# set our job id
		$self->job_id($job_id);
		
		# Wait until we can get our Jobs from doing qacct
		push @Tasks, $self->wait_for_tasks($job_id);
	}
	
	# Set our list of tasks for this job
	$self->job_tasks(\@Tasks);
	
	return @Tasks;
}

=head2 num_failed_tasks()
	
	Returns number of failed tasks
	
=cut
sub num_failed_tasks {
	my $failed_tasks = shift->failed_tasks;
	return scalar @$failed_tasks;
}

=head2 failed_tasks()
	
	Returns a list of failed tasks
	
=cut
sub failed_tasks {
	my $self = shift;
	my $Tasks = $self->job_tasks;
	my @Failed_Tasks;
	foreach my $T (@$Tasks) {
		my $status = $T->status;
		push @Failed_Tasks, $T unless $status && $status eq "complete";
	}
	return \@Failed_Tasks;
}

=head2 failed_tasks_report()

	Return a report for each job-task that failed for any run
	
=cut
sub failed_tasks_report {
	my $self = shift;
	my $Tasks = $self->failed_tasks;
	my $report = join("\n",
			"# ============================================================",
			"# Failed Tasks",
			"# ============================================================\n"
	);
	
	foreach my $T (@$Tasks) {
		$report .= $T->to_string();
	}
	return $report;
}

=head2 wait_for_tasks()
	
	Poll qacct for a job and wait for all tasks to complete
	
=cut

sub wait_for_tasks {
	my ($self, $job_id) = @_;
	
	# Get the status of our jobs
	my @Tasks = $self->check_jobs_complete($job_id);
	
	# This is the total number of tasks for this job
	my $num_tasks = $self->num_tasks;
	
	# This is the number of complete jobs (failed/success)
	my $num_tasks_complete = scalar @Tasks;
	
	if ($num_tasks) {
		if ($num_tasks == $num_tasks_complete) {
			return @Tasks;
		}
	} else {
		if ($num_tasks_complete) {
			return @Tasks;
		}
	}
	
	# Otherwise, wait and try again
	sleep($self->poll_delay);
	
	# Recursively call function until complete
	return $self->wait_for_tasks($job_id);
}


=head2 submit_command()

	Execute the qsub command and return the job number

=cut

sub submit_command {
	my ($self, $command) = @_;
	
	# Execute the qsub command and check the return values
	print STDERR "Submitting command...\n$command\n" if $self->verbose;
	my $qsub_results = `$command`;
	print STDERR "Results...\n$qsub_results\n" if $self->verbose;
	
	# Get our Job id
	my $job_id;
	
	# Job array: $qsub_results =~ /Your\sjob\-array\s(\d+)\.1\-/
	if ($qsub_results =~ /Your job (\d+)/i) {
		$job_id = $1;
	
	# Your job-array 5158966.1-2:
	} elsif ($qsub_results =~ /Your\sjob\-array\s(\d+)\.[1-9]\-/) {
		$job_id = $1;
		
	# problem running our qsub command
	} else {
		print STDERR "WARNING: No job id from qsub results: $qsub_results\n";
	}
	
	if ($job_id) {
		# Update the status
		$self->status('running');
	
		# set and return the job id
		return $self->job_id($job_id);
	}
}


=head2 _dieout()

Die nicely, with some kind of warning

=cut

sub _dieout {
 my ($self, $val)=@_;
 if ($val eq "command") {
  print STDERR <<EOF;

  You did not specify a command to run and so we are cowardly quitting.

EOF
 }
 elsif ($val =~ /qsub|qstat|qacct/i) {
  print STDERR <<EOF;

  $0 could not find a $val executable. Please check to make sure that it is in your path, and you are running SGE.

EOF
 }
 else {
  print STDERR <<EOF;

  $0 died out for an unexplained reason. Sorry.

EOF
 }

 exit(-1);
}


=head1 Grid::SGE

Interact with the Sun Grid Engine. This module locates the 
executables (qstat, qsub, etc), sets the level of verbosity, 
and other general commands related to the using the SGE.

=head1 ABSTRACT

Grid::SGE is a suite of modules for interacting with the Sun Grid Engine. 
The base module Grid::SGE handles locating the executables 
and making sure everything works fine. 


=head1 AUTHOR

Based on earlier work by
Rob Edwards
rob@salmonella.org
March 2005
http://search.cpan.org/~linsalrob/ScheduleSGE-0.02/

Redesigned by
Nelson Axelrod
naxelrod@jcvi.org
January 2010

=cut

=head2 new()

Instantiate the object, and preload some data. For example

my $sge = Grid::SGE->new(
			project		=> '810001',
			command		=> 'blastwrapper.sh data.in'
			mailto 		=> 'naxelrod@jcvi.org',
			mailon 		=> 'abe',
			notify 		=> 0,
			verbose		=> 1,
			priority		=> 1,
			out 			=> '/path/to/redirect/stdout',
			err 			=> '/path/to/redirect/stderr',
);


=cut

1;

=head1 Grid::SGE::Status

	Check on the status of the Grid::SGE queues. 
	You should not use this method directly, rather you should use the Grid::SGE method that inherits from this, then all the methods herein are available to you.

=head1 AUTHORS

	Nelson Axelrod,
	naxelrod@jcvi.org
	
	Based on earlier work on Schedule::SGE module by
	Rob Edwards,
	rob@salmonella.org
	3/24/05

=cut

package Grid::SGE::Status;
use strict;
use Exporter;
use Grid::SGE::Node;
use Grid::SGE::Job;
use IPC::Open3;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Grid::SGE Exporter);
@EXPORT_OK = qw(
	check_grid_status 
	check_jobs_status 
	check_jobs_status_by_user
	check_jobs_complete 
	_job_by_line
	_node_by_line
);
our $VERSION = '0.01';

=head2 check_jobs_qacct()

	Gets all tasks for a given job and set
	the exit status of each job-task.

=cut
sub check_jobs_complete {
	my ($G, $job_id) = @_;
	my $qacct = $G->executable('qacct') || $G->_dieout('qacct');
	my (@Jobs, %job);
	
	# Make sure job is found
	open READER, "$qacct -j $job_id 2>&1 |";
	
	# The first line is typically error: job id <$id> is not found
	# until the job is complete/failed
	my $first_line;
	
	# Get all of the jobs
	while (<READER>) {
		chomp;
		if (!$first_line) {
			return if $_ =~ /^error/;
			$first_line = 1;
		}
		
		if ($_ =~ /^=/) {
			if (scalar keys %job) {
				push @Jobs, new Grid::SGE::Job( \%job );
				%job = ();
			}
			next;
		} 
		elsif ($_ =~ /^error/) {
			next;
		}
		else {
			my @fields = split /\s+/, $_;
			my ($key, $value) = @fields;
			$value =~ s/^\s+//; $value =~ s/\s+$//;
			
			if ($key =~ /jobnumber/) {
				$job{'job_id'} = $value;
			} elsif ($key =~ /jobname/) {
				$job{'name'} = $value;
			} elsif ($key =~ /exit_status/) {
				$job{'exit_status'} = $value;
			} elsif ($key =~ /owner/) {
				$job{'user'} = $value;
			} elsif ($key =~ /taskid/) {
				$value = "" if $value =~ /undefined/;
				$job{'tasks'} = $value;
			} elsif ($key =~ /hostname/) {
				$job{'host'} = $value;
			} elsif ($key =~ /start_time/) {
				$job{'start_time'} = join(" ", @fields[1..scalar @fields -1]);
			} elsif ($key =~ /end_time/) {
				$job{'end_time'} = join(" ", @fields[1..scalar @fields -1]);
			} elsif ($key =~ /qsub_time/) {
				$job{'submit_time'} = join(" ", @fields[1..scalar @fields -1]);
			} elsif ($key =~ /failed/) {
				# Job is successful
				if (!$value) {
					$job{'status'} = "complete";
				
				# This is a worrysome case
				# eg. When a fuzznuc job dies, there is no evidence of it failing
				# except this failed row is present with a value of 0.
				# The exit_status=0 (bug in fuzznuc?)
				} else {
					$job{'status'} = "failed";
					$job{'error_reason'} = join(" ", @fields[3..scalar @fields -1])
				}
			} elsif ($key =~ /project/) {
				$job{'end_time'} = $value;
			} elsif ($key =~ /^error/) {
				$job{'error_reason'} .= $_."\n";
			}
		}
	}
	close READER;
	
	# Get last job
	push @Jobs, new Grid::SGE::Job( \%job )
		if scalar keys %job;
		
	return @Jobs;
}



=head2 check_jobs_status_by_user()

	Get the status of the grid jobs submitted by a user.
	Returns a list of Grid::SGE::Job objects with refreshed 
	status info

=cut
sub check_jobs_status_by_user {
	my ($G, $user) = @_;
	my $qstat = $G->executable('qstat') || $G->_dieout('qstat');
	$user ||= $G->user;
	foreach (`$qstat -u $user`) {
		$G->_job_by_line($_);
	}
	return $G->jobs;
}

sub check_jobs_status {
	my ($G, @jobs) = @_;
	my $user = $G->user;
	my $qstat = $G->executable('qstat') || $G->_dieout('qstat');
	my @Jobs;
	foreach (`$qstat -u $user`) {
		my $J = $G->_job_by_line($_);
		next unless $J && $J->isa("Grid::SGE::Job");
		my $is_match = 0;
		foreach my $want_job_id (@jobs) {
			if ($J->job_id == $want_job_id) {
				$is_match = 1;
				push @Jobs, $J;
				$G->add_or_update_job($J);
			}
		}
	}
	return @Jobs;
}


sub check_jobs_status_bak {
	my $G = shift;
	my @jobs = (@_) ? shift : ();
	my $qstat = $G->executable('qstat') || $G->_dieout('qstat');
	my %job;
	foreach (`$qstat -j @jobs`) {
		#$G->_job_by_line($_);
		
		chomp;
		if ($_ =~ /^=/) {
			next;
			if (scalar keys %job) {
				$G->add_or_update_job(new Grid::SGE::Job(%job));
			}
		} 
		else {
			my @fields = split /\:/, $_;
			my ($key, $value) = @fields;
			$value =~ s/^\s+//; $value =~ s/\s+$//;
			
			if ($key =~ /job_number/) {
				$job{'job_id'} = $value;
			} elsif ($key =~ /job_name/) {
				$job{'name'} = $value;
			} elsif ($key =~ /owner/) {
				$job{'user'} = $value;
			} elsif ($key =~ /job\-array tasks/) {
				$job{'tasks'} = "$value:$fields[2]";
			} elsif ($key =~ /^error/) {
				print STDERR "error: $_\n";
				$job{'error'} .= $_."\n";
			}
		}
	}
	
	# Get last job
	$G->add_or_update_job(new Grid::SGE::Job( \%job ))
		if scalar keys %job;
	
	return $G->jobs;
}

=head2 check_grid_status()

Get the status of all grid jobs. This will return a listref
of Job objects

=cut

sub check_grid_status {
	my ($G) = @_;
	my $qstat=$G->executable('qstat') || $G->_dieout('qstat');
	my $type;
	foreach (`$qstat -f`) {
		chomp;
		$_ =~ s/^\s+//; $_ =~ s/\s+$//;
		if ($_ =~ /queuename/) {
			$type = "node"; 
			next;
		} elsif ($_ =~ /PENDING/) {
			$type = "job";
			next;
		}
		next if !$_ || $_ =~ /^\-|^\#/;
		
		if ($type =~ /node/) {
			$G->_node_by_line($_);
		} elsif ($type =~ /job/) {
			$G->_job_by_line($_);
		} else {
			print STDERR "Error: $G does not know how to parse line: $_\n";
		}
	}
	
	return $G->jobs;
}

=head2 _node_by_line()

	Internal use only 
	Returns a Grid::SGE::Node object based on the expected
	string of nodes in a qstat -f command

=cut

sub _node_by_line {
	my ($G, $line) = @_;
	return if !$line || $line =~/^job-ID|^\-/;
	$line =~ s/^\s+//; $line =~ s/\s+$//;
	my ($name, $type, $usage, $load, $arch, $states) = split /\s+/, $line;
	warn "Uknown grid queue type $type on line $line\n"
			if $type =~ /[^BICPTN]/;
	warn "Malformed grid node usage $usage on line $line\n"
			if $usage !~ m/\d+\/\d+/;
	
	# Create our node object
	my $N = new Grid::SGE::Node({
				name			=> $name,
				type			=> $type,
				usage		=> $usage,
				load			=> $load,
				arch			=> $arch,
				states		=> $states
	});
	
	# Add/update our hash of nodes
	if ($N && $N->isa("Grid::SGE::Node")) {
		$G->add_or_update_node($N);
	} else {
		warn "Error: Unable to parse node from line: $line\n";
	}
	return $N;
}

=head2 task_id()

	Internal use only 
	Returns a Grid::SGE::Job object based on the expected
	qstat job string

=cut

sub _job_by_line {
	my ($G, $line) = @_;
	return if !$line || $line =~/^job-ID|^\-/;
	$line =~ s/^\s+//; $line =~ s/\s+$//;
	my ($job_id, $load, $name, $user, $states, $date, $time, $slots, $task_id) = split /\s+/, $line;
	$date .= " $time";
	my $J = new Grid::SGE::Job({
				job_id		=> $job_id,
				task_id		=> $task_id,
				name			=> $name,
				load			=> $load,
				user			=> $user,
				states		=> $states,
				date			=> $date
	});
	warn "Error: Unable to parse job from line: $line\n"
		unless $J && $J->isa("Grid::SGE::Job");
	
	return $J;
}

1;

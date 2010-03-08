package Grid::SGE::Job;
use strict;
use Time::Interval;
use base qw(Class::Accessor::Fast);

Grid::SGE::Job->mk_accessors(qw( 
		job_id
		task_id
		process_id
		name
		user
		priority
		states
		host
		load
		submit_time
		start_time
		end_time
		status
		error_reason
		exit_status
));

our %STATES = (
	"u" 	=> "unknown",
	"a"	=> "alarm",
	"c"	=> "Calendar suspended",
	"s"	=> "suspended",
	"S"	=> "subordinate",
	"d"	=> "disabled",
	"D"	=> "disabled",
	"E"	=> "error",
	"q"	=> "queue",
	"w"	=> "waiting",
	"h"	=> "hold",
	"r"	=> "running",
	"R"	=> "Restarted"
);

sub new {
	my $self = shift;
	$self = $self->SUPER::new(@_);
	$self->status("failed") if $self->exit_status;
	die "A job id is required to create a new instance of $self"
		unless defined $self->job_id;
	
	return $self;
}

sub elapsed_time {
	my ($J, $units) = @_;
	# validate units
	if ($units && $units ne "seconds" && $units ne "minutes") {
		warn "Error: elapsed_time() called with invalid units $units.\n";
		$units = "seconds";
	}
	
	# Get the start and end times in hh:mm:ss
	my (undef, $mon1, $day1, $time1, $yr1) = split /\s+/, $J->start_time;
	my (undef, $mon2, $day2, $time2, $yr2) = split /\s+/, $J->end_time;
	
	# Calculate the time interval between end and start time
	my ($hr1, $min1, $sec1) = split /:/, $time1;
	my ($hr2, $min2, $sec2) = split /:/, $time2;
	my $num_secs1 = convertInterval(
                days		=> 0,
                hours		=> $hr1,
                minutes		=> $min1,
                seconds		=> $sec1,
                ConvertTo     => "seconds"
      );
	my $num_secs2 = convertInterval(
                days		=> 0,
                hours		=> $hr2,
                minutes		=> $min2,
                seconds		=> $sec2,
                ConvertTo     => "seconds"
      );
     return $num_secs2 - $num_secs1;
}

sub state_names {
	my $J = shift;
	my $states = $J->states;
	return unless $states;
	my @states = split //, $J->states;
	my @state_names;
	foreach my $s (@states) {
		my $name = $STATES{$s} || warn "Error: uknown state $s\n";
		push @state_names, $name if $name;
	}
	return join(",", @state_names);
}

sub to_string {
	my $J = shift;
	my $attr_line = "%" . 15 . "s %s\n";
	
	return join("",
		"# ============================================================\n",
		sprintf($attr_line, "jobname", $J->name || ""),
		sprintf($attr_line, "jobnumber", $J->job_id || ""),
		sprintf($attr_line, "taskid", $J->task_id || ""),
		sprintf($attr_line, "owner", $J->user || ""),
		sprintf($attr_line, "hostname", $J->host || ""),
		sprintf($attr_line, "submit_time", $J->submit_time || ""),
		sprintf($attr_line, "start_time", $J->start_time || ""),
		sprintf($attr_line, "end_time", $J->end_time || ""),
		sprintf($attr_line, "elapsed_time", $J->elapsed_time || ""),
		sprintf($attr_line, "status", $J->status || ""),
		sprintf($attr_line, "error_reason", $J->error_reason || ""),
		sprintf($attr_line, "exit_status", $J->exit_status),
	);
}

sub to_csv_string {
	my $J = shift;
	my @fields = ($J->job_id, $J->task_id, $J->load, $J->name, $J->host, 
		$J->user, $J->states, $J->state_names, $J->exit_status, $J->start_time, $J->end_time, $J->error_reason);
	for (my $i=0; $i<scalar @fields; $i++) {
		$fields[$i] = "" unless defined $fields[$i];
	}
	return join("\t", @fields);
}

1;


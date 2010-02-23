package Grid::Tools::Generic;
use strict;
use base qw(Class::Accessor);

# The only requirement to create a new instance of 
# a generic tool is the path to the executable
Grid::Tools::Generic->mk_accessors(qw(
	executable
	command
	num_seqs
	split_fastas
));

# Any Grid::Tool::* has two attributes
# 1) executable 	string	path of executable
# 2) options		hash		any options that the executable takes

sub new {
	# Use Class::Accessor to get our blessed object
	my $self = shift->SUPER::new(@_);
	
	# but make sure the options are set properly
	my $options = $self->options;
	if ($options && ref($options) ne "HASH") {
		$self->{options} = {};
		$self->options($options);
	}
	return $self;
}

=head2 executable_command()
	
	Generate system command to execute this tool with
	prescribed options
	
=cut

sub executable_command() {
	my $self = shift;
	return join(" ", $self->executable, $self->get_options_string);
}

=head2 run()
	
	Runs the executable command on the grid, and waits
	for all jobs to complete
	
	Returns 1 on success, 0 if any task fails
	
=cut

sub run() {
	my ($self, $grid) = @_;
	unless ($grid && $grid->isa("Grid::SGE")) {
		die "Error: $self->run method expects a Grid::SGE object\n";
	}
	
	$grid->add_command($self->command($grid));
	
	# Submit our Fuzznuc jobs to the grid and wait until complete
	my @Jobs = $grid->submit_and_wait();
	
	# Report any failed tasks
	my $is_success = 1;
	foreach my $Job (@Jobs) {
		unless ($Job->status eq "complete") {
			print STDERR $Job->to_string() if $grid->verbose;
			$is_success = 0;
		}
	}
	return $is_success;
}

=head2 get_options_string()
	
	Generate options string for any generic tool
	
=cut

sub get_options_string {
	my $T = shift;
	my $options = $T->options;
	my $command;
	foreach my $key (reverse keys %{ $options }) {
		$command.= "-$key ";
		$command.= "$$options{$key} " if defined $options->{$key};
	}
	return $command;
}

=head2 options()

There are three way to set options

(1) Use a hash to set all options
$Tool->options( { 
		pattern		=> "[CG](5)TG{A}N(1,5)C", 
		mismatches 	=> 2,
		length		=> 10
});

(2) Getters and setters for specifying individual options
my $cutoff 	= $Tool->options('mismatches');
my $new_cutoff = $Tool->options('mismatches', 3);

(3) set options using an options string such as
my %options = $Tool->options("-mismatches 2 -length 10 -pattern [CG](5)TG{A}N(1,5)C");

=cut

sub options {
	my ($self, $option, $value) = @_;
	
	# Set options using hash, and
	# returns the options hash
	if (ref($option) eq "HASH") {
		foreach my $key (keys %$option) {
			$self->{options}{$key} = $option->{$key};
		}
	}
	
	# Set option-value pair, and
	# returns the value
	elsif ($option && $value) {
		$self->{options}{$option} = $value;
		return $value;
	}
	
	# Returns the value for an option
	elsif ($option) {
		
		# typical getter usage
		return $self->{options}->{$option}
			if defined $self->{options}->{$option};
			
		# Eg. $self->options("-mismatches 2 -length 10");
		# Set options using an options string, and
		# returns the options hash
		# -pattern \"[CG](5)TG{A}N(1,5)C\" -pmismatch 2 -filter -rformat excel -stdout true -complement Yes" --verbose 1
		my @opts = split /\s+/, $option;
		if (scalar @opts > 1) {
			for (my $i=0; $i< scalar @opts; $i++) {
				
				# Get the next two parameters
				my ($param1, $param2) = ($opts[$i], $opts[$i+1]);
				
				# Unexpected
				if ($param1 !~ /^-/) {
					die "Error: Unexpected option $param1. $self expects options attributes to begin with a '-' character.\n";
				}
				
				# boolean parameter (i.e. not an option-value pair)
				elsif (!$param2 || $param2 =~ /^-/) {
					$param1 =~ s/^\-+//g;
					$self->{options}{$param1} = undef;
				}
				
				# standard option-value pair
				else {
					$param1 =~ s/^\-+//g;
					$self->{options}{$param1} = $param2;
					$i++;
				}
			}
		}
	}
	
	# returns the options hash, otherwise
	return $self->{options};
}

=head2 command()

	Set the options based on a given options string
	This is the reverse of normal practice, where we use
	our options hash to define the options string.
	
=cut

sub set_options_by_string {
	my ($self, $options_str) = @_;
	my @options = split /\s+/, $options_str;
	while (@options) {
		my ($option, $value) = (shift, shift);
		$option =~ s/^\-+//g; # strip dashes
		$self->{options}{$option} = $value;
	}
}

=head2 split_fasta_files()

	Write temporary fasta files with $self->num_seqs per file
	Used by Grid::Tools::Fuzznuc, and Blast, so moved to Generic class
	
=cut
sub split_fasta_files {
	my ($self, $grid) = @_;
	my $num_seqs = $self->num_seqs;
	my $fasta_files = $self->fasta;
	my $i=0;
	
	my $tmpdir = $grid->tmpdir;
	my $guid = $grid->guid;
	
	# Write temp fasta files with $num_seqs per file
	my @tmp_files;
	my $tmp_fasta = "$tmpdir/split.$guid.$i.fasta";
	open (TMP_FASTA, ">$tmp_fasta") || die "Could not open temporary fasta file $tmp_fasta for writing.\n";
	push @tmp_files, $tmp_fasta;
	
	foreach my $fasta_file (@$fasta_files) {
		# Read in our fasta file
		open(FASTA, $fasta_file) || die "Could not open fasta file $fasta_file for reading.\n";
		
		while (<FASTA>) {
			# Create a new file every $num_seqs
			if ($_=~/^>/) {
				$i++;
				unless ($i % $num_seqs) {
					my $set = int($i/$num_seqs);
					$tmp_fasta = "$tmpdir/split.$guid.$set.fasta";
					#print STDERR "split fasta $tmp_fasta\n";
					close TMP_FASTA;
					open (TMP_FASTA, "> $tmp_fasta") || die "Could not open temporary fasta file $tmp_fasta for writing.\n";
					push @tmp_files, $tmp_fasta;
				}
			}
			print TMP_FASTA $_;
		}
		close FASTA;
	}
	$self->split_fastas(\@tmp_files);
}


1;

__END__
=head1 NAME

Grid::Tools::Generic - Run any executable on the grid

=head1 SYNOPSIS

	use Grid::Tools::Generic;
	use Grid::SGE;

	my $Tool = new Grid::Tools::Generic({
				executable => "/usr/local/devel/mytool.pl",
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
	
	$grid->add_request($Tool->get_request());
	
	$grid->execute();


=head1 DESCRIPTION

	Any Grid::Tool::* has two attributes
	1) executable 	string	path of executable
	2) options	hash		any options that the executable takes

	Any implementation of a tool, just needs to implement
	a get_command() method that does whatever work is required
	to make a grid job(s) for the tool.
	
	The typical implementation of get_command() should prep the input 
	data (eg. split up Fasta files for Blast, Fuzznuc, etc), and 
	possibly write a shell script that the get_command() method calls 
	on the processed input data set.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut


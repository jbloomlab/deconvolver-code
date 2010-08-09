package Grid::Tools::Barcode::Deconvolver;
use strict;
use base qw(Grid::Tools::Generic Class::Accessor);
use Grid::Tools::Barcode::Trimmer;
use Grid::Tools::Fuzznuc;
use File::Basename;
use FileIO::BarcodeUtils;
use FileIO::FastaUtils;
use FileIO::FuzznucUtils;
use Bio::Seq::Quality;
use Time::HiRes qw(gettimeofday);
use Cwd;
use IO::Handle;

######################## DEFAULTS  ##########################

# Default Settings
my $TRIM_READS 	= 1; 					# trim the barcodes off of the read sequences
my $TRIM_POINTS_ONLY = 0;					# output only the trim points
my $MIN_READ_LENGTH = 50; 					# Minimum acceptable read length after barcode trimming
my $NUM_MISMATCHES	= 2;						# Default number of mismatches
my $CLAMP_LENGTH 	= 6; 					# Length of barcode clamp to trim off
my $KEY_LENGTH 	= 4; 					# Length of 454 keys, used for defining the trim points file.
my $CWD 			= cwd(); 					# current working directory
my $TMP_DIR 		= "$CWD/tmp"; 				# temp directory
my $OUT_DIR 		= "$CWD/out";				# output directory
my $MULTIBARCODE_FILE = "report_multibarcode.log";
my $NUM_SEQS_PER_SEARCH = 50000; 				# Number of sequences per fuzznuc search (50k by default)
my $SFFFILE_MODE = 0;						# Run deconvolution in sfffile mode (default: false)
my $CLEANUP 		= 1;						# Boolean value to determine whether we should
										# delete temp files on completion.  True by default.

##############################################################

Grid::Tools::Barcode::Deconvolver->mk_accessors(qw( 
	grid
	infile
	informat
	sequence_file
	base_sequence_file
	fasta_file
	fastq_file
	quals_file
	sff_file
	pattern
	mismatches
	tmpdir
	outdir
	outformat
	barcode_table
	fasta_table
	assignments_table
	multibarcode_table
	barcode_distr_table
	hit_table
	trim_table
	key
	keylength
	trim_points_only 
	verbose
	guid
	num_barcodes
	num_seqs_with_hits
	num_seqs_validated
	sfffile_mode
	cleanup
	logfilehandle
));

sub new {
	my $self = shift->SUPER::new(@_);
	$self->init();
	return $self;
}


=item init
 Sets default values for unspecified parameters and
 creates temp and output dirs
=cut

sub init {
	my ($self, $args) = @_;
	
	# Check for required input parameters
	die "An input file is required." unless defined $self->infile;
	die "A pattern or pattern file is required." unless defined $self->pattern;
	die "Unable to read input file $$self{infile}\n" unless -f $self->infile && -r $self->infile;
	die "Unable to read barcode file $$self{pattern}\n" unless -f $self->pattern && -r $self->pattern;
	
	# Set defaults for any unspecified, optional parameters
	$self->set_defaults();
	$self->validate_parameters();
	
	# Confirm an input sff file is provided when running in sfffile mode
	# die "An input sff file is required to run in sfffile mode.\n"
	#		if $self->sfffile_mode && $self->informat ne "sff";
	
	# Validate the key sequence when provided an sff file
	if ($self->key) {
		if ($self->informat eq "sff") {
			my $key_sequence = $self->get_key_sequence_from_sffinfo();
			die "Error: Input key sequence (", $self->key, ") does not match key sequence ($key_sequence) determined by sffinfo.\n"
				if $self->key ne $key_sequence;
		}
		
	# Set the key automatically using sffinfo when running in sfffile mode
	# and no key sequence is provided
	} else {
		if ($self->sfffile_mode) {
			if ($self->informat eq "sff") {
				$self->key($self->get_key_sequence_from_sffinfo);
			} else {
				die "Error: A key sequence is required to run in sfffile mode for a ".$self->informat." input file.\n";
			}
		}
	}
	
	# Make output directories
	mkdir $self->outdir unless -e $self->outdir;
	mkdir $self->tmpdir unless -e $self->tmpdir;

	# Read in our pattern file to set barcode-specific clamplength and readlength 
	$self->read_pattern_file();
	
	return $self;
}

=item get_key_sequence_from_sffinfo
 
 Uses sffinfo to get the key sequence from an sff file
 
=cut
sub get_key_sequence_from_sffinfo {
	my $self = shift;
	my $infile = $self->infile;
	my $key = `sffinfo $infile | head -n 100 | grep 'Key Sequence:' | cut -d ':' -f 2`;
	$key =~ s/[\s]//g;
	return $key;
}

=item set_defaults
 
 Set default values based on various package variables
 
=cut

sub set_defaults {
	my $self = shift;
	
	# Set our options hash for these variables
	my $clamplength = (defined $self->{clamplength}) ? $self->{clamplength} : $CLAMP_LENGTH;
	my $readlength = (defined $self->{readlength}) ? $self->{readlength} : $MIN_READ_LENGTH;
	delete $self->{clamplength};
	delete $self->{readlength};
	$self->clamplength($clamplength);
	$self->readlength($readlength);
	
	# Specify the number of mismatches
	my $mismatches = (defined $self->{mismatches}) ? $self->{mismatches} : $NUM_MISMATCHES;
	$self->mismatches($mismatches);
	
	# Run in sfffile mode
	$self->sfffile_mode($SFFFILE_MODE) unless defined $self->sfffile_mode;
	
	# Use the default directories unless specified
	$self->outdir($OUT_DIR) unless $self->outdir;
	$self->tmpdir($OUT_DIR) unless $self->tmpdir;
	
	# Use fully qualified paths for input files and directories
	$self->outdir($CWD ."/". $self->outdir) if $self->outdir && $self->outdir !~ /^\//;
	$self->tmpdir($CWD ."/". $self->tmpdir) if $self->tmpdir && $self->tmpdir !~ /^\//;
	$self->infile($CWD ."/". $self->infile) if $self->infile && $self->infile !~ /^\//;
	$self->pattern($CWD ."/". $self->pattern) if $self->pattern && $self->pattern !~ /^\//;
	
	$self->trim_points_only($TRIM_POINTS_ONLY) unless defined $self->trim_points_only;
	$self->fasta_table({}) unless defined $self->fasta_table;
	$self->barcode_distr_table({}) unless defined $self->barcode_distr_table;
	$self->multibarcode_table({}) unless defined $self->multibarcode_table;
	$self->assignments_table({});
	$self->trim_table({});
	$self->num_seqs_with_hits(0);
	
	$self->cleanup($CLEANUP) unless defined $self->cleanup;
	my $is_verbose = (defined $self->verbose && $self->verbose) ? 1 : 0;
	$self->verbose($is_verbose);
	
	# Set key sequence (use uppercase as convention)
	if (defined $self->key) {
		die "Deconvolution requires either providing the key sequence or length, but not both.\n"
			if $self->keylength;
		$self->key(uc($self->key));
		$self->keylength(length($self->key));
	
	# Otherwise, we set our default key length
	} else {
		$self->keylength($KEY_LENGTH) unless defined $self->keylength;
	}
	
	# Log results to STDOUT unless a log filehandle is provided
	$self->logfilehandle(*STDOUT) unless defined $self->logfilehandle;
	
}

=head2 validate_parameters()

	Verify the input files are readable and output dirs are writeable
=cut
sub validate_parameters {
	my $self = shift;
	my $error_mssg;
	$error_mssg.= "Error: Unable to read input file $$self{infile}\n" unless -f $self->infile && -r $self->infile;
	$error_mssg.= "Error: Unable to read barcode file $$self{pattern}\n" unless -f $self->pattern && -r $self->pattern;
	$error_mssg.= "Error: Unable to write to output directory $$self{outdir}\n" unless -w $self->outdir;
	$error_mssg.= "Error: Unable to write to temp directory $$self{tmpdir}\n" unless -w $self->tmpdir;
	$error_mssg.= "Error: Invalid output format $$self{outformat}.  Choose fasta or fastq as the output format.\n"
					if $self->outformat && $self->outformat ne "fasta" && $self->outformat ne "fastq" && $self->outformat ne "sff";
	
	# validate input file, and define its format
	my @suffixlist = ("fastq", "fasta", "fa", "sff");
	if ($self->infile) {
		my ($name, $path, $ext) = fileparse($self->infile, @suffixlist);
		$ext = "fasta" if $ext eq "fa";
		if ($ext eq "fasta" || $ext eq "fastq" || $ext eq "sff") {
			$self->informat($ext);
			$self->outformat($ext) unless defined $self->outformat;
			
			if ($self->outformat eq "sff" && $self->informat ne "sff") {
				$error_mssg.= "Error: SFF output format requires an SFF input file.\n"
			}
			
		} else {
			$error_mssg.= "Error: Unrecognized extension ($ext) for input file $$self{infile}.  Input file must have a .fa, .fasta, .fastq or .sff extension.\n"
		}
	}
	
	# verify fuzznuc binary 
	my $guess_exec = `which fuzznuc`;
	$error_mssg.= "Error: Failed test attempt to guess path of fuzznuc executable\n"
			unless $guess_exec;
	
	return $error_mssg;
}

=head2 grid()
	
	Getter & Setter of SGE::Grid object
	
=cut
sub grid {
	my ($self, $grid) = @_;
	if ($grid) {
		if ($grid->isa("Grid::SGE")) {
			$self->{grid} = $grid;
		} else {
			print STDERR "Expecting an SGE::Grid object but got $grid instead\n";
		}
	}
	return $self->{grid};
}
=head2 fuzznuc()
	
	Getter & Setter of Grid::Tools::Fuzznuc object
	
=cut
sub fuzznuc {
	my ($self, $fuzznuc) = @_;
	if ($fuzznuc) {
		if ($fuzznuc->isa("Grid::Tools::Fuzznuc")) {
			$self->{fuzznuc} = $fuzznuc;
		} else {
			print STDERR "Expecting a Grid::Tools::Fuzznuc object but got $fuzznuc instead\n";
		}
	}
	return $self->{fuzznuc};
}

=head2 print_runtime_settings()

	Return the runtime settings for this program
	(useful for debugging purposes)
	
=cut
sub print_runtime_settings {
	my $self = shift;
	my $key = (defined $self->key) ? $self->key : " ";
	my $mismatches = (defined $self->mismatches) ? $self->mismatches : " ";
	return join("\n",
			"# ============================================================",
			"# Program: $0",
			"# RunTime: ". scalar localtime(),
			"# RunParameters: ",
			"#\t-infile $$self{infile}",
			"#\t-pattern $$self{pattern}",
			"#\t-informat $$self{informat}",
			"#\t-mismatches $mismatches",
			"#\t-readlength ".$self->readlength,
			"#\t-clamplength ".$self->clamplength,
			"#\t-keylength ".$self->keylength,
			"#\t-key $key",
			"#\t-trim_points_only $$self{trim_points_only}",
			"#\t-tmpdir $$self{tmpdir}",
			"#\t-outdir $$self{outdir}",
			"#\t-outformat $$self{outformat}",
			"#\t-informat $$self{informat}",
			"#\t-outformat $$self{outformat}",
			"# ============================================================\n"
	);	
}

=head2 readlength()

	(1) Set readlength options
	$self->readlength({ 
			DA3DFDA13		=> 3, 
			DA3DFDA25 	=> 5,
	});
	
	(2) Getters and setters for specifying individual options
	my $readlength = $Tool->readlength('DA3DFDA13');
	my $readlength = $Tool->readlength('DA3DFDA13', 3);
	my $default_readlength = $Tool->readlength();
	
=cut

sub readlength {
	return shift->_options('readlength', @_);
}

=head2 clamplength()
	Same approach as readlength
=cut
sub clamplength {
	return shift->_options('clamplength', @_);
}

=head2 _options()
	
	General mechanism used to set hash-based attributes,
	used for the purpose of enabling barcode-specific options
	
=cut

sub _options {
	my ($self, $attribute, $option, $value) = @_;
	
	# This is only used to manage barcode-specific readlength and clamplength values
	if ($attribute !~ /readlength|clamplength/) {
		die "Error: unexpected use of the _options function for attribute $attribute\n";
	}
	
	# Set options using hash, and returns the options hash
	if (ref($option) eq "HASH") {
		foreach my $key (keys %$option) {
			$self->{$attribute}{$key} = $option->{$key};
		}
	}
	
	# Set option-value pair, and returns the value
	elsif ($option && defined $value) {
		$self->{$attribute}{$option} = $value;
		return $value;
	}
	
	# Returns the value for an option (typical getter use)
	elsif (defined $option) {
		
		# standard getter, ie. $self->clamplength('BC009CG');
		if (defined $self->{$attribute}{$option}) {
			return $self->{$attribute}{$option};
		
		# Setter case: Setting default values
		# eg. $self->clamplength(5), ie. $self->_options('clamplength', 5);
		} elsif ($option =~ /\d/ && $option !~ /\D/) { # ints only
			$self->{$attribute}{'DEFAULT'} = $option;
			return $option;
		}
	}
	
	# returns the default value for the attribute
	return $self->{$attribute}{'DEFAULT'};
}


=item set_fuzznuc
 
 Creates our Grid::Tools::Fuzznuc object using the 
 temporary fasta file created during write_temp_fasta_files()
 
=cut

sub set_fuzznuc {
	my $self = shift;
	
	# Set our fuzznuc options
	my $options = $self->options;
	
	# Add the barcode pattern to our fuzznuc options
	$self->options('pattern', "\@$$self{pattern}");
	
	# Add the user-specified number of mismatches, if provided
	$self->options('pmismatch', $self->mismatches)
		if defined $self->mismatches;
	
	# Requires csv output
	$self->options('rformat', 'excel');
	$self->options('filter', 'true');
	
	# sfffile does not search the complementary direction
	if ($self->sfffile_mode) {
		$self->options('complement', 'No');
	}
	
	# Create our Grid::Tools::Fuzznuc object 
	# (uses the directory structure defined during set_defaults)
	my $Fuzznuc = new Grid::Tools::Fuzznuc({
			fasta	=> $self->fasta_file,
			options 	=> $self->options,
			outdir 	=> $self->outdir,
			tmpdir 	=> $self->tmpdir,
	});
	
	return $self->fuzznuc($Fuzznuc);
}

=item run
 
 Runs the deconvolution pipeline
 Returns number of sequence assignments on success, 0 on failure
 
=cut
sub run {
	my $self = shift;
	my $grid = $self->grid;
	
	# Validate our input parameters
	my $errorstr = $self->validate_parameters;
	die $errorstr if $errorstr;
	
	# Verify that grid job is run from a submit host
	if ($grid) {
		my $qstat = $grid->executable('qstat') || $grid->_dieout('qstat');
		if (system($qstat)) {
			die "Error: Attempt to run deconvolution on the grid from an invalid submit host.";
		}
	}
	
	# Log runtime settings
	print $self->print_runtime_settings();
	
	# Generate Fasta sequence and quality files
	$self->write_temp_fasta_files;
	
	# Run Fuzznuc searches, and get an iterator of FileIO::FuzznucHit objects
	my ($is_success, $iter);
	if ($grid) {
		
		# Run Fuzznuc on the grid
		$is_success = $self->run_grid_fuzznuc();
		
		# exit run pipeline if grid submission fails
		return $is_success if !$is_success;
		
		# Get an iterator of FileIO::FuzznucHit objects from our output files
		$iter = $self->fuzznuc->get_hits_iterator($grid, $self->verbose);
		
	# Otherwise, run a system command
	} else {
		
		# Run fuzznuc, returns 1 or 0 on success, and output file
		my $output_file;
		($is_success, $output_file) = $self->run_fuzznuc();
		
		# Return a hit iterator, from our fuzznuc results file
		$iter = $self->fuzznuc->get_hits_iterator_by_files([ $output_file ]);
	}
	
	# deconvolve if the fuzznuc searches completed successfully
	$self->deconvolve($iter) if $is_success;
	
	# Returns the number of sequence assignments
	return $self->num_assignments;
	
}

sub cleanup_files {
	my $self = shift;
	my $grid = $self->grid;
	return unless $grid;
	
	# Get dirs and job info
	my ($outdir, $tmpdir) = ($grid->outdir, $grid->tmpdir);
	my ($guid, $job_id, $user) = ($grid->guid, $grid->job_id, $grid->user);

	# Get the generic name for each of the job tasks
	my $Jobs = $grid->job_tasks;
	my $job_name = $Jobs->[0]->name if defined $Jobs->[0];
	
	# remove outdir/gridDeconvolve.o* files
	`rm -f $outdir/$job_name.o*` if $job_name;
	
	# Delete fuzznuc files
	`rm -f $outdir/fuzznuc.$user.$job_id.*`;
	
	if ($guid) {
		# Delete shell script
		`rm -f $outdir/grid-fuzznuc.$guid.sh`;
		
		# Delete split files
		`rm -f $outdir/split.$guid.*.fasta`;
	}
}

=head2 deconvolve()

	Runs the deconvolution pipeline, given an iterator of hits
	This is separated out in order to more easily swap out the
	search implementation (i.e. a better search implementation?)
	
	Returns number of sequence assignments (0 is failure)
	
=cut

sub deconvolve {
	my ($self, $iter) = @_;
	
	# Store the fasta sequences in a hash table
	my $fasta_table = $self->make_fasta_table;
	
	# Assign sequences to barcodes, and make our multicode table
	my $assignments_table = $self->make_assignment_table($iter);
	
	# Log our multicoded sequences
	$self->write_multicode_report();
	
	# Writes barcode fasta files
	$self->write_barcode_fasta;
	
	# Sanity check the results
	$self->validate_results;
	
	# print our log report on completion
	$self->print_log_report;
	
	# Cleanup after grid jobs
	$self->cleanup_files() if $self->cleanup;
	
	# Returns number of sequence assignments
	return scalar keys %$assignments_table;
}

=head2 run_grid_fuzznuc()

	Runs Fuzznuc searches on the grid using our
	Grid::Tools::Fuzznuc object
	
	Returns 1 on success, 0 if any task fails
	
=cut
sub run_grid_fuzznuc_and_wait {
	my $self = shift;
	print STDERR "Running fuzznuc on the grid\n" if $self->verbose;
	return $self->set_fuzznuc->run($self->grid);
}

=head2 run_grid_fuzznuc()

	Runs Fuzznuc searches on the grid using our
	Grid::Tools::Fuzznuc object
	
	Returns 1 on success, 0 if fails.
	Success is defined by a defined $grid->job_id
=cut
sub run_grid_fuzznuc {
	my $self = shift;
	print STDERR "Running fuzznuc on the grid\n" if $self->verbose;
	return $self->set_fuzznuc->run($self->grid);
	#return $self->set_fuzznuc->run_without_waiting($self->grid);
}


=head2 run_fuzznuc()

	Run Fuzznuc searches in parallel, if a Grid::SGE object is provided,
	or in serial otherwise.
=cut

sub run_fuzznuc {
	my ($self, $output_file) = @_;
		
	# Run on the grid if possible
	return $self->run_grid_fuzznuc() if $self->grid;
	
	# Set a unique name for our output file
	unless ($output_file) {
		my ($secut,$musec) = gettimeofday;
		my ($guid) = sprintf("%010d%06d%05d", $secut, $musec, $$);
		$self->guid($guid);
		$output_file = $self->outdir."/fuzznuc_$guid.csv";
	}
	
	# Get our fuzznuc executable command, and redirect stdout to our output file
	my $exec_fuzznuc = $self->set_fuzznuc->executable_command() 
			. " -sequence ".$self->fasta_file." | grep -v SeqName > $output_file";
	
	# Run fuzznuc and check the status code that is returned
	print STDERR "Running fuzznuc (not on the grid)\n$exec_fuzznuc\n" if $self->verbose;
	my $status_code = system($exec_fuzznuc);
	
	# fuzznuc should return a status code of 0, if it completes successfully
	if ($self->verbose) {
		if ($status_code) {
			print STDERR "\nError: fuzznuc returns with a status code $status_code\n";
		} else {
			print STDERR "fuzznuc is complete.\n";
		}
	}
	
	my $is_success = ($status_code) ? 0 : 1;
	return ($is_success, $output_file);
}


sub write_temp_fasta_files {
	my $self = shift;
	
	# Set path in tmp directory where we will write our fasta file
	my ($base_file) = fileparse($self->infile, qr/\.[^.]*/);
	my $fasta_file = $self->fasta_file($self->tmpdir."/$base_file.fasta");
	my $fastq_file = $self->fastq_file($self->tmpdir."/$base_file.fastq");
	my $quals_file = $self->quals_file($self->tmpdir."/$base_file.quals");
	my $infile = $self->infile;
	
	if ($self->informat eq "fasta") {
		
		# write a clean fasta file
		my $sc = "sed -s 's/\:/\_/g' $infile > $fasta_file";
		die "Error: Problem with writing clean fasta file from $infile\n" 
				if system($sc);
		
		# Prepend the key sequence to the fasta file
		$self->pre_and_post_append_fasta_with_key($fasta_file, $self->key) if $self->key;
		
	} elsif ($self->informat eq "fastq") {
	
		# write a clean fastq file
		print STDERR "Writing a clean fastq file... $fasta_file\n" if $self->verbose;
		my $sc_fastq = "sed -s 's/\:/\_/g' $infile > $fastq_file";
		die "Error: Problem with writing clean fastq file from $infile\n" 
				if system($sc_fastq);
		
		# Don't need quals file if we only want trim points
		print STDERR "Converting fastq to fasta... $fasta_file\n" if $self->verbose;
		my $num_seqs = ($self->trim_points_only) 
			? FileIO::FastaUtils->write_fasta_from_fastq($fastq_file, $fasta_file)
			: FileIO::FastaUtils->write_fasta_from_fastq($fastq_file, $fasta_file, $quals_file);
		
		die "Error: No sequences parsed from input fastq file $fastq_file\n" 
				unless $num_seqs;
	
		# Prepend the key sequence to the fasta file
		$self->pre_and_post_append_fasta_with_key($fasta_file, $self->key) if $self->key;
		
	} elsif ($self->informat eq "sff") {
		
		my $sff_file = $self->sff_file($self->tmpdir."/$base_file.sff");
		
		# (1) copy sff file to temp directory
		print STDERR "Copying sff input file to temporary directory\n" if $self->verbose;
		my $sc_sff = "cp $infile $sff_file";
		die "Error: Problem with copying input sff file to $sff_file\n" 
				if system($sc_sff);
		
		# (2) convert to fasta sequence file
		print STDERR "Writing fasta sequences from sff file\n" if $self->verbose;
		
		# Use the untrimmed sequence file (with key sequence prepended) for our searches
		my $sc_fasta = "sffinfo -s $sff_file > $fasta_file";
		die "Error: Problem with writing fasta from sff file $sff_file\n" 
				if system($sc_fasta);
		
		# We do not do the following since we need the lowercase sequence to be trimmed, and then we
		# can append the key sequence
		# ($self->key) ? "sffinfo -notrim -s $sff_file > $fasta_file"
		
		# Prepend the key sequence to the fasta file
		$self->pre_and_post_append_fasta_with_key($fasta_file, $self->key) if $self->key;
		
		# Qualities files are not required to define trim points
		unless ($self->trim_points_only) {
			# (3) convert to fasta qualities file
			print STDERR "Writing fasta qualities from sff file\n" if $self->verbose;
			my $sc_fasta_quals = "sffinfo -q $sff_file > $quals_file";
			die "Error: Problem with writing fasta qualies from sff file $sff_file\n" 
					if system($sc_fasta_quals);
		}
	}
}

=head2 pre_and_post_append_fasta_with_key()

	Takes a fasta file and prepends the key sequence to it, and replaces
	the fasta file
	
=cut
sub pre_and_post_append_fasta_with_key {
	my ($self, $fasta_file, $key) = @_;
	
	# Read fasta sequences into tmp_table
	my $tmp_table = FileIO::FastaUtils->parse_fasta_by_file($fasta_file);
	                                          
	# Verify that the key sequence is not already prepended
	my $is_key_prepended = 1;
	foreach my $F (values %$tmp_table) {
		if (substr($F->seq, 0, length($key)) !~ /$key/i) {
			$is_key_prepended = 0;
		}
	}
	
	# We do not expect the key to be prepended, but we can handle this situation 
	# simply by not prepending it
	if (!$is_key_prepended) {
		
		# Get the reverse complement of our key
		my $revcomp_key = $self->revcomp($key);
		
		# Write a new fasta file using this tmp_table
		open (PREPENDED_FASTA, "> $fasta_file") || die "Could not open fasta file for writing\n";
		foreach my $F (values %$tmp_table) {
			$F->seq($key.$F->seq.$revcomp_key);
			print PREPENDED_FASTA join(" ", ">".$F->id, $F->desc), "\n", $F->fasta_seq;
		}
		close PREPENDED_FASTA;
	}
}

sub revcomp {
	my ($self, $seq) = @_;
	my $revcomp_seq = reverse($seq);
	$revcomp_seq =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp_seq;
}

=head2 read_pattern_file()

	Reads in a barcode file, as a generic fasta file
	This is the standard format used for Fuzznuc
=cut

sub read_pattern_file {
	my $self = shift;
	my $pattern_file = $self->pattern;
	my $barcode_table = $self->barcode_table;
	
	# We read in the pattern file to do the following
	# 	- store barcodes in barcode_table lookup
	#	- Check and/or set barcode-specific values for 
	# 		- readlength
	#		- clamplength
	open(IN, "< $pattern_file") || die "Could not open fasta file $pattern_file for reading.\n";
	
	# We output a temporary pattern file with prepended key sequences if provided (in tmpdir)
	my ($base_file) = fileparse($pattern_file, qr/\.[^.]*/);
	my $tmp_pattern_file = $self->pattern($self->tmpdir."/$base_file.pat");
	open(OUT, "> $tmp_pattern_file") || die "Could not open fasta file $tmp_pattern_file for writing.\n";
	if ($self->verbose) {
		print STDERR "Reading barcode file $pattern_file\n";
		print STDERR "Writing temp barcode file to $tmp_pattern_file\n";
	}
	
	# delimit sequences by header line
	local $/ = "\n>";
	
	while (<IN>) {
		chomp;
		my ($header, @sequence_lines) = split /\n/, $_;
		
		# strip initial > char
		$header = substr($header,1) if $header=~/^>/; 
		
		# Separate the identifier from the rest of the header line
		my ($id, $desc);
		if ($header =~ /^\s*(\S+)\s*(.*)/) {
			($id, $desc) = ($1, $2);
			
			# Check for any parameters in the description line
			my @attributes = split /\s/, $desc;
			foreach (@attributes) {
				$_ =~ s/[\<\>]//g; # remove <> characters
				my ($key, $value) = split /\=/, $_;
				if ($key eq "readlength" || $key eq "clamplength") {
					$self->_options($key, $id, $value);
				}
			}
			
			# We need to do this to solve the problem of integer barcode ids
			# Set the clamplength and readlength to the default values for this barcode
			$self->{'clamplength'}{$id} = $self->{'clamplength'}{'DEFAULT'}
				unless defined $self->{'clamplength'}{$id};
			$self->{'readlength'}{$id} = $self->{'readlength'}{'DEFAULT'}
				unless defined $self->{'readlength'}{$id};
		}
		
		# Prepend the key to the barcode sequence
		$sequence_lines[0] = $self->key . $sequence_lines[0] if $self->key;
		
		# Output the header information to our temporary pattern file
		# Always include the mismatch parameter in the barcode file
		$header .= " <mismatch=".$self->mismatches.">" unless $header =~ /mismatch\=\d+/;
		print OUT join("\n", ">$header", @sequence_lines), "\n";
		
		# Get the sequence
		my $seq = join("", @sequence_lines);
		$seq =~ s/\s//g; # strip whitespace and newlines
		
		if (exists $barcode_table->{$id}) {
			# Get the object and set its quality values
			my $F = $barcode_table->{$id};
			$F->seq($seq);
		} else {
			# Add a new FileIO::Fasta object to our hash table
			$barcode_table->{$id} = new FileIO::Fasta({
					id => $id,
					desc => $desc,
					seq => $seq,
			});
		}
	}
     
	# Store number of barcodes
	$self->num_barcodes( scalar keys %{ $barcode_table } );
	$self->barcode_table($barcode_table);
	
     close IN;
     close OUT if $self->key;
     
     return $barcode_table;
}

=head2 write_pattern_file_from_glk()

	To be deprecated.
	
	This function reads in barcode data using the output format
	generated by the GLK, and writes a pattern file required by Fuzznuc
=cut

sub write_pattern_file_from_glk {
	my ($self, $glk_barcode_file) = @_;
	my $pattern_file = $self->pattern_file($self->tmpdir."/barcode.pat");
	my $barcode_table = $self->barcode_table(FileIO::BarcodeUtils->parse_barcode_by_file($glk_barcode_file));
	open(FH_PATTERN, "> $pattern_file") || die "Could not open barcode file $pattern_file for writing.\n";
	foreach my $B (values %{ $barcode_table }) {
		print FH_PATTERN ">" . $B->id . "\n" . $B->seq."\n";
	}
}


=head2 make_hit_table()
	
	Build our hash table of Fuzznuc hits
	
=cut

sub make_hit_table {
	my ($self, $iter) = @_;
	my $hit_table = $self->hit_table;
	while (my $Hit = $iter->()) {
		push @{ $hit_table->{$Hit->seq_id} }, $Hit;
	}
	$self->hit_table($hit_table);
}

=head2 make_fasta_table()
	
	Build our hash table of Sequences from fasta file
=cut

sub make_fasta_table {
	my $self = shift;
	
	# Read fasta sequences
	print STDERR "Reading fasta file $$self{fasta_file}\n" if $self->verbose;
	FileIO::FastaUtils->parse_fasta_by_file($self->fasta_file, $self->fasta_table);
	
	# Read fasta quality values
	unless ($self->trim_points_only || !(-f $self->quals_file)) {
		print STDERR "Reading quals file $$self{quals_file}\n" if $self->verbose;
		FileIO::FastaUtils->parse_quals_by_file($self->quals_file, $self->fasta_table);
	}
	
	return $self->fasta_table;
}

=head2 print_log_report()
	
	Return a report of some statistics/numbers on our barcode hits
=cut

sub print_log_report {
	my $self = shift;
	
	# Log to stdout unless a log filehandle is defined
	my $logfh = $self->logfilehandle || *STDOUT;
	
	# Get our hash table of barcode-sequence distribution
	my $barcode_distr_table = $self->barcode_distr_table;
	my $num_seqs_deconvolved = $barcode_distr_table->{total}{seqs} || 0;
	my $deconvolved_bp = $barcode_distr_table->{total}{bp} || 0;
	delete $barcode_distr_table->{total};
	
	# Header with basic stats on our search results
	my $multibarcode_table = $self->multibarcode_table;
	my $num_seqs = ($self->fasta_table) ? scalar keys %{ $self->fasta_table } : 0;
	my $num_barcodes = $self->num_barcodes;
	my $num_seqs_with_hits = $self->num_seqs_with_hits || 0;
	my $num_seqs_validated = $self->num_seqs_validated || 0;
	my $num_multicoded_seqs = ($multibarcode_table) ? scalar keys %{ $multibarcode_table } : 0;
	my $perc_seqs_with_hits = ($num_seqs) ? sprintf("%.2f", (($num_seqs_with_hits/$num_seqs)*100)) : "0.0";
	my $perc_seqs_deconvolved = ($num_seqs) ? sprintf("%.2f", (($num_seqs_deconvolved/$num_seqs)*100)) : "0.0"; 
	my $perc_seqs_validated = ($num_seqs_deconvolved) ? sprintf("%.2f", (($num_seqs_validated/$num_seqs_deconvolved)*100)) : "0.0";
	my $perc_multicoded_seqs = ($num_seqs) ? sprintf("%.2f", (($num_multicoded_seqs/$num_seqs)*100)) : "0.0";
	
	print $logfh join("\n",
		"# Barcodes: $num_barcodes",
		"# Sequences: $num_seqs",
		"# Sequences with barcode hits: $num_seqs_with_hits ($perc_seqs_with_hits\%)",
		"# Sequences with multiple barcodes: $num_multicoded_seqs ($perc_multicoded_seqs\%)",
		"# Sequences successfully deconvolved: $num_seqs_deconvolved ($perc_seqs_deconvolved\%)",
		"# Sequences successfully validated: $num_seqs_validated of $num_seqs_deconvolved ($perc_seqs_validated\%)",
	), "\n";
	
	# Sort our barcodes by the number of sequences in descending order
	print $logfh "# Distribution of barcodes: <barcode_id> <num_sequences> <perc_sequences> <num_bp> <percent_bp>\n";
	for my $barcode_id ( sort { $barcode_distr_table->{$b}{num_seqs} <=> $barcode_distr_table->{$a}{num_seqs} } keys %$barcode_distr_table) {
		my $num_seqs = $barcode_distr_table->{$barcode_id}{num_seqs};
		my $num_bp = $barcode_distr_table->{$barcode_id}{num_bp};
		my $perc_seqs = ($num_seqs_deconvolved) ? sprintf("%.1f", (($num_seqs/$num_seqs_deconvolved)*100)) : "0.0";
		my $perc_bp = ($deconvolved_bp) ? sprintf("%.1f", (($num_bp/$deconvolved_bp)*100)) : "0.0";
		print join(" ", "barcode", $barcode_id, $num_seqs, $perc_seqs, $num_bp, $perc_bp), "\n";
	}
}

=head2 write_multicode_report()
	
	Return a report of some statistics/numbers on our barcode hits
=cut

sub write_multicode_report {
	my $self = shift;
	my $multicode_file = $self->outdir."/report_multicode.log";
	open(FH, ">$multicode_file") || die "Could not open multicode file $multicode_file for writing.\n";
	
	# Report each multicoded sequence
	my $multibarcode_table = $self->multibarcode_table;
	foreach my $seq_id (keys %$multibarcode_table) {
		my @hits = @{ $multibarcode_table->{$seq_id} };
		my $report = "$seq_id\t";
		foreach my $H (@hits) {
			my $hit_loc = $H->min."..".$H->max.":".$H->strand;
			$report .= "\t".join(" ", $H->pattern, $hit_loc);
		}
		print FH $report, "\n";
	}
	close FH;
}

=head2 check_multicoded_hits()
	
	returns 1 or 0 dependending on whether the sequences have hits to 
	more than one barcode
=cut

sub check_multicoded_hits {
	my ($self, @Hits) = @_;
	
	# Sort the hits to a given sequence by number of mismatches
	my @Hits_sorted = sort { $a->num_mismatches <=> $b->num_mismatches } @Hits;
	
	# Set the first hit as our best hit
	my $BestHit = $Hits_sorted[0];
	
	# Assume its not, and test if it is
	my $is_multicoded = 0;
	
	# Get all of the hits to the barcode with the fewest mismatches
	my (@BestHits, $is_ambiguous);
	for (my $i=1; $i<scalar @Hits_sorted; $i++) {
		
		# Add this hit our list
		my $Hit = $Hits_sorted[$i];
		
		if ($BestHit->pattern ne $Hit->pattern) {
			$is_multicoded = 1;
			last;
		}
	}
	return ($is_multicoded, \@Hits_sorted);
}

=head2 make_assignment_table()
	
	Assign sequences to barcodes, and identify
	any multi-barcode sequences
	
	returns $self->assignment_table, a hash table
=cut
sub make_assignment_table {
	
	my ($self, $Hit_Iterator) = @_;
	print STDERR "Assigning sequences to barcodes...\n"
		if $self->verbose;
	
	die "Expecting a hit iterator, but got $Hit_Iterator instead.\n"
		unless $Hit_Iterator && ref($Hit_Iterator) eq "CODE";
	
	# Iterate through our hits, which we assume are sorted by sequence
	my ($seq_id, $prev_id, @Hits, $num_hits, $num_seqs_with_hits);
	while (my $Hit = $Hit_Iterator->()) {
		
		$seq_id = $Hit->seq_id;
		if ($prev_id) {
			
			# Build a list of hits for this sequence
			if ($prev_id eq $seq_id) {
				push @Hits, $Hit;
				
			# Reached a new sequence
			} else {
				$num_seqs_with_hits++;
				
				# Handle the previous set of Hits (to $prev_id)
				$self->assign_sequence_by_hits($prev_id, @Hits) if scalar @Hits;
				
				# Start a new list of hits for this current sequence
				@Hits = ($Hit);
				$prev_id = $seq_id;
			}
		
		# Start a list of hits for this new sequence
		} else {
			# $num_seqs_with_hits++;
			@Hits = ($Hit);
			$prev_id = $seq_id;
		}
	}
	
	# Handle the last hit
	if ($prev_id) {
		$num_seqs_with_hits++;
		$self->assign_sequence_by_hits($prev_id, @Hits) if scalar @Hits;
	}
	
	# Report count of multicoded sequences
	if ($self->verbose) {
		my $num_multibarcode_seqs = scalar keys %{ $self->multibarcode_table };
		print STDERR "Sequences with multi barcode hits: $num_multibarcode_seqs\n";
	}
	
	# Keep track of number of sequences with hits and assignments
	$self->num_seqs_with_hits($num_seqs_with_hits);
	
	# Set our multibarcode table
	$self->multibarcode_table();
	
	# Set and return our assignments table
	return $self->assignments_table;
}

sub assign_sequence_by_hits {
	my ($self, $prev_id, @Hits) = @_;
	
	if ($prev_id) {
		my $num_hits = scalar @Hits;
		
		# Assign sequence to its unique barcode hit
		if ($num_hits == 1) {
			$self->assignments_table->{$Hits[0]->pattern}{$prev_id} = [ $Hits[0] ];
			
		# Otherwise, check if this sequence has multiple barcode hits
		} else {
			my ($is_multicoded, $Hits_Sorted) = $self->check_multicoded_hits(@Hits);
			if ($is_multicoded) {
				$self->multibarcode_table->{$prev_id} = $Hits_Sorted;
			
			# Assign the sorted hits
			} else {
				$self->assignments_table->{$Hits[0]->pattern}{$prev_id} = $Hits_Sorted;
			}
		}
	}
}

=head2 num_assignments()
	
	Calculates the number of sequences assigned to barcodes
	
=cut
sub num_assignments {
	my $self = shift;
	my $assignments_table = $self->assignments_table;
	my $num_assignments;
	foreach my $pattern (keys %{ $assignments_table }) {
		$num_assignments += scalar keys %{ $assignments_table->{$pattern} };
	}
	return $num_assignments;
}

=head2 write_barcode_fasta()
	
	Writes a Fasta file of sequences for each barcode
	
	Additionally, writes the trim reports if requested
	
=cut

sub write_barcode_fasta {
	my $self = shift;
	my $assignments_table = $self->assignments_table;
	my $barcode_table = $self->barcode_table;
	my $fasta_table = $self->fasta_table;
	my $trim_table = $self->trim_table;
	my ($outdir, $outformat) = ($self->outdir, $self->outformat);
	my $key_length = $self->keylength || 0;
	
	# Trim our sequences, and write the trim logfile
	print STDERR "Writing barcode reports\n";
	
	# Create a file of untrimmable reads
	my $untrimmed_file = "$outdir/report_untrimmable.log";
	open(FH_UNTRIMMED, "> $untrimmed_file") || die "Could not open file $untrimmed_file for writing.\n";
	
	# Record the distribution of barcodes across our sequences (num_reads, num_bp per barcode)
	my $barcode_distr_table = $self->barcode_distr_table;
	
	# Use the barcode table to produce a directory for ALL barcodes, even if no assignments
	#foreach my $barcode_id (sort keys %$assignments_table) {
	foreach my $barcode_id (sort keys %$barcode_table) {
	
		# Put all of our trim files in this directory
		my $barcodedir = "$outdir/$barcode_id";
		mkdir $barcodedir unless -e $barcodedir;
		
		# Create a fasta file of barcodes
		# <barcode id>_<barcode sequence>_<num allowed mismatches>_<blinded_id>.fna
		my $barcode_fasta_file = "$barcodedir/$barcode_id.fasta";
		
		# Create a trim report file
		my $trim_file = "$barcodedir/$barcode_id.trim";
		open(FH_TRIMPOINTS, "> $trim_file") || die "Could not open file $trim_file for writing.\n";
		
		# write fasta output directly, or use Bio::SeqIO for writing fastq output
		my $seqio_fastq;
		unless ($self->trim_points_only) {
			print STDERR "Writing barcode results $barcode_fasta_file\n";
			if ($outformat eq "fastq") {
				$seqio_fastq = new Bio::SeqIO( -format => "fastq", -file => ">$barcodedir/$barcode_id.fastq" );
			} else {
				open(FH_BARCODE_FILE, "> $barcode_fasta_file") || die "Could not open results file $barcode_fasta_file for writing.\n";
			}
		}
		
		# Get the sequences assigned to this barcode (if any)
		my @seq_ids = (exists $assignments_table->{$barcode_id}) ? keys %{ $assignments_table->{$barcode_id} } : ();
		
		foreach my $seq_id (@seq_ids) {
			
			# Get the hits for this barcode and sequence pair
			my $Hits = $assignments_table->{$barcode_id}{$seq_id};
			
			# Get the read sequence from our fasta lookup table
			die "Bug: No sequence in fasta table for id: $seq_id\n" unless defined $fasta_table->{$seq_id};
			my $F = $fasta_table->{$seq_id};
			my $seq = $F->seq; # untrimmed sequence
			
			# Get the clear range of our sequence (extracts barcode and clamp)
			my ($clear_start, $clear_end, $reason) = Grid::Tools::Barcode::Trimmer->trim_clear_range($self, $barcode_id, $seq, $Hits);
			
			# Log the reason for throwing out this read if we can't define where the clear range start begins
			if (!defined $clear_start) {
				print FH_UNTRIMMED "$seq_id $reason\n";
				next;
			}
		
                        # Store the results in our trim table
                        $trim_table->{$seq_id} = [ $barcode_id, $clear_start + 1, $clear_end, $reason ];
	
			# Throw out the read if it does not fulfill our minimum length requirements
			my $clear_seqlen = $clear_end - $clear_start;
			if ($self->readlength && $clear_seqlen < $self->readlength) {
				print FH_UNTRIMMED "$seq_id BELOW_MIN_LENGTH $clear_seqlen\n";
				next;
			}
			
			# Report the distribution of barcodes across sequences
			$barcode_distr_table->{$barcode_id}{num_seqs}++;
			$barcode_distr_table->{$barcode_id}{num_bp} += $clear_seqlen;
			$barcode_distr_table->{total}{seqs}++;
			$barcode_distr_table->{total}{bp} += $clear_seqlen;
			
			# Write trimpoints in residue coordinates (add 1 to start position)
			print FH_TRIMPOINTS join("\t", $seq_id, $clear_start + 1, $clear_end), "\n";
			next if $self->trim_points_only;
			
			# Since the clear range coordinates are offset by the key_length 
			# (in order to produce the correct trim points file), then we need
			# to offset the coordinates by the key_length to obtain the correct 
			# sequence and quality values.
			my ($qv_start, $qv_end) 	= ($clear_start - $key_length, $clear_end - $key_length);
			my ($seq_start, $seq_end)= ($clear_start - $key_length, $clear_end - $key_length);
			
			# The above comments are true except when the sequence is prepended
			# with the key which only occurs when $self->key is defined
			# Quality values are never prepended with the key
			if ($self->key) {
				($seq_start, $seq_end) = ($clear_start, $clear_end);
			}
			
			# Update our fasta record with the "clear range" trimmed sequence
			$F->seq(substr($seq, $seq_start, $seq_end - $seq_start));
			die "Bug: No sequence", join("\n", $seq_id, $seq, "trimmed_seq clear range $clear_start..$clear_end seq range $seq_start..$seq_end ", $F->seq, $F->qual), "\n" if !$F->seq;
			
			# Update our quality values with the "clear range" trimmed values
			if (defined $F->qual) {
				my @quality_arr = split / /, $F->qual; 
				$F->qual(join(" ", @quality_arr[$qv_start..$qv_end - 1]));
			}
			
			# Report the locations of our barcode hits
			my $locs_string;
			my @sorted_hits = sort { $a->{min} <=> $b->{min} } @$Hits;
			foreach my $Hit (@sorted_hits) {
				$locs_string .= $Hit->min."..".$Hit->max.":".$Hit->strand.",";
			}
			$locs_string =~  s/\,+$//; # Remove trailing comma
			
			# After trimming is complete, we can write out the barcode-trimmed sequence report
			if ($outformat eq "fastq") {
				
				# Use Bio::Seq::Quality to write out fastq
				my $SeqWithQuality = new Bio::Seq::Quality( 
										-id => $F->id,
										-seq => $F->seq,
										-qual => $F->qual,
										-desc => ""
				);
				$seqio_fastq->write_fastq($SeqWithQuality);
				
			} else {
				my $desc = join(" ", "barcode=$barcode_id", "length=$clear_seqlen", $F->desc); # "fuzznuc_barcode_hits=$locs_string" 
				print FH_BARCODE_FILE join(" ", ">".$F->id, $desc), "\n", $F->fasta_seq;
				print FH_BARCODE_FILE join(" ", ">".$F->id, $desc), "\n", $F->fasta_qual if defined $F->qual;
			}
		}
		close FH_BARCODE_FILE;
		
		# Output a trimmed sff file
		if ($outformat eq "sff" && !$self->trim_points_only) {
			
			# Use our input sff file (-i $sff_file) and trim file (-t $trim_file)
			my $sff_file = $self->sff_file; # $self->infile 
			
			# Write a trimmed sff file (-o $trim_sff_file)
			my $trim_sff_file = "$barcodedir/$barcode_id.sff";
			
			# Produce a list of sequence ids to trim (-i $seqid_file)
			my $seqid_file = "$barcodedir/$barcode_id.ids";
			my $sc_seqidfile = "cut -f 1 $trim_file > $seqid_file";
			`$sc_seqidfile`; my $status_code = $?;
			#my $status_code = system($sc_seqidfile);
			print STDERR "Error: Failed system call to write seqid file with status code $status_code\n"
				if $status_code;			
			
			# Write the trimmed sff files using our trim file and sff file (-t $barcode_fasta_file $sff_file)
			my $sc_sff = "sfffile -o $trim_sff_file -i $seqid_file -t $trim_file $sff_file";
			print STDERR "Writing trimmed sff...\n$sc_sff\n" if $self->verbose;
			`$sc_sff`; $status_code = $?;
			#$status_code = system($sc_sff);
			print STDERR "Error: sfffile returns with status code $status_code while generating trimmed sff for barcode $barcode_id\n$sc_sff\n"
				if $status_code;
		}

	}
	
	my $num_seqs_deconvolved = $barcode_distr_table->{total}{seqs};
	return $num_seqs_deconvolved;
}

=head2 validate_results()
	
	Verifies that all sequence-to-barcode assignments have a corresponding
	hit in the Fuzznuc search output files.
	
=cut

sub validate_results {
	my $self = shift;
	my $outdir = $self->outdir;
	print STDERR "Validating results of deconvolution...\n" if $self->verbose;
	
	# Read in fuzznuc results and build a table of sequences-barcodes that is unique
	my %seqs;
	
	opendir (OUTDIR, $outdir) || die "Could not open directory $outdir\n";
	my @entries = readdir(OUTDIR);
	
	my $fuzznuc_file_prefix = ($self->grid) ? join(".", "fuzznuc", $self->grid->user, $self->grid->job_id) : "fuzznuc_" . $self->guid . ".csv";
	foreach my $fuzznuc_file (@entries) {
		if ($fuzznuc_file =~ /^$fuzznuc_file_prefix/) {
			open (FUZZNUC, "< $outdir/$fuzznuc_file") || die "Could not open file $fuzznuc_file\n";
			while (<FUZZNUC>) {
				chomp;
				next if $_=~/SeqName/;
				my ($seq_id, $start, $end, $length, $strand, $barcode_header, $num_mismatches) = split /\t/, $_;
				
				# Get barcode id, and exclude any other key-value pairs
				my ($barcode) = split /\s/, $barcode_header;
				
				# Split up barcode
				$seqs{$barcode}{$seq_id}++;
			}
			close FUZZNUC;
		}
	}
	my $is_fuzznuc_results_found = ( scalar keys %seqs ) ? 1 : 0;
	print STDERR "Validation of results failed because no fuzznuc results were found in $outdir\n" unless $is_fuzznuc_results_found;	
	
	# Read in results and validate the assignments
	my $num_invalid = 0; my $num_valid = 0;
	
	# TODO: iterate through barcodes here
	foreach my $barcode_dir (sort keys %{ $self->barcode_table }) {
		if (-d "$outdir/$barcode_dir") {
			next if $barcode_dir eq "." || $barcode_dir eq "..";
			opendir (BARCODE_DIR, "$outdir/$barcode_dir") || die "Could not open barcode directory $outdir/$barcode_dir\n";
			my @trim_files = readdir(BARCODE_DIR);
			foreach my $trim_file (@trim_files) {
				if ($trim_file =~ /\.trim/) {
					open (ASSIGNMENTS, "< $outdir/$barcode_dir/$trim_file" ) || die "Could not open trim file $outdir/$barcode_dir/$trim_file\n";
					
					# get barcode from the name of the file
					my $barcode = $trim_file;
					$barcode =~ s/\.trim//g;
					
					while (<ASSIGNMENTS>) {
						chomp;
						my ($seq_id) = split /\s+/, $_;
						if ($seq_id) {
							if (!exists $seqs{$barcode}{$seq_id}) {
								$num_invalid++;
							} else {
								$num_valid++;
							}
						}
					}
					close ASSIGNMENTS;
				}
			}
			close BARCODE_DIR;
		}
	}
	
	$self->num_seqs_validated($num_valid);
	if ($num_invalid && $self->verbose) {
		print STDERR "Failed validation with number of invalid assignments: $num_invalid\n";
	} else {
		print STDERR "Complete validation with valid:$num_valid invalid:$num_invalid sequences.\n";
	}
	
	my $is_success = ($is_fuzznuc_results_found && $num_valid && !$num_invalid) ? 1 : 0;
	return $is_success;
}

1;
__END__

=head1 NAME

Grid::Tools::Barcode::Deconvolver

=head1 SYNOPSIS

	use Grid::Tools::Barcode::Deconvolver;
	use Grid::SGE;
 	 
	# Get our SGE object based on options
	my $grid = new Grid::SGE({
				project 	=> 810001,
				name		=> "myDeconvolveJob",
				tmpdir	=> $tmpdir,
				errdir	=> $errdir,
				outdir	=> $outdir,
				verbose	=> $verbose,
				poll_delay => 30,
	});

	my $Deconvolver = new Grid::Tools::Barcode::Deconvolver({
				grid		=> $grid,
				pattern	=> $pattern,
				infile	=> $fasta,
				tmpdir	=> $tmpdir,
				outdir	=> $outdir,
				options 	=> "-pattern \@$pattern -pmismatch 2 -filter -rformat excel -stdout true -complement Yes"
	});
	
	my $is_success = $Deconvolver->run();
	
	die $grid->failed_tasks_report() if !$is_success;
	
	
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

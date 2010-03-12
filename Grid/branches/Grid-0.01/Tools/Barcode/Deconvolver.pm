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

######################## DEFAULTS  ##########################

# Default Settings
my $TRIM_READS 	= 1; 					# trim the barcodes off of the read sequences
my $TRIM_POINTS_ONLY = 0;					# output only the trim points
my $MIN_READ_LENGTH = 50; 					# Minimum acceptable read length after barcode trimming
my $CLAMP_LENGTH 	= 6; 					# Length of barcode clamp to trim off
my $KEY_LENGTH 	= 0; 					# Length of 454 keys, used for defining the trim points file.
my $CWD 			= cwd(); 					# current working directory
my $TMP_DIR 		= "$CWD/tmp"; 				# temp directory
my $OUT_DIR 		= "$CWD/out";				# output directory
my $MULTIBARCODE_FILE = "report_multibarcode.log";
my $NUM_SEQS_PER_SEARCH = 50000; 				# Number of sequences per fuzznuc search (50k by default)
my $CLEANUP 		= 1;						# Boolean value to determine whether we should
										# delete temp files on completion.  True by default.

##############################################################

Grid::Tools::Barcode::Deconvolver->mk_accessors(qw( 
	grid
	fuzznucs
	infile
	informat
	sequence_file
	base_sequence_file
	fasta_file
	fastq_file
	quals_file
	sff_file
	pattern
	tmpdir
	outdir
	outformat
	clamplength 
	keylength 
	min_readlength
	barcode_table
	fasta_table
	assignments_table
	multibarcode_table
	barcode_distr_table
	hit_table
	trim_points_only 
	verbose
	guid
	num_seqs_with_hits
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
	
	# Set defaults for any unspecified, optional parameters
	$self->set_defaults();
	
	# Make output directories
	mkdir $self->outdir unless -e $self->outdir;
	mkdir $self->tmpdir unless -e $self->tmpdir;
	
	return $self;
}

=item set_defaults
 
 Set default values based on various package variables
 
=cut

sub set_defaults {
	my $self = shift;
	my $cwd = cwd();
	$self->clamplength($CLAMP_LENGTH) unless defined $self->clamplength;
	$self->keylength($KEY_LENGTH) unless defined $self->keylength;
	$self->min_readlength($MIN_READ_LENGTH) unless defined $self->min_readlength;
	$self->outdir($OUT_DIR) unless defined $self->outdir;
	$self->tmpdir($TMP_DIR) unless defined $self->tmpdir;
	$self->trim_points_only($TRIM_POINTS_ONLY) unless defined $self->trim_points_only;
	$self->fasta_table({}) unless defined $self->fasta_table;
	$self->barcode_distr_table({}) unless defined $self->barcode_distr_table;
	$self->num_seqs_with_hits(0);
	$self->cleanup($CLEANUP) unless defined $self->cleanup;
	
	# Log results to STDOUT unless a log filehandle is provided
	$self->logfilehandle(*STDOUT) unless defined $self->logfilehandle;
	
}

=head2 validate_parameters()

	Verify the input files are readable and output dirs are writeable
=cut
sub validate_parameters {
	my $self = shift;
	my $error_mssg;
	$error_mssg.= "Error: Unable to read input file $$self{infile}\n" unless -r $self->infile;
	$error_mssg.= "Error: Unable to read barcode file $$self{pattern}\n" unless -r $self->pattern;
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
	return join("\n",
			"# ============================================================",
			"# Program: $0",
			"# RunTime: ". scalar localtime(),
			"# RunParameters: ",
			"#\t-infile $$self{infile}",
			"#\t-pattern $$self{pattern}",
			"#\t-informat $$self{informat}",
			"#\t-min_readlength $$self{min_readlength}",
			"#\t-clamplength $$self{clamplength}",
			"#\t-keylength $$self{keylength}",
			"#\t-trim_points_only $$self{trim_points_only}",
			"#\t-tmpdir $$self{tmpdir}",
			"#\t-outdir $$self{outdir}",
			"#\t-outformat $$self{outformat}",
			"#\t-grid ".ref($$self{grid}),
			"#\t-informat $$self{informat}",
			"#\t-outformat $$self{outformat}",
			"# ============================================================\n"
	);	
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
	
	# Requires csv output
	$self->options('rformat', 'excel');
	$self->options('filter', 'true');
	
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
	
	# Read our barcode file (stores results in $self->barcode_table)
	$self->read_pattern_file;
	
	# Run Fuzznuc searches, and get an iterator of FileIO::FuzznucHit objects
	my ($is_success, $iter);
	if ($grid) {
		
		# Run Fuzznuc on the grid
		$is_success = $self->run_grid_fuzznuc();

		# exit run pipeline if grid submission fails
		return $is_success if !$is_success;
		
		# Get an iterator of FileIO::FuzznucHit objects from our output files
		$iter = $self->fuzznuc->get_hits_iterator($grid, $self->verbose);
		# $iter = $self->fuzznuc->get_hits_pump_iterator($grid, $self->verbose);
		
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
	my ($outdir, $tmpdir) = ($grid->outdir, $grid->tmpdir, $grid->errdir);
	my ($guid, $job_id, $job_name, $user) = ($grid->guid, $grid->job_id, $grid->name, $grid->user);
	
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
	
	# Read in our barcode file
	my $barcode_table = $self->read_pattern_file;
	
	# Store the fasta sequences in a hash table
	my $fasta_table = $self->make_fasta_table;
	
	# Assign sequences to barcodes, and make our multicode table
	my $assignments_table = $self->make_assignment_table($iter);
	
	# Log our multicoded sequences
	$self->write_multicode_report();
	
	# Writes barcode fasta files
	$self->write_barcode_fasta;
	
	# print our log report on completion
	print $self->print_log_report;
	
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
	my $self = shift;
		
	# Run on the grid if possible
	return $self->run_grid_fuzznuc() if $self->grid;
	
	# Set a unique name for our output file
	my ($secut,$musec) = gettimeofday;
	my ($guid) = sprintf("%010d%06d%05d", $secut, $musec, $$);
	my $output_file = $self->outdir."/fuzznuc_$guid.csv";
	print STDERR "Running fuzznuc (not on the grid)\n" if $self->verbose;
		
	# Get our fuzznuc executable command, and redirect stdout to our output file
	my $exec_fuzznuc = $self->executable_command() . " > $output_file";
	
	# Run fuzznuc and check the status code that is returned
	my $status_code = system($exec_fuzznuc);
	
	# fuzznuc should return a status code of 0, if it completes successfully
	if ($status_code) {
		print STDERR "\nError: fuzznuc returns with a status code $status_code\n";
	} else {
		print STDERR "ok\n";
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
	
	# return if -f $fasta_file && -f $quals_file;
	
	if ($self->informat eq "fasta") {
		
		# write a clean fasta file
		my $sc = "sed -s 's/\:/\_/g' $infile > $fasta_file";
		die "Error: Problem with writing clean fasta file from $infile\n" 
				if system($sc);
		
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
	
	} elsif ($self->informat eq "sff") {
		
		my $sff_file = $self->sff_file($self->tmpdir."/$base_file.sff");
		
		# (1) copy sff file to temp directory
		print STDERR "Copying sff input file to temporary directory\n" if $self->verbose;
		my $sc_sff = "cp $infile $sff_file";
		die "Error: Problem with copying input sff file to $sff_file\n" 
				if system($sc_sff);
		
		# (2) convert to fasta sequence file
		print STDERR "Writing fasta sequences from sff file\n" if $self->verbose;
		my $sc_fasta = "sffinfo -s $sff_file > $fasta_file";
		die "Error: Problem with writing fasta from sff file $sff_file\n" 
				if system($sc_fasta);
		
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

=head2 read_pattern_file()

	Reads in a barcode file, as a generic fasta file
	This is the standard format used for Fuzznuc
=cut

sub read_pattern_file {
	my $self = shift;
	my $pattern_file = $self->pattern;
	print STDERR "Reading barcode file $pattern_file\n" if $self->verbose;
	return $self->barcode_table(FileIO::FastaUtils->parse_fasta_by_file($pattern_file));
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
	unless ($self->trim_points_only) {
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
	
	# Header with basic stats on our search results
	my $multibarcode_table = $self->multibarcode_table;
	my $num_seqs = ($self->fasta_table) ? scalar keys %{ $self->fasta_table } : 0;
	my $num_barcodes = ($self->barcode_table) ? scalar keys %{ $self->barcode_table } : 0;
	my $num_seqs_with_hits = $self->num_seqs_with_hits;
	my $num_multicoded_seqs = ($multibarcode_table) ? scalar keys %{ $multibarcode_table } : 0;
	my $perc_seqs_with_hits = ($num_seqs) ? sprintf("%.2f", (($num_seqs_with_hits/$num_seqs)*100)) : "0.0";
	my $perc_multicoded_seqs = ($num_seqs) ? sprintf("%.2f", (($num_multicoded_seqs/$num_seqs)*100)) : "0.0";
	
	print $logfh join("\n",
		"# Number of barcode entries: $num_barcodes",
		"# Number of fasta entries: $num_seqs",
		"# Sequences with barcode hits: $num_seqs_with_hits ($perc_seqs_with_hits\%)",
		"# Sequences with multiple barcodes: $num_multicoded_seqs ($perc_multicoded_seqs\%)",
	), "\n";
	
	# Get our hash table of barcode-sequence distribution
	my $barcode_distr_table = $self->barcode_distr_table;
	my ($total_seqs, $total_bp) = ($barcode_distr_table->{total}{seqs}, $barcode_distr_table->{total}{bp});
	delete $barcode_distr_table->{total};

	# Sort our barcodes by the number of sequences in descending order
	print $logfh "# Distribution of barcodes: <barcode_id> <num_sequences> <perc_sequences> <num_bp> <percent_bp>\n";
	for my $barcode_id ( sort { $barcode_distr_table->{$b}{num_seqs} <=> $barcode_distr_table->{$a}{num_seqs} } keys %$barcode_distr_table) {
		my $num_seqs = $barcode_distr_table->{$barcode_id}{num_seqs};
		my $num_bp = $barcode_distr_table->{$barcode_id}{num_bp};
		my $perc_seqs = ($total_seqs) ? sprintf("%.1f", (($num_seqs/$total_seqs)*100)) : "0.0";
		my $perc_bp = ($total_bp) ? sprintf("%.1f", (($num_bp/$total_bp)*100)) : "0.0";
		print join(" ", "# barcode", $barcode_id, $num_seqs, $perc_seqs, $num_bp, $perc_bp), "\n";
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
		my ($hit1, $hit2) = @{ $multibarcode_table->{$seq_id} };
		my @hits = @{ $multibarcode_table->{$seq_id} };
		my $report = "$seq_id\t";
		foreach my $H (@hits) {
			my ($pattern, $mismatches) = ($H->pattern, $H->num_mismatches);
			$report .= "\t$pattern ($mismatches)";
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
	my ($self, $Hits) = @_;
	
	# Sort the hits to a given sequence by number of mismatches
	my @Hits_sorted = sort { $a->num_mismatches <=> $b->num_mismatches } @$Hits;
	
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
	
	# Build our hash of sequences with their barcode assignments
	my $assignments_table;
	
	# Create our hash of sequences with multiple barcodes
	my $multibarcode_table;
	
	# Iterate through our hits, which we assume are sorted by sequence
	my ($curr_id, @Hits, $num_seqs_with_hits);
	while (my $Hit = $Hit_Iterator->()) {
		
		my $seq_id = $Hit->seq_id;
		
		if ($curr_id) {
			
			# Build a list of hits for this sequence
			if ($curr_id eq $seq_id) {
				push @Hits, $Hit;
				
			# Reached a new sequence
			} else {
				$num_seqs_with_hits++;
				my $num_hits = scalar @Hits;
				
				# Handle the previous set of Hits
				if ($num_hits) {
					
					# Assign sequence to its unique barcode hit
					if ($num_hits == 1) {
						$assignments_table->{$Hit->pattern}{$seq_id} = [ $Hit ];
						
					# Otherwise, check if this sequence has multiple barcode hits
					} else {
						my ($is_multicoded, $Hits_Sorted) = $self->check_multicoded_hits(\@Hits);
						if ($is_multicoded) {
							$multibarcode_table->{$seq_id} = $Hits_Sorted;
						
						# Assign the sorted hits
						} else {
							$assignments_table->{$Hits[0]->pattern}{$curr_id} = $Hits_Sorted;
						}
					}
				}
				
				# Start a list of hits for this new sequence
				@Hits = ($Hit);
				$curr_id = $seq_id;
			}
		
		# Start a list of hits for this new sequence
		} else {
			$num_seqs_with_hits++;
			@Hits = ($Hit);
			$curr_id = $seq_id;
		}
	}
	
	# Report count of multicoded sequences
	if ($self->verbose) {
		my $num_multibarcode_seqs = scalar keys %$multibarcode_table;
		print STDERR "Sequences with multi barcode hits: $num_multibarcode_seqs\n";
	}
	
	# Keep track of number of sequences with hits and assignments
	$self->num_seqs_with_hits($num_seqs_with_hits);
	
	# Set our multibarcode table
	$self->multibarcode_table($multibarcode_table);
	
	# Set and return our assignments table
	return $self->assignments_table($assignments_table);
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
	my ($outdir, $outformat) = ($self->outdir, $self->outformat);
	
	# Trim our sequences, and write the trim logfile
	print STDERR "Writing barcode reports\n";
	
	# Record the distribution of barcodes across our sequences (num_reads, num_bp per barcode)
	my $barcode_distr_table = $self->barcode_distr_table;
	
	foreach my $barcode_id (keys %$assignments_table) {
		
		# Put all of our trim files in this directory
		my $barcodedir = "$outdir/$barcode_id";
		mkdir $barcodedir unless -e $barcodedir;
		
		# Create a fasta file of barcodes
		# <barcode id>_<barcode sequence>_<num allowed mismatches>_<blinded_id>.fna
		my $Barcode = $barcode_table->{$barcode_id};
		my $barcode_fasta_file = "$barcodedir/$barcode_id.fasta";
		
		# Create a trim report file
		my $trim_file = "$barcodedir/$barcode_id.trim";
		open(FH_TRIMPOINTS, "> $trim_file") || die "Could not open file $trim_file for writing.\n";
		
		# write fasta output directly, or use Bio::SeqIO for writing fastq output
		my $seqio_fastq;
		unless ($self->trim_points_only) {
			print STDERR "Writing barcode results $barcode_fasta_file\n";
			if ($outformat eq "fastq") {
				$seqio_fastq = new Bio::SeqIO( -format => "fastq", -file => ">$barcode_fasta_file" );
			} else {
				open(FH_BARCODE_FILE, "> $barcode_fasta_file") || die "Could not open results file $barcode_fasta_file for writing.\n";
			}
		}
		
		# Get the sequences assigned to this barcode
		my @seq_ids = keys %{ $assignments_table->{$barcode_id} };
		
		foreach my $seq_id (@seq_ids) {
			
			# Get the hits for this barcode and sequence pair
			my $Hits = $assignments_table->{$barcode_id}{$seq_id};
			
			# Get the read sequence from our fasta lookup table
			die "Bug: No sequence in fasta table for id: $seq_id\n" unless defined $fasta_table->{$seq_id};
			my $F = $fasta_table->{$seq_id};
			my $seq = $F->seq; # untrimmed sequence
	
			# Get the clear range of our sequence (extracts barcode and clamp)
			my ($clear_start, $clear_end) = Grid::Tools::Barcode::Trimmer->trim_clear_range($seq, $Hits, $self->clamplength);
			
			# Throw out this read if we can't defined where the clear range start begins
			next unless defined $clear_start;
			
			# Throw out the read if it does not fulfill our minimum length requirements
			my $clear_seqlen = $clear_end - $clear_start;
			next if $self->min_readlength && $clear_seqlen < $self->min_readlength;
			
			# Report the distribution of barcodes across sequences
			$barcode_distr_table->{$barcode_id}{num_seqs}++;
			$barcode_distr_table->{$barcode_id}{num_bp} += $clear_seqlen;
			$barcode_distr_table->{total}{seqs}++;
			$barcode_distr_table->{total}{bp} += $clear_seqlen;
			
			# Include 454 key offset in writing the trim points
			# Write trimpoints in residue coordinates (add 1 to start position)
			print FH_TRIMPOINTS join("\t", $seq_id, $clear_start + $self->keylength + 1, $clear_end + $self->keylength), "\n";
			next if $self->trim_points_only;
			
			# Update our fasta record with the "clear range" trimmed sequence
			$F->seq(substr($seq, $clear_start, $clear_end-$clear_start));
			
			# Update our quality values with the "clear range" trimmed values
			my @quality_arr = split / /, $F->qual; 
			$F->qual(join(" ", @quality_arr[$clear_start..$clear_end-1]));
			
			# Verify that we have sequence and quality values
			die "Bug: No sequence", join("\n", $seq_id, $seq, "trimmed_seq ($clear_start..$clear_end) ",$F->seq, $F->qual), "\n" if !$F->seq;
			die "Bug: No qualities", join(" ", $seq_id, "trimmed_qual ($clear_start..$clear_end)", "length=", scalar @quality_arr, $F->qual), "\n" if !$F->qual;
			
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
				print FH_BARCODE_FILE join(" ", ">".$F->id, $desc), "\n", $F->fasta_qual;
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
			my $status_code = system($sc_seqidfile);
			print STDERR "Error: Failed system call to write seqid file with status code $status_code\n"
				if $status_code;			
			
			# Write the trimmed sff files using our trim file and sff file (-t $barcode_fasta_file $sff_file)
			my $sc_sff = "sfffile -o $trim_sff_file -i $seqid_file -t $trim_file $sff_file";
			print STDERR "Writing trimmed sff...\n$sc_sff\n" if $self->verbose;
			$status_code = system($sc_sff);
			print STDERR "Error: sfffile returns with status code $status_code while generating trimmed sff for barcode $barcode_id\n$sc_sff\n"
				if $status_code;
		}

	}

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
				errdir	=> $errdir,
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

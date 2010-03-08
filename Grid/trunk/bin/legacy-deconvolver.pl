#!/usr/local/bin/perl -w
use strict;
use Pod::Usage;
use Getopt::Long;
use Cwd;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../../";
use FileIO::BarcodeUtils;
use FileIO::FastaUtils;
use FileIO::FuzznucUtils;
use FileIO::Fasta;
use Grid::Tools::Barcode::Trimmer;

# Use DB_File DBM module and Storable for serializing and deserializing the data
use DB_File;
use MLDBM qw( DB_File Storable );

# BioPerl for reading and writing fastq files
use Bio::SeqIO;
use Bio::Seq::Quality;

######################## GLOBALS ############################

# Default Settings
my $FUZZNUC_BIN	= "/usr/local/packages/EMBOSS-5.0.0/bin/fuzznuc";
my $NUM_MISMATCHES 	= 2;
my $CWD 			= cwd(); 					# current working directory
my $TMP_DIR 		= "$CWD/tmp"; 				# temp directory
my $TRIM_READS 	= 1; 					# trim the barcodes off of the read sequences
my $TRIM_POINTS_ONLY = 0;					# output only the trim points
my $MIN_READ_LENGTH = 50; 					# Minimum acceptable read length after barcode trimming
my $CLAMP_LENGTH 	= 6; 					# Length of barcode clamp to trim off
my $KEY_LENGTH 	= 0; 					# Length of 454 keys, used for defining the trim points file.
my ($USE_DB, $REUSE_DB) = (0, 0);				# Use and reuse a DB File

#############################################################

# Get command-line options
my %options = ();
my $results = GetOptions( \%options, 
                          'barcode|b=s',
                          'sff|s=s',
                          'fasta|f=s',
                          'quals|q=s',
                          'fastq|fq=s',
                          'outformat|out=s',
                          'mismatches|m:i',
                          'readlength|r:i',
                          'clamplength|c:i',
                          'keylength|k:i',
                          'outdir|o:s',
                          'tmpdir|t:s',
					 'trim!',
					 'trimpoints!',
					 'usedb!',
					 'reusedb!',
                          'help|h'
) || &_pod;

&_pod if $options{help};

# Override default setting on optional parameters
my $OUTPUT_FORMAT = lc($options{outformat}) if defined $options{outformat};
my $FASTA_FORMAT = "fasta"; # assume fasta unless fastq is provided
$NUM_MISMATCHES = $options{mismatches} if defined $options{mismatches};
$TRIM_READS = 0 if defined $options{trim} && !$options{trim}; # true by default
$TRIM_POINTS_ONLY = 1 if defined $options{trimpoints} && $options{trimpoints}; # false by default
$CLAMP_LENGTH = $options{clamplength} if defined $options{clamplength};
$KEY_LENGTH = int($options{keylength}) if defined $options{keylength};
$MIN_READ_LENGTH = $options{readlength} if defined $options{readlength} && $options{readlength} > 10;
$TMP_DIR = $options{tmpdir} if defined $options{tmpdir};
$USE_DB = 1 if defined $options{usedb} && $options{usedb};
$REUSE_DB = 1 if defined $options{reusedb} && $options{reusedb};

# Verify input files are readable, and output directory is writable
my $BARCODE_FILE = $options{barcode} || &_pod("Error: A barcode metadata file is required.\n");
my $FASTA_FILE = $options{fasta} || "";
my $FASTQ_FILE = $options{fastq} || "";
my $FASTA_QUALS_FILE = $options{quals} || "";;
my $SFF_FILE = $options{sff} || "";
my ($SEQUENCE_FILE, $FORMAT) = ($FASTA_FILE, "fasta"); # assume fasta
if ($FASTQ_FILE) {
	$SEQUENCE_FILE = $FASTQ_FILE;
	$FORMAT = "fastq";
} elsif ($SFF_FILE) {
	$SEQUENCE_FILE = $SFF_FILE;
	$FORMAT = "sff";
}
my ($BASE_FASTA_FILE) = fileparse($FASTA_FILE, qr/\.[^.]*/);
my ($BASE_FASTQ_FILE)= fileparse($FASTQ_FILE, qr/\.[^.]*/);
my ($BASE_SFF_FILE)= fileparse($SFF_FILE, qr/\.[^.]*/);
my ($BASE_SEQUENCE_FILE) = fileparse($SEQUENCE_FILE, qr/\.[^.]*/);

# Ouput the results in the same format as input, unless user requests otherwise
$OUTPUT_FORMAT = $FORMAT unless $OUTPUT_FORMAT;

# Setup the temp and results directories unless they already exist
my $OUTPUT_DIR = (defined $options{outdir}) ? $options{outdir} : "$CWD/deconvolved_".lc($FORMAT);
mkdir $TMP_DIR unless -e $TMP_DIR;
mkdir $OUTPUT_DIR unless -e $OUTPUT_DIR;
my $MULTIBARCODE_FILE = "$OUTPUT_DIR/report_multibarcode.log";

# Veryify the input files are readable and output dirs are writeable
&_pod("Error: Unable to read input $FORMAT file $SEQUENCE_FILE\n") unless -r $SEQUENCE_FILE;
&_pod("Error: A read FASTA, FASTQ, or SFF file is required.\n") unless $FASTA_FILE || $FASTQ_FILE || $SFF_FILE;
&_pod("Error: Unable to read input barcode file $BARCODE_FILE\n") unless -r $BARCODE_FILE;
&_pod("Error: Unable to write to output directory $OUTPUT_DIR\n") unless -w $OUTPUT_DIR;
&_pod("Error: Unable to write to temp directory $TMP_DIR\n") unless -w $TMP_DIR;
&_pod("Error: Invalid output format $OUTPUT_FORMAT.  Choose fasta or fastq as the output format.\n") 
		unless $OUTPUT_FORMAT =~ /fasta|fastq|sff/;

# Log runtime settings
print join("\n",
		"########################################",
		"# Program: barcode_deconvolver.pl",
		"# RunTime: ". scalar localtime(),
		"# RunParameters: ",
		"#\t-fastq $FASTQ_FILE",
		"#\t-fasta $FASTA_FILE",
		"#\t-quals $FASTA_QUALS_FILE",
		"#\t-sff $SFF_FILE",
		"#\t-barcode $BARCODE_FILE",
		"#\t-mismatches $NUM_MISMATCHES",
		"#\t-readlength $MIN_READ_LENGTH",
		"#\t-clamplength $CLAMP_LENGTH",
		"#\t-keylength $KEY_LENGTH",
		"#\t-trim $TRIM_READS",
		"#\t-trimpoints $TRIM_POINTS_ONLY",
		"#\t-usedb $USE_DB",
		"#\t-reusedb $REUSE_DB",
		"#\t-outdir $OUTPUT_DIR",
		"#\t-tmpdir $TMP_DIR",
		"########################################\n"
);

# Get our barcodes
print STDERR "Reading barcodes $BARCODE_FILE\n";
my $barcode_table = FileIO::BarcodeUtils->parse_barcode_by_file($BARCODE_FILE);

if ($FORMAT =~ /fasta|fastq/) {
	# Clean and save input fasta or fastq sequence files to our tmp directory
	if ($FORMAT eq "fasta") {
	
		# write a clean fasta file
		$FASTA_FILE = "$TMP_DIR/$BASE_SEQUENCE_FILE.fasta";
		my $sc = "sed -s 's/\:/\_/g' $SEQUENCE_FILE > $TMP_DIR/$BASE_SEQUENCE_FILE.$FORMAT";
		&_pod("Error: Problem with cleaning input sequence file $SEQUENCE_FILE\n") 
				if system($sc);
		
	} elsif ($FORMAT eq "fastq") {
	
		# write a clean fastq file
		$FASTQ_FILE = "$TMP_DIR/$BASE_SEQUENCE_FILE.fastq";
		my $sc = "sed -s 's/\:/\_/g' $SEQUENCE_FILE > $FASTQ_FILE";
		&_pod("Error: Problem with cleaning input sequence file $SEQUENCE_FILE\n") 
				if system($sc);
		
		# Write FASTA sequence and qualities files from our input FASTQ file, required for fuzznuc searches
		$FASTA_FILE = "$TMP_DIR/$BASE_SEQUENCE_FILE.fasta";
	
		$FASTA_QUALS_FILE = "$TMP_DIR/$BASE_SEQUENCE_FILE.quals";
		print STDERR "Converting FASTQ to FASTA... $FASTA_FILE\n";
		my $num_seqs = FileIO::FastaUtils->write_fasta_from_fastq($FASTQ_FILE, $FASTA_FILE, $FASTA_QUALS_FILE);
		&_pod("Error: No sequences parsed from input fastq file $FASTQ_FILE\n") unless $num_seqs;
	}
} 

# Write a FASTQ file from our input SFF file
elsif ($FORMAT eq "sff") {
	$FASTA_FILE 		= "$TMP_DIR/$BASE_SEQUENCE_FILE.fasta";
	$FASTA_QUALS_FILE 	= "$TMP_DIR/$BASE_SEQUENCE_FILE.quals";
	
	# skip conversions if files exists
	unless (-f $FASTA_FILE && -f $FASTA_QUALS_FILE) {
		
		# (1) convert to fasta sequence file
		my $sc_fasta = "sffinfo -s $SEQUENCE_FILE > $FASTA_FILE";
		&_pod("Error: Problem with writing fasta from sff file $SEQUENCE_FILE\n") 
				if system($sc_fasta);
		
		# (2) convert to fasta qualities file
		my $sc_fasta_quals = "sffinfo -q $SEQUENCE_FILE > $FASTA_QUALS_FILE";
		&_pod("Error: Problem with writing fasta qualies from sff file $SEQUENCE_FILE\n") 
				if system($sc_fasta_quals);
	}
	# Note: process for sff file is now the same as input fasta + quals file
	$FORMAT = "fasta";
}

# $FORMAT is either "fastq" or "fasta"
# Guaranteed to have the fasta sequence and qualities files (either as input, or fastq->fasta, or sff->fasta)
# update path of our sequence file to our tmp directory
$SEQUENCE_FILE = "$TMP_DIR/$BASE_SEQUENCE_FILE.$FORMAT"; 

# Stores the best hits by barcode
my $results_table;

# Tracks reads with hits to multiple barcodes
my $multibarcode_table;

# Stores the fuzznuc hits by sequence
my $hit_table;

# Run fuzznuc against our barcode sequences
foreach my $Barcode (values %$barcode_table) {
	my ($barcode_id, $barcode_seq) = ($Barcode->id, $Barcode->seq);
	
	# Set the barcode clamp length based on user-defined parameter
	$Barcode->clamp_length($CLAMP_LENGTH);
	
	print STDERR "Barcode $barcode_id $barcode_seq\n";
	
	my $fuzznuc_file = "$TMP_DIR/fuzznuc_$barcode_id.csv";
	
	# system call to run fuzznuc
	my $exec_fuzznuc = join " ", (
			$FUZZNUC_BIN, 
			"-stdout true",
			"-filter", 
			"-complement Yes",
			"-pmismatch $NUM_MISMATCHES",
			"-sequence $FASTA_FILE",
			"-pattern $barcode_seq",
			
			# tab-delimited format
			"-rformat excel",
			
			# Skip the redundant header lines
			"| grep -v '^SeqName'", 
			
			# Reverse sort by number of mismatches (not necessary)
			# "| sort -rnk7",
			
			# output to our tmp directory
			"> $fuzznuc_file"
	);
	
	# Check if fuzznuc was already run
	unless (-r $fuzznuc_file) {
		print STDERR "\tRunning fuzznuc...";
		
		# Run fuzznuc and check the status code that is returned
		my $status_code = system($exec_fuzznuc);
		
		# fuzznuc should return a status code of 0, if it completes successfully
		if ($status_code) {
			print STDERR "\nError: fuzznuc returns with a status code $status_code\n";
		} else {
			print STDERR "ok\n";
		}
	}
	
	# Parse fuzznuc output
	print STDERR "\tParsing fuzznuc...";
	my $Hits = FileIO::FuzznucUtils->parse_fuzznuc_by_csv_file($fuzznuc_file);
	print STDERR "ok (".scalar @$Hits." hits)\n";
	
	# Add hits to our hit table
	foreach my $H (@$Hits) {
		$H->{barcode} = $Barcode;
		push @{ $hit_table->{$H->seq_id} }, $H;
	}
}

# Get a lookup table from our fasta file
print STDERR "Reading $FORMAT file: $SEQUENCE_FILE\n";

# Store fasta sequences in Berkeley DB_File
my $fasta_table;
if ($USE_DB || $REUSE_DB) {
	my $fasta_db = "$TMP_DIR/_fasta.db";
	unlink $fasta_db unless $REUSE_DB;
	tie %$fasta_table, 'MLDBM', $fasta_db, O_RDWR|O_CREAT, 0666 
		|| die "Cannot open file '$fasta_db': $!\n";
}

# Parse fasta (faster than BioPerl's code for parsing fastq)
my $num_seqs = scalar keys %$fasta_table; print STDERR "fasta.db: $num_seqs sequences\n";
unless ($num_seqs) {
	FileIO::FastaUtils->parse_fasta_by_file($FASTA_FILE, $fasta_table);
	FileIO::FastaUtils->parse_quals_by_file($FASTA_QUALS_FILE, $fasta_table)
		unless $TRIM_POINTS_ONLY;
	$num_seqs = scalar keys %$fasta_table;
}

# Stats for our log file
my $num_barcodes = scalar keys %$barcode_table;
my $num_seqs_with_hits = scalar keys %$hit_table;
my $perc_seqs_with_hits = ($num_seqs) ? sprintf("%.1f", (($num_seqs_with_hits/$num_seqs)*100)) : "0.0";
print join("\n",
	"# Number of barcode entries: $num_barcodes",
	"# Number of fasta entries: $num_seqs",
	"# Number of sequences with barcode hits: $num_seqs_with_hits",
	"# Percent of sequences with barcode hits: $perc_seqs_with_hits",
), "\n";

# Assign sequences to barcodes
print STDERR "Assigning sequences to barcodes...\n";
foreach my $seq_id (keys %$hit_table) {
	
	# Get the list of fuzznuc hits for this sequence
	my $Hits = $hit_table->{$seq_id};
	
	# Sort the hits to a given sequence by number of mismatches
	my @Hits_sorted = sort { $a->num_mismatches <=> $b->num_mismatches } @$Hits;
	
	# Get all of the hits to the barcode with the fewest mismatches
	my (@BestHits, $is_ambiguous);
	for (my $i=0; $i<scalar @Hits_sorted; $i++) {
		my $Hit = $Hits_sorted[$i];
		
		# Set the first hit as our best hit
		if (!scalar @BestHits) {
			push @BestHits, $Hit;
			
		# Then, confirm that its our unambiguous best hit
		} else {
			my $BestHit = $BestHits[0];
			
			# Add any hits to this same barcode to our list
			if ($BestHit->barcode->id eq $Hit->barcode->id) {
				push @BestHits, $Hit;
			
			# Throw out any sequences with hits to multiple barcodes within $NUM_MISMATCHES
			} else {
				$multibarcode_table->{$seq_id} =[$BestHit, $Hit];
				last;
			}
		}
	}
	# Assign the list of best hits to this barcode-sequence pair
	unless (defined $multibarcode_table->{$seq_id}) {
		# Avoid redundant best hits
		foreach my $Hit (@BestHits) {
			$results_table->{$Hit->barcode->id}{$seq_id} = \@BestHits;
		}
	}
}

print STDERR "Number of sequences with ambiguous best barcode hits... ", scalar keys %$multibarcode_table, "\n";

# Trim our sequences, and write the trim logfile
print STDERR "Writing barcode reports\n";

# Record the distribution of barcodes across our sequences (num_reads, num_bp per barcode)
my $barcode_distr_table;

foreach my $barcode_id (keys %$results_table) {
	
	# Create a fasta file of barcodes
	# <barcode id>_<barcode sequence>_<num allowed mismatches>_<blinded_id>.fna
	my $Barcode = $barcode_table->{$barcode_id};
	my $barcode_results_file = "$OUTPUT_DIR/".join("_", $barcode_id, $Barcode->seq, $Barcode->blinded_id).".$OUTPUT_FORMAT";
	print STDERR "Writing barcode results $barcode_results_file\n";
	
	# Create a trim report file
	my $trim_file = "$OUTPUT_DIR/trim_$barcode_id.txt";
	open(FH_TRIMPOINTS, "> $trim_file") || die "Could not open file $trim_file for writing.\n";
	
	# write fasta output directly, or use Bio::SeqIO for writing fastq output
	my $seqio_fastq;
	unless ($TRIM_POINTS_ONLY) {
		if ($OUTPUT_FORMAT eq "fasta") {
			open(FH_BARCODE_FILE, "> $barcode_results_file") || die "Could not open results file $barcode_results_file for writing.\n";
		} elsif ($OUTPUT_FORMAT eq "fastq") {
			$seqio_fastq = new Bio::SeqIO( -format => "fastq", -file => ">$barcode_results_file" );
		} else {
			open(FH_BARCODE_FILE, "> $barcode_results_file") || die "Could not open results file $barcode_results_file for writing.\n";
		}
	}
	
	# Get the sequences assigned to this barcode
	my @seq_ids = keys %{ $results_table->{$barcode_id} };
	
	foreach my $seq_id (@seq_ids) {
		
		# Get the hits for this barcode and sequence pair
		my $Hits = $results_table->{$barcode_id}{$seq_id};
		
		# Get the read sequence from our fasta lookup table
		die "Bug: No sequence in fasta table for id: $seq_id\n" unless defined $fasta_table->{$seq_id};
		my $F = $fasta_table->{$seq_id};
		my $seq = $F->seq; # untrimmed sequence

		# Get the clear range of our sequence, after trimming off barcodes
		my ($clear_start, $clear_end) = Barcode::Trimmer->trim_clear_range($seq, $Hits, $CLAMP_LENGTH);
		
		# Throw out this read if we can't defined where the clear range start begins
		next unless defined $clear_start;
		
		# Throw out the read if it does not fulfill our minimum length requirements
		my $clear_seqlen = $clear_end - $clear_start;
		next unless $clear_seqlen >= $MIN_READ_LENGTH;
		
		# Report the distribution of barcodes across sequences
		$barcode_distr_table->{$barcode_id}{num_seqs}++;
		$barcode_distr_table->{$barcode_id}{num_bp} += $clear_seqlen;
		$barcode_distr_table->{total}{seqs}++;
		$barcode_distr_table->{total}{bp} += $clear_seqlen;
		
		# Include 454 key offset in writing the trim points
		# Write trimpoints in residue coordinates (add 1 to start position)
		print FH_TRIMPOINTS join("\t", $barcode_id, $seq_id, $clear_start + $KEY_LENGTH + 1, $clear_end + $KEY_LENGTH), "\n";
		next if $TRIM_POINTS_ONLY;
		
		# Update our fasta record with the "clear range" trimmed sequence
		$F->seq(substr($seq, $clear_start, $clear_end-$clear_start));
		
		# Update our quality values with the "clear range" trimmed values
		my @quality_arr = split / /, $F->qual; 
		$F->qual(join(" ", @quality_arr[$clear_start..$clear_end-1]));
		
		die "Bug: No sequence", join("\n", $seq_id, $seq, "trimmed_seq ($clear_start..$clear_end) ",$F->seq, $F->qual), "\n" if !$F->seq;
		die "Bug: No qualities", join(" ", $seq_id, "trimmed_qual ($clear_start..$clear_end)", "length=", scalar @quality_arr, $F->qual), "\n" if !$F->qual;
		
		# Report the locations of our barcode hits
		my $locs_string;
		my @sorted_hits = sort { $a->{min} <=> $b->{min} } @$Hits;
		foreach my $Hit (@sorted_hits) {
			$locs_string .= $Hit->min."..".$Hit->max.":".$Hit->strand.",";
		}
		$locs_string =~  s/\,+$//; # Remove trailing comma
		
		# Throw out any reads that do not meet the minimum length requirements
		my $seqlen = length($F->seq);
		
		# After trimming is complete, we can write out the barcode-trimmed sequence report
		if ($OUTPUT_FORMAT eq "fasta") {
			my $desc = join(" ", "barcode=$barcode_id", "length=$seqlen", "barcode_loc=$locs_string", $F->desc);
			print FH_BARCODE_FILE join(" ", ">".$F->id, $desc), "\n", $F->fasta_seq;
			print FH_BARCODE_FILE join(" ", ">".$F->id, $desc), "\n", $F->fasta_qual;
		
		# Use Bio::Seq::Quality to write out fastq
		} elsif ($OUTPUT_FORMAT eq "fastq") {
			my $SeqWithQuality = new Bio::Seq::Quality( 
									-id => $F->id,
									-seq => $F->seq,
									-qual => $F->qual,
									-desc => ""
			);
			$seqio_fastq->write_fastq($SeqWithQuality);
			
		} elsif ($OUTPUT_FORMAT eq "sff") {
			# Write the trimmed sff files
			my $sc_sff = "sfffile -o $$OUTPUT_DIR/trimmed_sff/$barcode_id.sff -i trimfiles/trim_$barcode_id -t trimfiles/trim_$barcode_id $SFF_FILE";
			my $status_code = system($sc_sff);
			print STDERR "Error: sfffile returns with status code $status_code while generating trimmed sff for barcode $barcode_id\n$sc_sff\n"
				if $status_code;
		}
	}
	close FH_BARCODE_FILE;
}

# Write the trimmed sff files based on the trim files
# for bc in `ls trimfiles | sed -s 's/trim_//g'`; do sfffile -o trimmed_sff/$bc.sff -i trimfiles/trim_$bc -t trimfiles/trim_$bc ../input_data/sff/input1.sff; done;

# Log any ambiguities
print STDERR "Writing multi-barcode reports\n";
print "# Number of multi-barcoded sequences: ", scalar keys %$multibarcode_table, "\n";
open(FH_MULTIBARCODE, ">$MULTIBARCODE_FILE") || die "Could not open trim file $MULTIBARCODE_FILE for writing.\n";
foreach my $seq_id (keys %$multibarcode_table) {
	my ($hit1, $hit2) = @{ $multibarcode_table->{$seq_id} };
	print FH_MULTIBARCODE join("\t", $seq_id, $hit1->barcode->id, $hit1->num_mismatches, $hit2->barcode->id, $hit2->num_mismatches), "\n";
}
close FH_MULTIBARCODE;

# Log our barcode distribution
my ($total_seqs, $total_bp) = ($barcode_distr_table->{total}{seqs}, $barcode_distr_table->{total}{bp});
delete $barcode_distr_table->{total};

# Sort our barcodes by the number of sequences in descending order
print "# Distribution of barcodes: <barcode_id> <num_sequences> <perc_sequences> <num_bp> <percent_bp>\n";
for my $barcode_id ( sort { $barcode_distr_table->{$b}{num_seqs} <=> $barcode_distr_table->{$a}{num_seqs} } keys %$barcode_distr_table) {
	my $num_seqs = $barcode_distr_table->{$barcode_id}{num_seqs};
	my $num_bp = $barcode_distr_table->{$barcode_id}{num_bp};
	my $perc_seqs = ($total_seqs) ? sprintf("%.1f", (($num_seqs/$total_seqs)*100)) : "0.0";
	my $perc_bp = ($total_bp) ? sprintf("%.1f", (($num_bp/$total_bp)*100)) : "0.0";
	print join(" ", "# barcode", $barcode_id, $num_seqs, $perc_seqs, $num_bp, $perc_bp), "\n";
}

######################## SUB ROUTINES ############################

sub _pod {
	print STDERR "\n", shift, "\n" if @_; # display a specific error message if provided
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR } );
}

=head1 NAME
	
barcode_deconvolver.pl - Deconvolves a fasta file by barcode, and optionally trims off 
					barcode sequences by default.
	
=head1 SYNOPSIS

USAGE:
	
	barcode_deconvolver.pl 
		--barcode /path/to/barcodes.txt 
	[	--fastq /path/to/read/fastq
	     --fasta /path/to/read.fasta 
	  	--sff /path/to/read.sff 
	  	--outformat fastq
		--mismatches 2
		--readlength 50
		--clamplength 6
		--keylength 0
		--trimpoints
		--usedb
		--reusedb
		--outdir /path/to/output_dir
		--tmpdir /path/to/tmp_dir
		--help 
	 ]


=head1 OPTIONS

B<--barcode,-b>
	REQUIRED. Path to barcode metadata file.

B<--fasta,-f>
	OPTIONAL.  Path to read sequence FASTA file.
	A fasta or fastq file is required.
	
B<--fastq,-fq>
	OPTIONAL.  Path to read sequence FASTQ file.
	The fastq parser use BioPerl and its slow!

B<--output,-o>
	OPTIONAL.  Path to output directory.  Directory will be created if not found.

B<--mismatches,-m>
	OPTIONAL.  Number of mismatches to allow in running fuzznuc searches.

B<--trim,-t>
	OPTIONAL.  Trim the barcode sequences off of each read.

B<--readlength,-l>
	OPTIONAL.  Minimal acceptable read length after barcode trimming is complete.

B<--clamplength,-l>
	OPTIONAL.  The length of the barcode clamp. A/K/A hexamer length.  Default is 6.

B<--keylength,-l>
	OPTIONAL.  The trim points for 454 sequences typically need to be offset by
			 the key length of 4bp.  This is used for defining the trim points
			 file.  Default is 0.

B<--usedb>
	OPTIONAL.  Use a Berkeley DB File to store datasets. Default is false.
		
B<--reusedb>
	OPTIONAL.  Reuse the DB File from a previous run <tmp_dir/fasta.db>..  
			 Default is false.  This is useful for rapid testing purposes.

B<--help,-h>
	Print this message

=head1  DESCRIPTION

This script deconvolves a fasta, fastq or sff file based on the barcode sequences, as 
determined by running fuzznuc to find the best hits of read to barcode sequences.
The results are a report of trim points, and a list of barcode fasta files with 
each entry representing a read sequence that has its unambiguous best hit to 
the bar code.  The bar code sequences are trimmed by default, unless otherwise 
specified.  

All results are in space-based coordinates.

The rules for trimming reads:
1) If barcode hits overlap, extend the hit and trim the extended region;
2) Allow any number of hits, as long as the trimmed read length >= $MIN_READ_LENGTH
3) Throw out and log any sequences with hits to multiple barcodes


=head1  CONTACT
	
	Nelson Axelrod
	naxelrod@jcvi.org

=head1  CAVEATS

	There's a bug in Fuzznuc where it reports the ids incorrectly.  
	Eg. seq_id SOLEXA1:1:100:1016:995#0/1 is reported as 995#0/1
	Requires sed -s 's/\:/\_/g' on fasta input files
	
=cut



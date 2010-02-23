package FileIO::FastaUtils;
use strict;
require FileIO::Fasta;
use Bio::SeqIO;

=item $obj->parse_fastq_by_file();

B<Description:> This utility method is used to parse a FASTQ file, and
return a lookup table of FileIO::Fasta objects, where each object represents a 
fasta sequence entry including sequence and quality information.

B<Parameters:> Accepts a path to a fasta file.  Accepts a fasta_table
to make the function reentrant.

B<Returns:> Hashref with keys as sequence ids and values of FileIO::Fasta objects

=cut

sub parse_fastq_by_file {
    my ($self, $fastq_file, $fasta_table) = @_;
    
    my $seq_in = Bio::SeqIO->new( -file => $fastq_file, -format => "fastq" );
    while (my $inseq = $seq_in->next_seq) {
    	    
    	    # Add this FileIO::Fasta object to our hash table
    	    my $id = $inseq->id;
    	    die "Error: sequence $id contains colon character in fastq file $fastq_file\n"
    	    		if $id =~ /:/;
    	    		
    	    $fasta_table->{$id} = new FileIO::Fasta({
    	    		    id => $id,
    	    		    desc => "", 
    	    		    seq => $inseq->seq, 			# sequence string
    	    		    qual => $inseq->qual_text	# string of quality values
    	    });
    	    #print STDERR "parse: $i\n" if !($i++ % 100000);
    }
    return $fasta_table;
}

=item $obj->parse_fasta_by_file();

B<Description:> This utility method is used to parse a fasta file, and
return a lookup table of FileIO::Fasta objects, where each object represents a 
fasta sequence entry.

B<Parameters:> Accepts a path to a fasta file.  Accepts a fasta_table
to make the function reentrant.

B<Returns:> Hashref with keys as sequence ids and values of FileIO::Fasta objects

=cut

sub parse_fasta_by_file {
    my ($self, $fasta_file, $fasta_table) = @_;
    open(FH_FASTA, "< $fasta_file") || die "Could not open fasta file $fasta_file for reading.\n";
    
    # delimit sequences by header line
    local $/ = "\n>";
    
    while (<FH_FASTA>) {
    	    chomp;
    	    my ($header, @sequence_lines) = split /\n/, $_;
    	    $header = substr($header,1) if $header=~/^>/; # strip initial > char
    	    
    	    # Separate the identifier from the rest of the header line
    	    my ($id, $desc);
    	    if ($header =~ /^\s*(\S+)\s*(.*)/) {
    	    	    ($id, $desc) = ($1, $2);
    	    }
    	    
    	    # Get the sequence
    	    my $seq = join("", @sequence_lines);
    	    $seq =~ s/\s//g; # strip whitespace and newlines
    	    
	    if (exists $fasta_table->{$id}) {
		    # Get the object and set its quality values
		    my $F = $fasta_table->{$id};
		    $F->seq($seq);
	    } else {
		    # Add a new FileIO::Fasta object to our hash table
		    $fasta_table->{$id} = new FileIO::Fasta({
						id => $id,
						desc => $desc,
						seq => $seq,
		    });
	    }
    }
     
     close FH_FASTA;
     return $fasta_table;
}

=item $obj->parse_quals_by_file();

B<Description:> Parses a FASTA qualities file, and adds the quality 
information to a lookup table of FileIO::Fasta objects, where each object represents 
a fasta sequence and qualities entry.

B<Parameters:> Accepts a path to a fasta quals file.  Accepts a fasta_table
to make the function reentrant.

B<Returns:> Hashref with keys as sequence ids and values of FileIO::Fasta objects

=cut

sub parse_quals_by_file {
    my ($self, $fasta_file, $fasta_table) = @_;
    open(FH_QUALS, "< $fasta_file") || die "Could not open fasta file $fasta_file for reading.\n";
    
    # delimit sequences by header line
    local $/ = "\n>";
    
    while (<FH_QUALS>) {
    	    chomp;
    	    my ($header, @quals_lines) = split /\n/, $_;
    	    $header = substr($header,1) if $header=~/^>/; # strip initial > char
    	    
    	    # Separate the identifier from the rest of the header line
    	    my ($id, $desc);
    	    if ($header =~ /^\s*(\S+)\s*(.*)/) {
    	    	    ($id, $desc) = ($1, $2);
    	    }
    	    
    	    # Get the quality values
    	    my $quals = join(" ", @quals_lines);
    	    $quals =~ s/[\n\s]+/ /g; # 1 space between values
    	    
	    if (exists $fasta_table->{$id}) {
		    my $F = $fasta_table->{$id};
		    $F->qual($quals);
		    
	    } else {
		    # Add a new FileIO::Fasta object to our hash table
		    $fasta_table->{$id} = new FileIO::Fasta({
						id => $id,
						desc => $desc,
						qual => $quals,
		    });
	    }
     }
     
     close FH_QUALS;
     return $fasta_table;
}

=item $obj->count_special_str();

B<Description:> Counts the number of lines that contain a given string.  
This is used during barcode deconvolution to identify fasta files with 
invalid special characters.

B<Parameters:> Accepts a path to a fasta file, and the special string.

B<Returns:> Number of lines that contain the string.

=cut

sub count_special_str {
    my ($self, $fasta_file, $special_str) = @_;
    my $badline_count = `grep ':' $fasta_file | wc -l`; 
    $badline_count =~ /^(\d+).*/g;
    return $1;
}

=item $obj->write_fasta_from_fastq();

B<Description:> Takes a fastq and writes a fasta file.

B<Parameters:> Accepts path to a fastq file, and a path to write a fasta file

B<Returns:> Number of sequences written.

=cut


sub write_fasta_from_fastq {
	my ($self, $fastq_file, $fasta_file, $fasta_quals_file) = @_;
	my $seq_in = Bio::SeqIO->new(-file => $fastq_file, -format => "fastq" );
	my $seq_out = Bio::SeqIO->new( -file => ">$fasta_file", -format => "fasta" );
	my $quals_out = Bio::SeqIO->new( -file => ">$fasta_quals_file", -format => "qual" );
	my $i=0;
	while (my $inseq = $seq_in->next_seq) {
		# write fasta sequence file $fasta_file
		$seq_out->write_seq($inseq); 
		
		# write fasta qualities file $fasta_quals_file
		# Note: quality values in BioPerl have the incorrect offset for quality values
		$quals_out->write_seq($inseq);
		
		$i++;
	}
	return $i;
}

=item $obj->write_fastq_from_fasta();

B<Description:> Takes a fasta and writes a fastq file.

B<Parameters:> Accepts path to a fasta file, and a path to write a fastq file

B<Returns:> Number of sequences written.

=cut


sub write_fastq_from_fasta {
	my ($self, $fastq_file, $fasta_file) = @_;
	my $seq_in = Bio::SeqIO->new(-file => $fasta_file, -format => "fasta" );
	my $seq_out = Bio::SeqIO->new( -file => ">$fastq_file", -format => "fastq" );
	my $i=0;
	while (my $inseq = $seq_in->next_seq) {
		$seq_out->write_seq($inseq);
		
		# Note: quality values in BioPerl have the incorrect offset for quality values
		$i++;
	}
	return $i;
}

1;

__END__

=head1 NAME

FileIO::FastaUtils - Provides utilities for reading fasta files

=head1 SYNOPSIS

  use FileIO::FastaUtils;
  my $fasta_table = FileIO::FastaUtils->parse_fasta_by_file($fasta_file);
  foreach my $seq_id (keys %$fasta_table) {
  	  my $F = $fasta_table->{$seq_id};
  	  print join(" ", ">", $F->id, $F->desc), "\n", $F->fasta_seq(), "\n";
  }

=head1 DESCRIPTION

This module provides some methods to read fasta files and return a lookup table
FileIO::Fasta objects, with sequence ids as keys.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut

#!/usr/local/bin/perl -w
use strict;
use FileIO::FastaUtils;
use Bio::SeqIO;
use Bio::Seq::Quality;

my $fasta_seqs_file = shift @ARGV || die "A fasta sequence is required.\n";
my $fasta_qual_file = shift @ARGV || die "A fasta qualities file is required.\n";

# Read fasta sequences into tmp_table
my $fasta_table;
$fasta_table = FileIO::FastaUtils->parse_fasta_by_file($fasta_seqs_file, $fasta_table);
$fasta_table = FileIO::FastaUtils->parse_quals_by_file($fasta_qual_file, $fasta_table);

# Define our fasta to fastq Bio::SeqIO converter
my $seqio_fastq = new Bio::SeqIO( -format => "fastq", -file => ">sequences.fastq" );

foreach my $F (values %$fasta_table) {
	# Use Bio::Seq::Quality to write out fastq
	my $SeqWithQuality = new Bio::Seq::Quality( 
							-id => $F->id,
							-seq => $F->seq,
							-qual => $F->qual,
							-desc => $F->desc
	);
	$seqio_fastq->write_fastq($SeqWithQuality);
}


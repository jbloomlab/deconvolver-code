#!/usr/local/bin/perl -w
use strict;
use IO::File;
use FileIO::FastaUtils;

my $fasta_file = "/usr/local/scratch/naxelrod/tmp/s_1_1_sequence.fasta";
my $fasta_table = {};

FileIO::FastaUtils->parse_fasta_by_file($fasta_file, $fasta_table);

my $id = 'SOLEXA1_1_36_1096_1964#0/1';
my $Fasta = $fasta_table->{$id};

print "found ", $fasta_table->{$id}->id, $fasta_table->{$id}->seq, "\n";



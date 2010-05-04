package FileIO::FuzznucUtils;
use strict;
use FileIO::FuzznucHit;
use IO::File;

=head2 hit_iterator()
	
	Returns an Iterator of FileIO::FuzznucHit objects,
	given a list of files from Fuzznuc csv output.
	
	Returns undef once we have finished iterating all files.
	
=cut

sub hit_iterator {
	my ($self, $files, $verbose) = @_;
	
	# Open the first file for reading
	my $file = shift @$files;
	my $fh = IO::File->new($file) || die "Could not open fuzznuc results file $file.\n";
	
	my ($Hit_Iterator, $i, $done);
	$Hit_Iterator = sub {
		return undef if $done;

		# Iterate thru our filehandle, and return the next Hit
		while (<$fh>) {
			chomp;
			next if $_ =~ /^SeqName/; # skip headers
			
			# my ($seq_id, $start, $end, $length, $strand, $pattern, $num_mismatches) = split /\t/, $_;
			my ($seq_id, $start, $end, $length, $strand, $pattern_line, $num_mismatches) = split /\t/, $_;
			my ($pattern) = split /\s+/, $pattern_line;
			
			$num_mismatches = 0 if $num_mismatches eq ".";
			# Use min, max, space-based coordinates
			my ($min, $max) = ($start < $end) ? ($start, $end) : ($end, $start);
			$min--;
			
			return new FileIO::FuzznucHit({ 
					seq_id => $seq_id, 
					min => $min, 
					max => $max,
					length => $length,
					strand => $strand,
					pattern => $pattern,
					num_mismatches => $num_mismatches
			});
		}
		close $fh;
		
		# Open the next file for reading
		$file = shift @$files;
		if ($file) {
			$fh = IO::File->new("< $file") || die "Could not open fuzznuc output file $file\n";
			return $Hit_Iterator->();
		} else {
			# Finished iterating thru all of the files
			$done = 1 unless scalar @$files;
			return undef;
		}
	};
	
	return $Hit_Iterator;
}

=head2 hit_pump_iterator()
	
	Returns an Iterator of FileIO::FuzznucHit objects,
	given a list of files from Fuzznuc csv output.
	
	Usage is like so
	while ($Iter->('has_next')) {
		my $hit = $Iter->(); # or $Iter->('next')
		next unless $hit;
		
	}
	
=cut

sub hit_pump_iterator {
	my ($self, $files, $num_expected_files, $verbose) = @_;
	
	# Open the first file for reading
	my $file = shift @$files;
	my $fh = IO::File->new($file) || die "Could not open fuzznuc results file $file.\n";
	
	my ($Hit_Iterator, $i, $done);
	$Hit_Iterator = sub {
		# Users can call $Iterator->('has_next') 
		# to check if we have exhausted all of the data files that fuel our 
		# iterator pump
		my $action = shift || 'next';
		
		return undef if $done;

		# Iterate thru our filehandle, and return the next Hit
		while (<$fh>) {
			chomp;
			next if $_ =~ /^SeqName/; # skip headers
			
			# my ($seq_id, $start, $end, $length, $strand, $pattern, $num_mismatches) = split /\t/, $_;
			my ($seq_id, $start, $end, $length, $strand, @barcode_desc) = split /\t/, $_;
			my $pattern = shift @barcode_desc;
			my $num_mismatches = pop @barcode_desc;
			$num_mismatches = 0 if $num_mismatches eq ".";
			# Use min, max, space-based coordinates
			my ($min, $max) = ($start < $end) ? ($start, $end) : ($end, $start);
			$min--;
			
			return new FileIO::FuzznucHit({ 
					seq_id => $seq_id, 
					min => $min, 
					max => $max,
					length => $length,
					strand => $strand,
					pattern => $pattern,
					num_mismatches => $num_mismatches
			});
		}
		close $fh;
		
		# Open the next file for reading
		$file = shift @$files;
		if ($file) {
			$fh = IO::File->new("< $file") || die "Could not open fuzznuc output file $file\n";
			return $Hit_Iterator->();
		
		} else {
			
			# Finished iterating thru all of the files
			$done = 1 unless scalar @$files;
			return undef;
		}
	};
	
	return $Hit_Iterator;
}


=item $obj->parse_fuzznuc_by_csv_file();

B<Description:> This utility method is used to parse a fuzznuc Excel (csv) file, 
and returns a list of Fuzznuc objects, where each object represents a fuzznuc entry.

B<Parameters:> Accepts a path to a fuzznuc csv file

B<Returns:> Listref of FuzznucHit objects

=cut

sub parse_fuzznuc_by_csv_file {
	my ($self, $file) = @_;
	my @Hits;
	
	open(FH_FUZZNUC, "< $file") || die "Could not open fuzznuc file $file for reading.\n";
	while (<FH_FUZZNUC>) {
		chomp;
		next if $_ =~ /^SeqName/; # skip headers
		my ($seq_id, $start, $end, $length, $strand, $pattern, $num_mismatches) = split /\t/, $_;
		$num_mismatches = 0 if $num_mismatches eq ".";
		# Use min, max, space-based coordinates
		my ($min, $max) = ($start < $end) ? ($start, $end) : ($end, $start);
		$min--;
		
		push @Hits, new FileIO::FuzznucHit({ 
				seq_id => $seq_id, 
				min => $min, 
				max => $max,
				length => $length,
				strand => $strand,
				pattern => $pattern,
				num_mismatches => $num_mismatches
		});
	}
	close FH_FUZZNUC;
	return \@Hits;
}

=item $obj->parse_fuzznuc_by_file();

B<Description:> This utility method is used to parse a fuzznuc file, and
returns a list of Fuzznuc objects, where each object represents a fuzznuc entry.

B<Parameters:> Accepts a path to a fuzznuc file

B<Returns:> Listref of FuzznucHit objects

=cut

sub parse_fuzznuc_by_file {
	my ($self, $file) = @_;
	my (@Hits, $seq_id);
	open(FH_FUZZNUC, "< $file") || die "Could not open fuzznuc file $file for reading.\n";
	while (<FH_FUZZNUC>) {
		chomp;
		if($_=~/^# Sequence: (\S+) /){
			$seq_id=$1;
		}elsif($_=~/^#/ || $_=~/^$/){
			# Noop
		}elsif($_=~/^  Start\s+End\s+Pattern_name\s+Mismatch\s+Sequence/){
			# Noop
		}else{
			$_=~s/^\s+//;
			
			my ($start, $end, $pattern, $num_mismatches) = split /\s+/, $_;
			$num_mismatches = 0 if $num_mismatches eq ".";
			# Use min, max coordinates
			my ($min, $max, $strand);
			if ($start < $end) {
				($min, $max) = ($start, $end);
				$strand = "+";
			} else {
				($min, $max) = ($end, $start);
				$strand = "-";
			}
			$min--; # space-based coordinates
			
			my $H = new FileIO::FuzznucHit({ 
					seq_id => $seq_id, 
					min => $min, 
					max => $max,
					length => $max - $min,
					strand => $strand,
					pattern => $pattern,
					num_mismatches => $num_mismatches
			});
			
			push @Hits, $H;
		}
	}
	close FH_FUZZNUC;
	return \@Hits;
}

=item $obj->filter_redundant_hits();

B<Description:> This method takes a list of FuzznucHit objects, and 
removes any redundant hits from the list.  It was developed because Fuzznuc
will report identical hits in its output.  Redundancy is defined as a given
hit from and to the same sequences at the same location.

B<Parameters:> Accepts a filehandle to a fuzznuc file

B<Returns:> Listref of FuzznucHit objects

=cut

sub filter_redundant_hits {
	my ($self, $Hits) = @_;
	
	# Make sure the hits are non-redundant
	my (@NonRedundantHits);
	
	my $hit_table;
	foreach my $H (@$Hits) {
		push @NonRedundantHits, $H unless defined $hit_table->{$H->seq_id}{$H->min}{$H->max};
		$hit_table->{$H->seq_id}{$H->min}{$H->max} = 1;
	}
	return \@NonRedundantHits;
}

=item $obj->merge_overlapping_hits();

B<Description:> This method takes a listref of FileIO::FuzznucHit objects,
merges any overlapping barcode hits and extends the hit region of each merged hit.
This is typically used prior to applying the barcode trimming of sequences.

B<Parameters:> Accepts a listref of FileIO::FuzznucHit objects

B<Returns:> a listref of merged FileIO::FuzznucHit objects sorted by
min position in ascending order.



=cut

sub merge_overlapping_hits {
	
	my ($self, $Hits) = @_;
	return $Hits unless scalar @$Hits >= 2;
	
	# Sort the hits by min position in descending order
	my @sorted_hits = sort { $a->min <=> $b->min } @$Hits;
	
	my $barcode_id = $sorted_hits[0]->barcode->id;
	
	# Merge our overlapping hits
	my @merged_hits;
	for (my $i=0; $i<scalar @sorted_hits - 1; $i++) {
		
		# Check for overlap with the next hit
		my ($H1, $H2) = ($sorted_hits[$i], $sorted_hits[$i+1]); 
		
		# This function should only be called on hits to the same barcode
		die "ERROR: Attempt to merge_overlapping_hits() on hits to multiple barcodes.\n"
			if $H2->barcode->id ne $barcode_id;
		
		# Merge the overlapping hits
		if ($H1->max >= $H2->min) {
			# Extend H2 to the start of H1
			$H2->min($H1->min);
			# and exclude H1 from our merged results
			push @merged_hits, $H2;
			
		# No merge, so add $H2 to the list
		} else {
			push @merged_hits, ($H1, $H2);
		}
	}
	
	return \@merged_hits;
}

=item $obj->extend_hits_with_clamp($clamp_length);

B<Description:> This method takes a listref of FileIO::FuzznucHit objects,
and extends the hit location on both 5' and 3' ends by the lenght of the barcode
clamp.

B<Parameters:> Accepts a listref of FileIO::FuzznucHit objects, table of fasta 
entries, and the clamp length

B<Returns:> a listref of merged FileIO::FuzznucHit objects extended in both directions
by the length of the clamp



=cut

sub extend_hits_with_clamp {
	my ($self, $Hits, $fasta_table, $clamp_length) = @_;
	my @Extended_Hits;
	foreach my $Hit (@$Hits) {
		my $seqlen = length($fasta_table->{$Hit->seq_id}->seq);
		my ($min, $max) = ($Hit->min - $clamp_length, $Hit->max + $clamp_length);
		$min = 0 if $min<0;
		$max = $seqlen if $max > $seqlen;
		$Hit->min($min);
		$Hit->max($max);
		push @Extended_Hits, $Hit;
	}
	return \@Extended_Hits;
}

1;

__END__

=head1 NAME

FileIO::FuzznucUtils - Provides utilities for reading fuzznuc files

=head1 SYNOPSIS

  use FileIO::FuzznucUtils;
  my $FuzznucHits = FileIO::FuzznucUtils->parse_fuzznuc_by_csv_file($fuzznuc_file);
  foreach my $Hit (@$FuzznucHits) {
	  print join("\t", $B->seq_id, $B->min, $B->max, $B->score, $B->strand, $B->pattern, $B->num_mismatches), "\n";
  	
  }

=head1 DESCRIPTION

This module provides some methods to read fuzznuc files and return a list of FuzznucHit
objects.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>

=head1 SEE ALSO

FileIO::FuzznucHit, FileIO::BarcodeUtils, FileIO::Barcode

=cut

package Grid::Tools::Barcode::Trimmer;
use strict;

=item $obj->trim_sequence;

B<Description:> This method returns the "clear range" trimmed sequence given a 
list of barcode hits and the length of the barcode clamp.

B<Parameters:> Accepts a sequence string, list of hits, and the length of the 
barcode clamp.

B<Returns:> a string representing the sequence after barcode trimming


=cut

sub trim_sequence {
	my ($self, $seq, $hits, $clamp_length) = @_;
	my ($start, $end) = $self->trim_clear_range($seq, $hits, $clamp_length);
	return substr($seq, $start, $end-$start);
}

=item $obj->trim_clear_range;

B<Description:> This method takes a listref of FileIO::FuzznucHit objects,
and extends the hit location on both 5' and 3' ends by the lenght of the barcode
clamp.

B<Parameters:> Accepts a barcode hashref, list of hits, and the length of the 
barcode clamp.

B<Returns:> a listref of merged FileIO::FuzznucHit objects extended in both directions
by the length of the clamp

Typical trimming cases

sff file:		sequence starts at base 33
| key (4bp) | barcode (22bp) | hexamer/clamp (6bp) | -------- sequence ------->
1           5                27                    33

fastq file:	sequence starts at base 29
| barcode (22bp) | hexamer/clamp (6bp) | -------- sequence ------->
1                23                    29

=cut

sub trim_clear_range {
	my ($self, $Deconvolver, $barcode, $seq, $hits) = @_;
	my $clamp_length = $Deconvolver->clamplength($barcode);
	$clamp_length = 0 unless $clamp_length > 0;
	my $key_length = $Deconvolver->keylength;
	$key_length = 0 unless $key_length > 0;
	
	# Define our trim points
	my ($start, $end);
	my $seqlen = length($seq);
	
	# If the user provides the key length but not the sequence,
	# then offset the sequence length by the key length
	# (since $F->seq is not prepended with the $Deconvolver->key)
	if (!$Deconvolver->key && $key_length) {
		$seqlen += $key_length;
		
	# Since we appended the key sequence, we have to decrement
	# the actual seqlen by the key_length
	} elsif ($Deconvolver->key) {
		$seqlen -= $key_length;
	}
	
	my $num_hits = scalar @$hits;
	if ($num_hits == 1) {
		my $hit = $hits->[0];
		
		# negative strand: we want the sequence 5' of the clamp-barcode
		# ==================================================clamp---Barcode---
		if ($hit->strand =~ /-/) {
			$start = 0;
			$end = $hit->min - $clamp_length;

			# Not sure why?
			if ($Deconvolver->key) {
				$start += $key_length;
			}
			
		# positive strand: we want the sequence 3' of the barcode-clamp
		# ---Barcode---Clamp==================================================
		} else {
			$start = $hit->max + $clamp_length;
			$end = $seqlen;
		}
		
	
	} elsif ($num_hits == 2) {
		
		# Sort hits by min position in ascending order
		my ($hit1, $hit2) = sort { $a->{min} <=> $b->{min} } @$hits;
		
		# we expect the first hit to be on the + strand
		#                   [******************************]
		# ---Barcode---Clamp================================Clamp---Barcode---
		if ($hit1->strand eq "+" && $hit2->strand eq "-") {
			$start = $hit1->max + $clamp_length;
			$end = $hit2->min - $clamp_length;
		
		# this is odd, but...
		#                                               [********************]
		# ---Barcode---Clamp==========---Barcode---Clamp======================
		} elsif ($hit1->strand eq "+" && $hit2->strand eq "+") {
			$start = $hit2->max + $clamp_length;
			$end = $seqlen;
		
		# this is odd, but...
		# [******************]
		# ====================Clamp---Barcode---============Clamp---Barcode---
		} elsif ($hit1->strand eq "-" && $hit2->strand eq "-") {
			$start = 0;
			$end = $hit1->min - $clamp_length;
			
			# Not sure why?
			if ($Deconvolver->key) {
				$start += $key_length;
			}
		}
		
	} elsif ($num_hits == 3) {
		
		# Sort hits by min position in ascending order
		my ($hit1, $hit2, $hit3) = sort { $a->{min} <=> $b->{min} } @$hits;
		
		# this is odd, but...
		#                                    [**************]
		# ---Barcode---Clamp---Barcode---Clamp==============Clamp---Barcode---
		if ($hit1->strand eq "+" && $hit2->strand eq "+" && $hit3->strand eq "-") {
			$start = $hit2->max + $clamp_length;
			$end = $hit3->min - $clamp_length;
		
		# this is odd, but...
		#                  [**************]
		# ---Barcode---Clamp==============Clamp---Barcode---Clamp---Barcode---
		} elsif ($hit1->strand eq "+" && $hit2->strand eq "-" && $hit3->strand eq "-") {
			$start = $hit1->max + $clamp_length;
			$end = $hit2->min - $clamp_length;
		}
		
	}
	
	# Offset start and end of our clear range by the key sequence
	if (!$Deconvolver->key && $key_length) {
		$start += $key_length;
		$end += $key_length;
	}
	
	# Let's provide an explanation for our being unable to trim off these barcodes
	my $reason;
	if (!defined $start) {
		if ($num_hits > 3) {
			$reason = "TOO_MANY_HITS";
		} else {
			$reason = "MISORIENTED_BARCODES";
		}
	} else {
		# Make sure we don't try to capture beyond the ends of the sequence
		$start = 0 if $start < 0; 
		$end = $seqlen if !defined $end || $end > $seqlen;
		
		if ($end < $start) {
			$reason = "NO_TRIMMED_SEQUENCE";
			$start = undef; # unable to provide trimmed sequence
		}
	}
	# Include the location of the barcode hits if trimming fails, for any reason
	if ($reason) {
		my @hits = sort { $a->{min} <=> $b->{min} } @$hits;
		my @locs;
		push @locs, $_->pattern.":".$_->min."..".$_->max.":".$_->strand.":".$_->num_mismatches foreach @hits;
		my $loc_info = join(" ", @locs);
		$reason .= "\t$loc_info";
	}
	
	return ($start, $end, $reason);
}


1;
__END__

=head1 NAME

Barcode::Trimmer

=head1 SYNOPSIS


  use Barcode::Trimmer;
  my $Trimmer = new Barcode::Trimmer ();
  
  # You can use FileIO::FuzznucHit and FileIO::Barcode objects
  $Trimmer->trim_clear_range($seq, $Barcode, $Hits)
  
  # Or, just give it hashrefs with the hit locations and Barcode information
  $Trimmer->trim_barcode_hits($seq, { }, [{min,max,strand}, ...}])
   
=head1 DESCRIPTION

This module provides function for barcode trimming of sequences.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut

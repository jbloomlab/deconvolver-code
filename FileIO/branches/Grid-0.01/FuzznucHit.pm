package FileIO::FuzznucHit;
use strict;
use base qw(Class::Accessor);
FileIO::FuzznucHit->mk_accessors(qw( seq_id min max length strand pattern num_mismatches barcode ));

=head2 to_csv_string()
	
	Returns a single-line string in the same format as the CSV output
	as the Fuzznuc executable.
	
=cut

sub to_csv_string() {
	my $Hit = shift;
	my @fields = ($Hit->seq_id, $Hit->min, $Hit->max, $Hit->length, $Hit->strand, $Hit->pattern, $Hit->num_mismatches);
	for (my $i=0; $i<scalar @fields; $i++) {
		$fields[$i] = "" unless defined $fields[$i];
	}
	return join("\t", @fields);
}

1;
__END__

=head1 NAME

FileIO::FuzznucHit - Models a fuzznuc hit, ie. an entry in the results of
a fuzznuc  results file. .

=head1 SYNOPSIS

  use FileIO::FuzznucHit;
  my $Hit = new FileIO::FuzznucHit ({ seq_id => "SOLEXA1_1_24_325_997#0/1", 
  						start => 67,
  						end => 47, 
 						length => 20
  						strand => "-",
  						pattern => "pattern1",
  						num_mismatches => 2
  });
  
  
=head1 DESCRIPTION

This module uses the Class::Accessor module to provide a quick solution to 
model a Fuzznuc hit

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut

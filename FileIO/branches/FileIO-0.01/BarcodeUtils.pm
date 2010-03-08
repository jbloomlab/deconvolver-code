package FileIO::BarcodeUtils;
use strict;
use FileIO::Barcode;

=item $obj->parse_barcode_by_file(");

B<Description:> This utility method uses the parse_barcode_by_filehandle, and
does the extra step of opening and closing the filehandle for you.

B<Parameters:> Accepts a path to a barcode file

B<Returns:> Hashref with keys as barcode ids and values of FileIO::Barcode objects

=cut


sub parse_barcode_by_file {
    my ($self, $barcode_file, $barcode_table) = @_;
    open(FH_BARCODE, "< $barcode_file") || die "Could not open barcode file $barcode_file for reading.\n";
    my $Barcodes = $self->parse_barcode_by_filehandle(*FH_BARCODE, $barcode_table);
    close FH_BARCODE;
    return $Barcodes;
}

=item $obj->parse_barcode_by_filehandle();

B<Description:> This utility method is used to parse a Solexa barcode file, and
return a list of FileIO::Barcode objects, where each object represents a barcode entry.

B<Parameters:> Accepts a filehandle to a barcode file.. Accepts a fasta_table
to make the function reentrant.

B<Returns:> Hashref with keys as barcode ids and values of FileIO::Barcode objects

=cut

sub parse_barcode_by_filehandle {
    my ($self, $fh, $barcode_table) = @_;
    while (<$fh>) {
    	    chomp;
    	    my @f = split /\t/, $_;
    	    $barcode_table->{$f[0]} = new FileIO::Barcode({ 
    	    		    id => $f[0], 
    	    		    seq => $f[1], 
    	    		    sample_id => $f[2], # also referred to as bac_id in the GLK
    	    		    blinded_id => $f[3],
    	    		    species_code => $f[4]
    	    });
    }
    return $barcode_table;
}


1;

__END__

=head1 NAME

FileIO::BarcodeUtils - Provides utilities for reading barcode files

=head1 SYNOPSIS

  use FileIO::BarcodeUtils;
  my $Barcodes = FileIO::BarcodeUtils->parse_barcode_by_file($barcode_file);
  foreach my $B (@$Barcodes) {
	  print join("\t", $B->id, $B->seq, $B->barcode_num, $B->project_id, $B->sample), "\n";
  	
  }

=head1 DESCRIPTION

This module provides some methods to read barcode files and return a list of 
FileIO::Barcode objects.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut

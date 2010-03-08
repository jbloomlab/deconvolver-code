package FileIO::Barcode;
use strict;
use base qw(Class::Accessor);
FileIO::Barcode->mk_accessors(qw( id seq sample_id blinded_id species_code clamp_length ));

1;
__END__

=head1 NAME

FileIO::Barcode - Models an entry in a Solexa barcode file.

=head1 SYNOPSIS

  use FileIO::Barcode;
  my $barcode = new FileIO::Barcode ({ id => "BC004CG", 
  						seq => "CGTAGTACACTCTAGAGCACTA",
  						sample_id => "26258", 
 						blinded_id => "NIGSP-HI-00031"
  						species_code => "Influenza A virus (A/northern shoveler/Washington/44249-604/2006(H N))"
  });
  
  # Uses Class::Accessor to provide setters for any of the attributes
  $barcode->clamp_length(6); 
  
  # Generate entry matching the format of the barcode file
  print join("\t", $barcode->id, $barcode->seq, $barcode->barcode_num, $barcode->project_id, $barcode->sample), "\n";
  
=head1 DESCRIPTION

This module uses the Class::Accessor module to provide a quick solution to 
model a Solexa barcode entry.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut

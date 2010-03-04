package FileIO::Fasta;
use strict;
use base qw(Class::Accessor);
FileIO::Fasta->mk_accessors(qw( id seq desc qual blinded_id ));

=item $obj->fasta_header();

B<Description:> This method returns a string of the header line.

B<Parameters:> None

B<Returns:> String

=cut

sub fasta_header {
    my $F = shift;
    return join(" ", ">", $F->id, $F->desc), "\n";
}
   
=item $obj->fasta_header();

B<Description:> This method returns a fasta formatted sequence.

B<Parameters:> None

B<Returns:> String

=cut

sub fasta_seq {
	my ($F, $seq, $width) = @_;
	$width = 60 unless $width && $width > 0;
	$seq = $F->seq unless $seq;
	my $length=length($seq);
	my $fasta_seq;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		$fasta_seq.= substr($seq, $pos, $out_width)."\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
	
	return $fasta_seq;
}


sub fasta_qual {
	my ($F, $qual, $width) = @_;
	$width = 60 unless $width && $width > 0;
	$width--;
	$qual = $F->qual unless $qual;
	my @quals = split / /, $qual;
	my $length = scalar @quals;
	
	# Add newlines every 60th value
	my @quals_block;
	for (my $i=0; $i<$length; $i++) {
		push @quals_block, "\n" if $i && !($i % 60);
		push @quals_block, "$quals[$i] ";
	}
	return join("", @quals_block, "\n");
}

=item $obj->fasta_record();

B<Description:> This method returns a string of the fasta record, including
the header and sequence..

B<Parameters:> None

B<Returns:> String

=cut

sub fasta_record {
    my $F = shift;
    return $F->header . $F->fasta_seq;
}

1;
__END__

=head1 NAME

FileIO::Fasta - Models an entry in a fasta file.

=head1 SYNOPSIS

  use FileIO::Fasta;
  my $F = new FileIO::Fasta({
	    id => "",
	    desc => "",
	    seq => "CGTAGTACACTCTAGAGCACTA",
  });

  # Generate a fasta record
  print $F->fasta_record;
  print join(" ", ">", $F->id, $fasta->desc), "\n", $F->fasta_seq;
  
=head1 DESCRIPTION

This module uses the Class::Accessor module to provide a quick solution to 
model a Solexa fasta entry.

=over 4

=cut

=back

=head1 BUGS

If you would like to report a problem with this module or would like to request
an enhancement, please submit a bug report to the author.

=head1 AUTHOR

Nelson Axelrod <naxelrod@jcvi.org>


=cut

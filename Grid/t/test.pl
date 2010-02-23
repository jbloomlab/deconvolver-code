#!/usr/local/bin/perl -w
use strict;
use IO::File;
use FileIO::FuzznucHit;
use FileIO::FuzznucUtils;

my @files = `ls /usr/local/scratch/naxelrod/out/fuzznuc.naxelrod.3474978.*`;

my $iter = FileIO::FuzznucUtils->hit_iterator(\@files);
#my $iter = hit_iterator(undef, \@files);

my $count=0;
while (my $h = $iter->()) {
	if ($h->seq_id) {
		$count++;
	}
}
print "total $count\n";
exit;


sub hit_iterator {
	my ($self, $files) = @_;
	
	# Open the first file for reading
	my $file = shift @$files;
	my $fh = IO::File->new("< $file") || die "Could not open fuzznuc results file $file.\n";

	my ($i, $done, $count);
	my $Iterator;
	$Iterator = sub {
		return undef if $done;
		
		# Iterate thru our filehandle, and return the next Hit
		while (<$fh>) {
			chomp;
			next if $_ =~ /^SeqName/; # skip headers
			$count++;
			
			my ($seq_id, $start, $end, $length, $strand, $pattern, $num_mismatches) = split /\t/, $_;
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
		
		print STDERR "$count $file";
		my $count = 0;
		
		# Open the next file for reading
		$file = shift @$files;
		if ($file) {
			$fh = IO::File->new("< $file") || die "Could not open fuzznuc results file $file.\n";
			return $Iterator->();
		} else {
			# Finished iterating thru all of the files
			$done = 1 unless scalar @$files;
			return undef;
		}
	};
	
	return $Iterator;
}



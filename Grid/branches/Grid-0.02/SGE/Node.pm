package Grid::SGE::Node;
use strict;
use base qw(Class::Accessor::Fast);
Grid::SGE::Node->mk_accessors(
	qw( 
		name
		type
		usage
		load
		arch
		states
	)
);

sub to_string {
	my $N = shift;
	# troubleshooting.q@dell-0-4-8.j I     0/2       0.00     lx26-eon64
	my @fields = ($N->name, $N->type, $N->usage, $N->load, $N->arch, $N->states);
	for (my $i=0; $i<scalar @fields; $i++) {
		$fields[$i] = "" unless defined $fields[$i];
	}
	return join("\t", @fields);
}

1;


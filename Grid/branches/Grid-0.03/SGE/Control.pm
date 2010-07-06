=head1 Grid::SGE::Control

	Wrapper for qdel 

=head1 AUTHORS

	Nelson Axelrod,
	naxelrod@jcvi.org

=cut

package Grid::SGE::Control;
use strict;
use Exporter;

use vars qw(@ISA @EXPORT_OK);
@ISA = qw( Exporter );
@EXPORT_OK = qw(delete_jobs);
our $VERSION = '0.01';


=head2 delete_jobs()

	delete all jobs given in list
	Returns 1 on success, 0 on failure

=cut
sub delete_jobs {
	my ($self, @jobs) = @_;
	if (!scalar @jobs) {
		warn "Error: delete_jobs function requires a list of specified jobs.\n";
		return;
	}
	my $qdel = $self->executable('qdel');
	my $command = join(" ", $qdel, @jobs);
	my $status_code = system($command);
	if ($status_code) {
		warn "Error: qdel returns with a status code $status_code (command: $command)\n";
	} else {
		
	}
	return ($status_code) ? 0 : 1;
}

1;

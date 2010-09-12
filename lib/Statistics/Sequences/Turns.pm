package Statistics::Sequences::Turns;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);
use Scalar::Util qw(looks_like_number);
use Statistics::Zed 0.02;
use Statistics::Sequences 0.05;
@ISA = qw(Statistics::Sequences);

$VERSION = '0.01';

=pod

=head1 NAME

Statistics::Sequences::Turns - Kendall's test for turning-points - peaks or troughs - in a numerical sequence

=head1 SYNOPSIS

 use Statistics::Sequences::Turns;
 $turns = Statistics::Sequences::Turns->new();
 $turns->load(0, 3, 9, 2 , 1, 1, 3, 4, 0, 3, 5, 5, 5, 8, 4, 7, 3, 2, 4, 3, 6);
 $turns->test()->dump();
 #Z = -0.0982471864864821, 2p = 0.92174

=head1 DESCRIPTION

This module implements a test of randomness that is suitable for data of the continuous numerical type - not static categories (like choices between a "banana" and "cheese"), but sequences of numerical values where it is meaningful to speak of values, on trial I<i>, that might be higher or lower than the value on trials I<i> - 1 and I<i> + 1; the trial's neighbours. It is particularly commended for time-series, testing if a numerical sequence shows systematic rather than random oscillations.

Specifically, the test concerns whether the value on trial I<i> (for I<i> is greater than zero and less than I<n>), with respect to its neighbours, is a peak (greater than both neighbours) or a trough (less than both neighbours), and if the frequencies of these turns, as peaks and troughs, is commensurate with what is expected for a randomly generated sequence. In this way, the Turns-test is always based on three consecutive values in a sequence.

If you have one or two sets of categorical data, or two groups of numerical data, you can firstly dichotomize them into an array of zeroes and ones (see L<Statistics::Sequences/Dichotomising data|Statistics::Sequences/Dichotomising data>), and then perform the Turns test to assess the randomness of their sequential association.

=head1 METHODS

=head2 new

 $turns = Statistics::Sequences::Turns->new();

Returns a new Turns object. Expects/accepts no arguments but the classname.

=head2 load

 $turns->load(@data);
 $turns->load(\@data);
 $turns->load('dist1' => \@data1, 'dist2' => \@data2)
 $turns->load({'dist1' => \@data1, 'dist2' => \@data2})

Loads data anonymously or by name. See L<load|Statistics::Sequences/load> in the Statistics::Sequences manpage.

=head2 test

 $turns->test();

Performs Kendall's turning-points test on the given or named distribution, yielding a I<Z> statistic.

The number of turns as peaks and troughs are then counted up, starting from element 1, checking if both its left/right (or past/future) neighbours are lesser than it (a peak) or greater than it (a trough). Wherever the values in successive indices of the list are equal, they are treated as a single observation/datum -so the following:

 0 0 1 1 0 1 1 1 0 1

is counted up for turns as

 0 1 0 1 0 1

So, e.g., there are four turns in the above example - two peaks (0 1 0) and two troughs (1 0 1). (This would not be picked up as a non-random sequence, but if it were repeated, it would be seen to significantly deviate from expectation, I<p> = .035.)

The observed number of turns is compared to the number expected, and this deviation is assessed against the expected deviation, i.e., as a Z-value; Kendall (1973) having observed that the statistic shows "a fairly rapid tendency of the distribution to normality" (p. 24). 

=cut

#-----------------------------------------------
sub test {
#-----------------------------------------------
    my $self = shift;

    # Initialize arguments as given or by default:
    my $args = ref $_[0] ? $_[0] : {@_};
    $self->_testdata_aref($args);
    my @data = @{$self->{'testdata'}};

    # Remove equivalent successors: e.g., strip 2nd 2 from (3, 2, 2, 7, 2) (Note: List::MoreUtils::uniq would strip the final 2 as well)
    # Also check that are elements are numeric:
    my ($i, $s, @data_u) = ();
    my %seen   = ();
    for ($i = 0; $i < scalar(@data); $i++) {
        croak __PACKAGE__, '::test All data must be numerical for testing turning-points' if ! looks_like_number($data[$i]);
        push @data_u, $data[$i] if !scalar(@data_u) || $data[$i] != $data_u[-1];
    }
    my $n = scalar(@data_u);
    croak __PACKAGE__, '::test Insufficient data to conduct a test of turning-points' if $n < 4;

    # Having passed all possible croaks, continue with default/offered attribs:
    $self->{'_tails'} = $args->{'tails'} || 2;
    
    # Compute observed number of turns:
    my ($obs_val) = (0);
    for ($i = 1; $i < $n - 1; $i++) {
        if ( ($data_u[$i - 1] > $data_u[$i]) && ($data_u[$i + 1] > $data_u[$i]) ) { # we have a trough at $i
            $obs_val++;
        }
        elsif ( ($data_u[$i - 1] < $data_u[$i]) && ($data_u[$i + 1] < $data_u[$i]) ) { # we have a peak at $i
            $obs_val++;
        }
    }
    
    # Compute expected number of turns and its variance:
    my $exp_val = 2/3*($n - 2);
    my $var = (16*$n - 29) / 90;  
  
    # Test significance and lump various values into the class object:
    $self->_expound($obs_val, $exp_val, $var, $args);

    return $self;
}

=head2 dump

 $turns->dump(flag => '1|0', text => '0|1|2');

Print test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details.

=cut

#-----------------------------------------------
sub dump {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $args->{'testname'} = 'Turns';
    $self->SUPER::_dump_pass($args);
    return $self;
}

__END__

=head1 EXAMPLE

=head2 Seating at the diner

This is the data from Swed and Eisenhart (1943) also given as an example for the L<Runs test|Statistics::Sequences::Runs/EXAMPLE> and L<Vnomes test|Statistics::Sequences::Vnomes/EXAMPLE>. It lists the occupied (O) and empty (E) seats in a row at a lunch counter.
Have people taken up their seats on a random basis? The Runs test suggested some non-random basis for people to take their seats, ouputting (as per C<dump>):

  Runs: observed = 11.00, expected = 7.88, Z = 1.60, 1p = 0.054834

That means there was more serial discontinuity than expected. What does the test of Turns tell us?

 use Statistics::Sequences::Turns;
 my $turns = Statistics::Sequences::Turns->new();
 my @seating = (qw/E O E E O E E E O E E E O E O E/);
 $turns->load(\@data);
 $turns->binate(); # transform Es and Os into 1s and 0s
 $turns->test(tails => 1)->dump();

This outputs, as returned by C<string>: 

 Z = 1.95615199108988, 1p = 0.025224

So each seated person is neighboured by empty seats, and/or each empty seat is neighboured by seated persons, more so than would be expected if people were taking their seats randomly.

=head1 REFERENCES

Kendall, M. G. (1973). I<Time-series>. London, UK: Griffin. [The test is described on pages 22-24. Note that in the Example 2.1 for this test, the variable used in the calculation of the expected number of turns should be 52 (i.e., I<n> - 2), not 54.]

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> for other tests of sequences, and for sharing data between these tests.

=head1 TO DO/BUGS

Implementation of the serial test for non-overlapping I<v>-nomes.

=head1 REVISION HISTORY

See CHANGES in installation dist for revisions.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2010 Roderick Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=back

=head1 DISCLAIMER

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=head1 END

This ends documentation of the Perl implementation of Kendall's turning-points test for randomness of a numerical sequence.

=cut

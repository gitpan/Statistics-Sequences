package Statistics::Sequences::Runs;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);
use Statistics::Sequences 0.051;
@ISA = qw(Statistics::Sequences);

$VERSION = '0.051';

=pod

=head1 NAME

Statistics::Sequences::Runs - The Runs-test (Wald-Walfowitz or Swed-Eisenhard Test)

=head1 SYNOPSIS

 use Statistics::Sequences::Runs;
 $runs = Statistics::Sequences::Runs->new();
 $runs->load(qw/1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 1 1 0 1/);
 $runs->test()->dump();

=head1 DESCRIPTION

The Runs-test assesses the difference between two independent distributions, or a difference within a single distribution of dichotomous observations, in terms of the frequency of the runs of states within them.

A run is a sequence of identical states on 1 or more consecutive trials. For example, in a signal-detection test, there'll be a series, over time, of hits (H) and misses (M), which might look like H-H-M-H-M-M-M-M-H. Here, there are 5 runs: 3 of hits, and 2 of misses. This number of runs can be compared with the number expected to occur by chance, given the number of observed hits and misses. More runs than expected ("negative serial dependence") generally indicates irregularity, or instability; fewer runs than expected ("positive serial dependence") indicates regularity, or stability. Both can indicate a sequential dependency: either negative (an extra-chance factor, or bias, to produce too many alternations), or positive (an extra-chance factor, or bias, to produce too many repetitions).

The distribution of runs is asymptotically normal - quite quickly, with probabilities well estimated by the normal distribution when both the numbers of H and M exceed 10 (e.g., Kelly, 1982). The deviation of the observed number of runs is therefore reliably assessed by way of a I<Z>-score.

=head1 METHODS

Methods are essentially as described in L<Statistics::Sequences>. See this manpage for how to handle non-dichotomous data, e.g., numerical data, or those with more than two categories.

=head2 new

 $runs = Statistics::Sequences::Runs->new();

Returns a new Runs object. Expects/accepts no arguments but the classname.

=head2 load

 $runs->load(@data);
 $runs->load(\@data);
 $runs->load('dist1' => \@data1, 'dist2' => \@data2)
 $runs->load({'dist1' => \@data1, 'dist2' => \@data2})

Loads data anonymously or by name. See L<load|Statistics::Sequences/load> in the Statistics::Sequences manpage.

=head2 test

 $runs->test();

Performs the runs test on the named distributions. If only one distribution name is given, the "one-sample" Runs test is performed, cutting the data at the median, or by the value given as cut. Observations that fall above and below the cut-value then constitute the "groups" to be searched for runs. Otherwise, with two named groups, runs are sought on the basis of the observations belonging to one or the other named group.

=cut

#-----------------------------------------------
sub test {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};

    $self->_testdata_aref($args);

    my @data = @{$self->{'testdata'}};
    my ($lim, $obs_val, $i, %states) = (scalar(@data), 0);
    
    # Count the observed number of runs - while initing a hash giving the frequency of each state:
    foreach ($i = 0; $i < $lim; $i++) {
        $states{$data[$i]}++;
        $obs_val++ if !$i || ( $data[$i] ne $data[$i - 1] );
    }

    #croak __PACKAGE__, '::test Less than two states were found in the data' if scalar keys %states < 2;
    return undef if ! scalar keys %states;
    ##croak __PACKAGE__, '::test No states/no test data appear to be loaded' if ! scalar keys %states;
    croak __PACKAGE__, '::test More than two states were found in the data: ' . join(' ', keys(%states)) if scalar keys %states > 2;

    # Get the probability of these many runs:
    my ($m, $n) = keys %states;
    $n ||= 0; # account for only 1 state having appeared (i.e, a run of 1)
    my ($n1, $n2) = ($states{$m}, $states{$n});
    $n2 ||= 0;
    my $sum = $n1 + $n2; # the sum of all observations

    $self->_expire($obs_val, 1, $args) if $sum > 1 && scalar keys %states > 1; # for a valid test, need at least 2 states in 2 trials

    # Calc expectation value and variance:
    my ($exp_val, $var) = ();
    $exp_val = ( ( 2 * $n1 * $n2 ) /  $sum ) + 1;
    $var = ( 2 * $n1 * $n2 * ( ( 2 * $n1 * $n2 ) - $sum) ) 
            /
           ( ( $sum**2 ) * ( $sum - 1 ) );

    # Test significance and lump various values into the class object:
    $self->_expound($obs_val, $exp_val, $var, $args);

    return $self;
}

=head2 dump

 $runs->dump(flag => '1|0', text => '0|1|2');

Print Runs-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details.

=cut

#-----------------------------------------------
sub dump {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $args->{'testname'} = 'Runs';
    if ($args->{'text'} and $args->{'text'} > 1) {
        $args->{'title'} = "Runs test results:";
        $self->SUPER::_dump_verbose($args);
    }
     else {
        $self->SUPER::_dump_sparse($args);
    }
    return $self;
}

__END__

=head1 EXAMPLE

=head2 Seating at the diner

Swed and Eisenhart (1943) list the occupied (O) and empty (E) seats in a row at a lunch counter.
Have people taken up their seats on a random basis?
There is no need to dichotomise these data: there is already a single sample, with dichotomous, categorical observations.

 use Statistics::Sequences::Runs;
 my $runs = Statistics::Sequences::Runs->new();
 my @seating = (qw/E O E E O E E E O E E E O E O E/);
 $runs->load(\@seating);
 $runs->test(ccorr => 1, tails => 1)->dump();

Suggesting some non-random basis for people taking their seats, this outputs: 

 Runs: observed = 11.00, expected = 7.88, Z = 1.60, 1p = 0.054834 

These data are also used as examples of the L<Turns test|Statistics::Sequences::Turns/EXAMPLE> and the L<Vnomes test|Statistics::Sequences::Vnomes/EXAMPLE>.

=head2 ESP runs

In a single run of a classic ESP test, there are 25 trials, each composed of a randomly generated state (typically, one of 5 possible geometric figures), and a human-generated state drawn from the same pool of alternatives. Tests of the synchrony between the random and human data are then made, typically in terms of the number of "hits" observed versus that expected. The runs of hits and misses can also be tested by dichotomising the data on the basis of the L<match|Statistics::Sequences/match> of the random "targets" with the human "responses", like so:

 use Statistics::Sequences::Runs;

 # Produce pseudo ESP targets and responses:
 my ($i, @targets, @responses);
 for ($i = 0; $i < 250; $i++) {
    $targets[$i] = (qw/circle plus square star wave/)[int(rand(5))];
    $responses[$i] = (qw/circle plus square star wave/)[int(rand(5))];
 }

 # Test for runs of matches between targets and responses:
 my $runs = Statistics::Sequences::Runs->new();
 $runs->load(targets => \@targets, responses => \@responses);
 $runs->match(data => [qw/targets responses/]);
 $runs->test();
 print "The probability of obtaining these $runs->{'observed'} runs is $runs->{'p_value'}\n";

 # But let's test (preferably, if predicted) that the responses were matched to the target on the trial one ahead (as if by "precognition"):
 $runs->match(data => [qw/targets responses/], lag => 1)->test();
 print "With responses synchronised to targets on the next (+1) sample,\n 
 $runs->{'observed'} runs in 250 samplings were produced when $runs->{'expected'} were expected,\n 
 a deviation with an associated probability of $runs->{'p_value'}\n"; 

=head1 REFERENCES

Kelly, E. F. (1982). On grouping of hits in some exceptional psi performers. I<Journal of the American Society for Psychical Research>, I<76>, 101-142.

Swed, F., & Eisenhart, C. (1943). Tables for testing randomness of grouping in a sequence of alternatives. I<Annals of Mathematical Statistics>, I<14>, 66-87. [Look in C<ex/checks.pl> in the installation dist for a few examples from this paper for testing.]

Wald, A., & Wolfowitz, J. (1940). On a test whether two samples are from the same population. I<Annals of Mathematical Statistics>, I<11>, 147-162.

Wolfowitz, J. (1943). On the theory of runs with some applications to quality control. I<Annals of Mathematical Statistics>, I<14>, 280-288. [Suggests some ways in which data may be dichotomised for testing runs.]

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> for other tests of sequences, and for sharing data between these tests.

=head1 TO DO/BUGS

Results are dubious if there are only two observations.

Testing not by I<z>-scores, and/or using poisson distribution for low number of observations

Fu's Markovian solution

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

This ends documentation of a Perl implementation of the Wald-Walfowitz Runs test for randomness and group differences within a sequence.

=cut

package Statistics::Sequences::Vnomes;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);
use Algorithm::Combinatorics qw(variations_with_repetition);
use Math::Cephes 0.43;
use Statistics::Sequences 0.051;
use Statistics::Lite qw(sum);
@ISA = qw(Statistics::Sequences);

$VERSION = '0.04';

=pod

=head1 NAME

Statistics::Sequences::Vnomes - The Serial Test (psi-square) and Generalized Serial Test (delta psi-square) for equiprobability of I<v>-nomes (or I<v>-plets/bits) (Good's and Kendall-Babington Smith's tests)

=head1 SYNOPSIS

 use Statistics::Sequences::Vnomes;
 $vnomes = Statistics::Sequences::Vnomes->new();
 $vnomes->load(\@data); # categorical
 $vnomes->test(length => 2)->dump();

=head1 DESCRIPTION

This module implements tests of the independence of successive elements of a sequence/series of data (list, vector, etc.) - specifically, "serial tests" for I<v>-nomes (a.k.a I<v>-plets or, for binary data, I<v>-bits) - what are call singlets/monobits, dinomes/doublets, trinomes/triplets, etc..

Serial tests tell us if all the variations of the states, of a certain sub-sequence length, I<v>, that would be possible in the population from which the series has been sampled, are equally represented in the sample. For example, a sequence sampled from a "heads'n'tails" (H and T) distribution can be tested for its equal representation of the trinomes HTH, HTT, TTT, THT, and so on. Counting up these I<v>-nomes at all points in the series, permitting overlaps, yields a statistic - psi-square - that is approximately distributed as I<chi>-square; the Kendall-Babington Smith statistic.

However, because these counts are not independent (given the overlaps), Good's Generalized Serial Test is the default test-statistic returned by this module's C<test> routine: It computes psi-square by differencing, viz., in relation to not only the specified C<length>, or value of I<v>, but also its value for the first two prior lengths of I<v>, yielding a statistic, delta-square-psi-square (the "second backward difference" measure) that is I<exactly> distributed as I<chi>-square.

The test is suitable for multi-state data, not only the binary, dichotomous series suitable for the Runs and Joins tests in this package. It can also be used to test that the individual elements in the list are uniformly distributed, that the states are equally represented, i.e., as a I<chi>-square-based frequency test (a.k.a. test of uniformity, equiprobability, equidistribution). This is done by giving a C<length> of 1, i.e., testing for mononomes.

Note that this is I<not> the so-called serial test described by Knuth (1998, Ch. 2), which involves non-overlapping pairs of sequences. Given this variety of definitions of what is a "serial test," this module - like that for Runs, Pot, etc. - is named after the basic construct tested - i.e., I<v>-nomes - rather than after any property of I<v>-nomes (seriality, successive independence, etc.).

=head1 METHODS

=head2 new

 $vnomes = Statistics::Sequences::Vnomes->new();

Returns a new Vnomes object. Expects/accepts no arguments but the classname.

=head2 load

 $vnomes->load(@data);
 $vnomes->load(\@data);
 $vnomes->load('dist1' => \@data1, 'dist2' => \@data2)
 $vnomes->load({'dist1' => \@data1, 'dist2' => \@data2})

Loads data anonymously or by name. See L<load|Statistics::Sequences/load> in the Statistics::Sequences manpage.

=head2 test

 $vnomes->test(length => ?integer?, delta => '1|0', circularize => '1|0', states => [qw/A C G T/]);

Performs the Kendall-Babington Smith or (by default) Good's Generalized Serial Test of I<v>-nomes on the given or named distribution, yielding a I<psi-square> statistic.

To test for the significance of the psi-square statistic, the raw psi-square value for sub-sequences of length I<v> is, by default, I<not> used - because, unless I<length> (I<v>) = 1, psi-square is not asymptotically distributed as I<chi>-square. However, the differences between psi-square values for backwardly adjacent values of I<length> (I<v>) are asymptotically distributed as I<chi>-square. By default, then, a "second backward differences" psi-square value is calculated, named (as per Good, 1953) as C<del2psi2>, which makes use of the psi-square values for sub-sequences of length I<v>, I<v> - 1, and I<v> - 2. This statistic is logically (and has been mathematically shown to be) not only I<chi>-square distributed, but to offer statistically independent counts of all the possible variations of sequences of I<length> for the series in question. 

Note that the "I<first> backward differences" of psi-square, which is the difference between the psi-square values for sub-sequences of length I<v> and length I<v> - 1, is not ordinarily returned. While it is I<chi>-square distributed, counts of such first-differences are not statistically independent (Good, 1953; Good & Gover, 1967). This value can, however, be returned in place of the default by specifying C<delta> => 1. But note: "the sequence of second differences forms a much better set of statistics for testing the hypothesis of flat-randomness" (Good & Gover, 1967, p. 104) [compared to the first differences].

For consistency with other modules in the L<Statistics::Sequences> package, the I<Z>-value associated with psi-square's probability is also computed, and this is what is printed in dumps. You can get the psi-square statistic itself, however, as $vnomes->{'psisq'}, noting that it will, according to your value of C<delta>, either be properly called I<del2psi2>, I<delpsi2> or I<psi2>.

The algorithm implemented for psi-square is that given by Good (1953, Eq. 1); benchmarking shows no reliable speed difference to alternative forms of the equation, as given by Good (1957, Eq. 2) and in the NIST test suite (Rukhin et al., 2001). Good's original algorithm can also be found in individual papers describing the application of the Serial Test (e.g., Davis & Akers, 1974).

By default, the I<p>-value associated with the test-statistic is 2-tailed. See the L<Statistics::Sequences|Statistics::Sequences/test> manpage for generic options other than the following Vnome test-specific ones. At the end of the test, the class object is lumped with the usual statistics; this time, however, the value of I<observed> is the average of the observed frequencies of each I<v>-nome, and an additional statistic, I<observed_stdev>, the standard deviation of the observed frequencies is also formed.

=head3 Options

=over 4

=item length

The length of the I<v>-nome, i.e., the value of I<v>. Must be an integer greater than or equal to 1, and smaller than than the sample-size.

What is a meaningful maximal value of C<length>? As a I<chi>-square test, it could be held that there should be an expected frequency of at least 5 for each I<v>-nome. This is "conventional wisdom" recommended by Knuth (1988) but can be judged to be too conservative (Delucchi, 1993). The NIST documentation on the serial test (Rukhin et al., 2001) recommends that length should be less than the floored value of log2 of the sample-size, minus 2. No tests are here made of these recommendations, but if you choose to "dump" your results with verbosity (see C<dump>), you will get a note if the NIST warning would apply.

=item circularize

By default, circularizes the data series; i.e., the datum after the last element is the first element. This affects (and slightly simplifies) the calculation of the expected frequency of each I<v>-nome, and so the value of each psi-square. Circularizing ensures that the expected frequencies are accurate; otherwise, they might only be approximate. As Good and Gover (1967) offer, "It is convenient to circularize in order to get exact checks of the arithmetic and also in order to simplify some of the theoretical formulae" (p. 103).

=item delta

By default, the statistics are based on the second backward difference of psi-squares, i.e., as the Generalized Serial Test, as described by Good, see L<REFERENCES|Statistics::Sequences::Vnomes/REFERENCES>. If I<delta> => 0, the original Kendall-Babington Smith statistic (psi-square) is used. If I<delta> => 1, the value of C<psisq> is the first backward difference (delpsi2), as described above.

=item states

A referenced array listing the unique states (or 'events' or 'letters') in the population from which the series was sampled. This is useful to specify if the series itself is likely not to include all the possible states; it might even include only one of them. If this array is not specified, the unique states are identified from the series itself - in which case there ought to be at least two states in the series. Having only one state specified is not permissible. If giving a list of states, a check in each test is made to ensure that the data series contains I<only> those elements in the list.

=back

=cut

#-----------------------------------------------
sub test {
#-----------------------------------------------
    my $self = shift;

    # Initialize arguments as given or by default:
    my $args = ref $_[0] ? $_[0] : {@_};
    my $v = delete $args->{'length'} || croak __PACKAGE__, '::test Must have a v-nome length greater than zero';
    $self->_testdata_aref($args);
    my @data = @{$self->{'testdata'}};
    my $n = scalar(@data);
    # Sanitize length in proportion to size:
    croak __PACKAGE__, '::test Sequence length v for testing should be no more than the size of the sample of data - 2' if $v >= scalar(@data) - 1;
    my $lim =  Math::Cephes::floor(Math::Cephes::log2($n)) - 2; # NIST recommmendation
    # Could warn here if $v is too small for $n amd $nstates if $v >= $lim ... At present, only use $lim for "verbose" dumps.

    # Having passed all possible croaks, continue with default/offered attribs:
    my $delta = defined $args->{'delta'} ? delete $args->{'delta'} : 2;
    my $circularize = defined $args->{'circularize'} ? delete $args->{'circularize'} : 1;
    $self->{'_tails'} = $args->{'tails'} || 2;

    # Get a list of unique states, as explicitly given or from the data sequence itself:
    my ($nstates, $states_aref) = _states(\@data, $args->{'states'});

    # Init a hash to keep the psi-square values for the v, v-1, and v-2 ( = $v_i) sequence lengths, where relevant:
    my (%stats, %stats_d) = ();
    $stats_d{$v}->{'psisq'} = 0;
    $stats_d{$v - 1}->{'psisq'} = 0 if $v >= 2;
    $stats_d{$v - 2}->{'psisq'} = 0 if $v >= 3;

    # While looping through tests of v, v-1 and v-2 (= $v_i) sequence lengths ...
    my ($v_i, $n_i, $x_i, $nx_i, @data_i, $freq) = ();
    foreach $v_i(keys %stats_d) {
        # Extend the data series by appending the first $v_i - 1 bits to the end of the series:
        @data_i = @data;
        push @data_i, @data[0 .. $v_i - 2] if $circularize; # e.g., if v = 3, append elements [0] & [1]; if v = 2, append only 1st value
        $n_i = scalar(@data);
 
        # Count up frequencies of each sequence of length $v_i:
        $freq = _frequencies(\@data_i, $v_i, $states_aref);
 
        # Compute expected number of any form of sequences of length $v_i:
        $nx_i = $circularize ? $n_i : ($n - $v_i + 1);
        $x_i = $nx_i * ( $nstates**(-1 * $v_i) );

        # Compute psi^2 itself for this v-nome length ($v_i) (Good, 1953, Eq. 1):
        $stats{$v_i}->{'psisq'} = ( ($nstates**$v_i) / $n_i) * sum( map{ ($_ - $x_i)**2 } values %{$freq});
 
        # Cache values that might be needed post-loop:
        $stats{$v_i}->{'expected'} = $x_i;
        $stats{$v_i}->{'obs_counts'} = $freq;
    }

    # Assess significance of the psi-square (or chi-square) value:
    my ($psisq, $df, $p_value) = _assess($v, $nstates, \%stats, $delta, $args->{'orig'});
    $p_value /= 2 if defined $args->{'tails'} and $args->{'tails'} == 1;

    # Lump values into class object:
    $self->{'lim'} = $lim;
    $self->{'length'} = $v;
    $self->{'samplings'} = $n;
    $self->{'nstates'} = $nstates;
    $self->{'delta'} = $delta;
    $self->{'counts'} = $stats{$v}->{'obs_counts'};
    $self->{'observed'} = Statistics::Lite::mean(values %{$self->{'counts'}});
    $self->{'observed_stdev'} = Statistics::Lite::stddev(values %{$self->{'counts'}});
    $self->{'expected'} = $stats{$v}->{'expected'};
    $self->{'obs_dev'} = $self->{'observed'} - $self->{'expected'};
    $self->{'p_value'} = $p_value;
    $self->{'df'} = $df;
    $self->{'_statname'} = $delta ? $delta == 2 ? 'del2psi2' : 'delpsi2' : 'psi2';
    $self->{'psisq'} = $psisq;
    
    $self->{'z_value'} = $self->{'zed'}->p2z(value => $p_value, tails => $args->{'tails'});

    $self->{'_tested'} = 1;

    return $self;
}

=head2 dump

 $vnomes->dump(flag => '1|0', text => '0|1|2');

Print Vnome-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details. For comparability with other modules in the Statistics::Sequences package, the I<Z>-value associated (by Math::Cephes::ndtri) with the obtained I<p>-value is reported. If C<text> => 2, then you get a verbose dump, including (1) the actual test-statistic depending on the value of C<delta> tested (C<del2psi2> for the second difference measure (default), C<delpsi2> for the first difference measure, and C<psi2> for the raw measure), followed by degrees-of-freedom in parentheses; and (2) a warning, if relevant, that your C<length> value might be too large with respect to the sample size (see NIST reference, above, in dicussing C<length>). If C<text> => 1, you just get the average observed and expected frequencies for each I<v>-nome, the I<Z>-value, and its associated I<p>-value.

If you want to access the actual psi-square value that was tested against the I<chi>-square distribution, you can retrieve it as $vnomes->{'psisq'}.

After testing, parameters named 'nstates' (the number of states), 'samplings' (the size of the sample), 'length' (what you requested) can be retrieved from the class object. You can retrieve the counts for each of the Vnomes in the series as a hash-reference named 'counts' in the class object, e.g.:

 print "No. of $vnomes->{'length'}-nome variations of $vnomes->{'nstates'} states among $vnomes->{'samplings'} samplings:\n";
 printf("\t%s\t%d\n", $_, $vnomes->{'counts'}->{$_}) foreach sort keys %{$vnomes->{'counts'}};

=cut

#-----------------------------------------------
sub dump {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};

    $args->{'testname'} = "Vnomes($self->{'length'})";
 
    if ($args->{'text'} and $args->{'text'} > 1) { # verbosity
        my $n_vars = scalar(keys(%{$self->{'counts'}}));
        $args->{'title'} = "Vnome test results:\n$n_vars $self->{'length'}-nome variations of $self->{'nstates'} states among $self->{'samplings'} observations:";
        $self->SUPER::_dump_verbose($args);
        printf(" %s (df = %d) = %f\n", $self->{'_statname'}, $self->{'df'}, $self->{'psisq'});
        # "Warn" if the length to test was too long for the number of observations (as per NIST recommendation):
        if ($self->{'length'} > $self->{'lim'}) {
            my $stub = $self->{'lim'} == 1 ? 'of' : 'less than or equal to';
            print " Note: The V-nome length ($self->{'length'}) could be too large for the sample-size ($self->{'samplings'}); a length $stub $self->{'lim'} is optimal\n" ;
        }
     }
     else {
        $self->SUPER::_dump_sparse($args);
    }
    return $self;
}

sub _states {
    my ($data, $states) = @_;
    if (! ref $states) { # Get states from the data themsevles:
        my %hash   = map { $_, 1 } @{$data};
        $states = [keys %hash];
    }
    else { # Ensure that the data only contain states in the given list:
        my ($g, $h) = ();
        DATA:
        foreach $g(@{$data}) {
            foreach $h(@{$states}) {
                if ($h eq $g) {
                    next DATA;
                }
            }
            croak __PACKAGE__, "::test The element $g in the data is not represented in the stipulated states"; 
        }
    }
    my $nstates = scalar(@{$states});
    croak __PACKAGE__, '::test At least two different values must be in the series to test its sub-sequences' if $nstates <= 1;
    return ($nstates, $states);
}

sub _frequencies {
    my ($data_i, $v_i, $states_aref) = @_;

    # Get a list of all possible combinations of states at the current length ($v_i):
    my @variations = variations_with_repetition($states_aref, $v_i);

    # Count up the frequency of each variation in the data:
    my ($i, $probe_str, $test_str, %freq) = ();
    foreach (@variations) {
        $probe_str = join'', @{$_};
        $freq{$probe_str} = 0;
        for ($i = 0; $i < scalar(@{$data_i}) - $v_i + 1; $i++) {
            $test_str = join'', @{$data_i}[$i .. ($i + $v_i - 1)];
            $freq{$probe_str}++ if $probe_str eq $test_str;
        }
    }
    return \%freq;
}
 
sub _assess {# Test significance
    my ($v, $nstates, $stats, $delta, $orig) = @_;
    my ($psisq, $df, $p_value) = ();
    
    if ($v == 1) { # psisq is asymptotically distributed chisq, can use psisq for chisq distribution:
        $psisq = $stats->{1}->{'psisq'};
        $df = $nstates - 1;
    }
    else {
        if ($delta == 2) { # Second backward difference (default):
            $psisq = $stats->{$v}->{'psisq'} - ( 2 * ($v - 1 <= 0 ? 0 : $stats->{$v - 1}->{'psisq'}) ) + ($v - 2 <= 0 ? 0 : $stats->{$v - 2}->{'psisq'});
            $df = ( $nstates**($v - 2) ) * ( $nstates - 1)**2;
        }
        elsif ($delta == 1) {# First backward difference:
            $psisq = $stats->{$v}->{'psisq'} - ($v - 1 <= 0 ? 0 : $stats->{$v - 1}->{'psisq'});
            $df = $nstates**$v - $nstates**($v - 1);
        }
        else { # Raw psisq:
            $psisq = $stats->{$v}->{'psisq'};
            $df = $nstates**$v - 1; #$v;
        }
    }

    $p_value = Math::Cephes::igamc($df/2, $psisq/2);

    return ($psisq, $df, $p_value);
}

__END__

=head1 EXAMPLE

=head2 Seating at the diner

This is the data from Swed and Eisenhart (1943) also given as an example for the L<Runs test|Statistics::Sequences::Runs/EXAMPLE> and L<Turns test|Statistics::Sequences::Turns/EXAMPLE>. It lists the occupied (O) and empty (E) seats in a row at a lunch counter.
Have people taken up their seats on a random basis? The Runs test suggested some non-random basis for people to take their seats, ouputting (per dump with text => 1):

  Runs: observed = 11.00, expected = 7.88, z = 1.60, 1p = 0.054834

So there were more runs, or greater serial discontinuity, than expected. What does the test of Vnomes tell us?

 use Statistics::Sequences::Vnomes;
 my $vnomes = Statistics::Sequences::Vnomes->new();
 my @seating = (qw/E O E E O E E E O E E E O E O E/);
 $vnomes->load(\@seating);
 $vnomes->test(length => 2)->dump();

This outputs, as returned by C<string>: 

 Z = -0.475232849247084, 2p = 0.317310507862914

That is, the observed frequency of each possible pair of seating arrangements (the dinomses OO, OE, EE, EO) did not differ significantly from that expected. Taking a bigger picture and widening our perspective, though, changing the value of C<length> to 3 so as to assess the equal representation of all possible I<tri>nomes, yields:

 Z = -1.70672129474387, 2p = 0.0439369336234074

=head1 REFERENCES

Davis, J. W., & Akers, C. (1974). Randomization and tests for randomness. I<Journal of Parapsychology>, I<38>, 393-407.

Delucchi, K. L. (1993). The use and misuse of chi-square: Lewis and Burke revisited. I<Psychological Bulletin>, I<94>, 166-176.

Good, I. J. (1953). The serial test for sampling numbers and other tests for randomness. I<Proceedings of the Cambridge Philosophical Society>, I<49>, 276-284.

Good, I. J. (1957). On the serial test for random sequences. I<Annals of Mathematical Statistics>, I<28>, 262-264.

Good, I. J., & Gover, T. N. (1967). The generalized serial test and the binary expansion of [square-root]2. I<Journal of the Royal Statistical Society A>, I<130>, 102-107.

Kendall, M. G., & Babington Smith, B. (1938). Randomness and random sampling numbers. I<Journal of the Royal Statistical Society>, I<101>, 147-166.

Knuth, D. E. (1998). I<The art of computer programming> (3rd ed., Vol. 2 Seminumerical algorithms). Reading, MA, US: Addison-Wesley.

Rukhin, A., Soto, J., Nechvatal, J., Smid, M., Barker, E., Leigh, S., et al. (2001). A statistical test suite for random and pseudorandom number generators for cryptographic applications. Retrieved September 4 2010, from L<http://csrc.nist.gov/groups/ST/toolkit/rng/documents/SP800-22b.pdf>.

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> for other tests of sequences, and for sharing data between these tests.

=head1 TO DO/BUGS

Implementation of the serial test for non-overlapping I<v>-nomes.

=head1 REVISION HISTORY

v.051: Default test-statistic reported in "dumps" is now the I<Z>-value associated with the psi-square I<p>-value; whatever the manner of calculating psi-square (by differencing or not), it can itself be retrieved as $vnomes->{'psisq'}, and is also dumped, by name, if C<text> => 2 in a call to C<dump>.

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

This ends documentation of the Perl implementation of the I<chi>-square, Kendall-Babington Smith, and Good's Generalized Serial Test for randomness in a sequence.

=cut

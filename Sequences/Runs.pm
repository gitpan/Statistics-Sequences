package Statistics::Sequences::Runs;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);

use Statistics::Deviation;
use Statistics::Sequences 0.01;
@ISA = qw(Statistics::Sequences);

$VERSION = '0.01';

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub test {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = Class::OOorNO::coerce_array(@_);
    
    $self->_check_testdata($args);

    my @data = @{$self->{'testdata'}};
    my ($obs_val, $i, %events) = ();
    
    # Count the number of runs:
    foreach ($i = 0; $i < scalar(@data); $i++) {
        ##print "$i = $data[$i]\n";next;
        $events{$data[$i]}++;
        $obs_val++ if !$i || ( $data[$i] ne $data[$i - 1] );
    }

    #croak __PACKAGE__, '::test Less than two events were found in the data' if scalar keys %events < 2;
    return undef if ! scalar keys %events;
    ##croak __PACKAGE__, '::test No events/no test data appear to be loaded' if ! scalar keys %events;
    croak __PACKAGE__, '::test More than two events were found in the data' if scalar keys %events > 2;

    # Now get the probability of these many runs:
    my ($m, $n) = keys %events;
    $n ||= 0; # account for only 1 event having appeared (i.e, a run of 1)
    my ($n1, $n2) = ($events{$m}, $events{$n});
    $n2 ||= 0;
    my $sum = $n1 + $n2; # the sum of all observations
    my $precision = defined $args->{'p_precision'} ? $args->{'p_precision'} : $self->{'p_precision'};

    if ($sum > 1 && scalar keys %events > 1) { # for a valid test of deviation, there should be at least 2 events in 2 trials
        my ($exp_val, $var, $tails, $z, $pz, $obs_dev, $std_dev) = ();
        $exp_val = ( ( 2 * $n1 * $n2 ) /  $sum ) + 1;              
        $var =   ( 2 * $n1 * $n2 * ( ( 2 * $n1 * $n2 ) - $sum) ) 
                /
                ( ( $sum**2 ) * ( $sum - 1 ) );
        
        #$self->{'dist'} = ($n1 <= 20 && $n2 <= 20) ? 'poisson' : 'normal';
        #$args->{'dist'} ||= $self->{'dist'};
        $args->{$_} ||= $self->{$_} foreach qw/tails ccorr/;

        my $ccorr = $sum <= 60 && $args->{'ccorr'} ? 1 : 0;
        my $dev = Statistics::Deviation->new(ccorr => $ccorr, tails => $args->{'tails'}, p_precision => $precision, distribution => $args->{'dist'});
        ($z, $pz, $obs_dev, $std_dev) = $dev->test(observed => $obs_val, expected => $exp_val, variance => $var);
        $self->{'expected'} = $exp_val;
        $self->{'z_value'} = $z;
        $self->{'p_value'} = $pz;
        $self->{'r_value'} = $dev->z_2_r($z);
        $self->{'obs_dev'} = $obs_dev;
        $self->{'std_dev'} = $std_dev;
        $self->{'variance'} = $var;
    }
    else {
        $self->{$_} = 0 foreach qw/observed z_value r_value obs_dev std_dev variance/;
        $self->{'p_value'} = $self->{'p_value'} = sprintf('%.' . $precision . 'f', 1);
        $self->{'expected'} = 1;
    }
    $self->{'observed'} = $obs_val;

    return $self;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub dump {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : Class::OOorNO::coerce_array(@_);
    $args->{'testname'} = 'Runs';
    $self->SUPER::_dump_pass($args);
    return $self;
}

__END__

=pod

=head1 NAME

Statistics::Runs - The Runs-test (Wald-Walfowitz or Swed-Eisenhard Test).

=head1 VERSION

This is documentation for version 0.01 of Statistics::Sequences::Runs, released 31 July 2006.

=head1 SYNOPSIS

 use Statistics::Sequences::Runs;
 $runs = Statistics::Sequences::Runs->new();
 $runs->load(qw/1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 1 1 0 1/);
 $runs->test()->dump();

=head1 DESCRIPTION

The Runs-test nonparametrically assesses the difference between two independent samples, or a difference within a single sample of dichotomous observations.

A run is a sequence of identical events on 1 or more consecutive trials. For example, in a signal-detection test, there'll be a series, over time, of hits (H) and misses (M), which might look like H-H-M-H-M-M-M-M-H. Here, there are 5 runs: 3 of hits, and 2 of misses. This number of runs can be compared with the number expected to occur by chance, given the number of observed hits and misses. More runs than expected generally indicates irregularity, or instability; fewer runs than expected indicates regularity, or stability.

=head1 METHODS

Methods are essentially as described in L<Statistics::Sequences>. See this manpage for how to handle non-dichotomous data, e.g., numerical data, or those with more than two categories; the relevant methods are not described here.

=head2 new

 $run = Statistics::Sequences::Runs->new();

Returns a new Runs object. Expects/accepts no arguments but the classname.

=head2 load

 $runs->load(@data);
 $runs->load(\@data);
 $runs->load('sample1' => \@data1, 'sample2' => \@data2)
 $runs->load({'sample1' => \@data1, 'sample2' => \@data2})

Loads data anonymously or by name. See L<load|Statistics::Sequences/load> in the Statistics::Sequences manpage.

=head2 test

 $runs->test();

Performs the runs test on the named samples. If only one sample name is given, the one-sample Runs test is performed, cutting the data at the median, or by the value given as cut. Observations that fall above and below the cut-point then constitute the "groups" to be searched for runs. Otherwise, with two named groups, runs are sought on the basis of the observations belonging to one or the other named group.

=head2 dump

 $runs->dump(data => '1|0', flag => '1|0', text => '0|1|2');

Print Runs-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details.

=head1 EXAMPLE

=head2 ESP runs

In a run of a classic ESP test, there are 25 trials where the responses are either hits (1) or misses (0) with respect to randomly selected targets, which can be any of 5 events (typically, geometric figures). Classic studies test the significance of the number of hits versus what is theoretically expected (usually, 5 hits in 25). Here we test the significance of the number of runs of hits (and misses) given the number actually produced, firstly by L<dichotomising the data|Statistics::Sequences/Dichotomising data> (which are based on 5, not the required 2, categories).

 use Statistics::Sequences::Runs;

 # Produce pseudo targets and responses:
 my ($i, @targets, @responses);
 for ($i = 0; $i < 25; $i++) {
    $targets[$i] = (qw/circle plus square star wave/)[int(rand(5))];
    $responses[$i] = (qw/circle plus square star wave/)[int(rand(5))];
 }

 # Do the run thing:
 my $runs = Statistics::Sequences::Runs->new();
 $runs->load(targets => \@targets, responses => \@responses);
 $runs->match(data => [qw/targets responses/]);
 $runs->test();
 print "The probability of obtaining these $runs->{'observed'} runs is $runs->{'p_value'}\n";

 # But what if the responses were actually guessed for the target on the trial one ahead?
 $runs->match(data => [qw/targets responses/], lag => 1)->test();
 print "With precognitive responses to the target for the next trial,\n 
 $runs->{'observed'} runs in 24 trials were produced when $runs->{'expected'} were expected,\n 
 for which the probability is $runs->{'p_value'}\n"; 

=head1 SEE ALSO

L<Statistics::Sequences|Statistics::Sequences> for other tests of sequences, and for sharing data between these tests.

=head1 AUTHOR/LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

Copyright (C) 2006 Roderick Garton  
This module may be modified, used, copied, and redistributed at your own risk.
Publicly redistributed modified versions must use a different name.

=head1 DISCLAIMER

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=cut

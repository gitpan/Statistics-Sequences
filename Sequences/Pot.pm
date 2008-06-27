package Statistics::Sequences::Pot;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);

use Statistics::Zed 0.01;
use Statistics::Sequences 0.02;
@ISA = qw(Statistics::Sequences);

$VERSION = '0.02';

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#sub new {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#    my $class = shift;
#    push @_, qw(samples states observed expected range obs_dev std_dev z p_value);
#    bless {}, $class;
#    $class->SUPER::new(@_);
#}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub test {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = Class::OOorNO::coerce_array(@_);
    
    $self->_check_testdata($args);
            
    my $m = scalar @{$self->{'testdata'}} || return undef; ##croak __PACKAGE__, '::test More than an empty list of data for pot-testing is needed';

    # Set the target-state, allowing for a zero value:
    my $state = defined $args->{'state'} ? 
        delete $args->{'state'} : 
        croak __PACKAGE__, '::test A state for pot-testing is needed';

	# Init an array holding the indices at which the state appears in the given data:
    # Meanwhile, build arrays of bunch and space frequencies, should this be requested:
    my ($k, $j, $i, @indices, @bunches, @spaces) = (0, 0);
    
    for ($i = 0; $i < $m; $i++) {
        # Allow for matching numerical or string values:
        if ($self->{'testdata'}->[$i] eq $state) {
             $k++ if $spaces[$k];
             $j++ if ( scalar @indices ) and ( $indices[-1] != ($i - 1) );
             $bunches[$j]++;
             push @indices, $i;
        }
        else {
            $spaces[$k]++;
        }
    }

    my $n = scalar @indices;
    # Assume scale = 1 if not specified or invalid
    $self->{'scale'} = (!$args->{'scale'} or $args->{'scale'} < 1) ? 1 : delete $args->{'scale'};

    # Lump relevant quantities into the Pot object:
    $self->{'state'} = $state;
    $self->{'samples'} = $m;
    $self->{'count'} = $n;
    $self->{'range'} = $n ? ( $self->{'scale'} * $m / $n ) : '';
    #$self->{'range_e'} = $r;

	# Provide optional descriptives for 'bunches' and their 'spaces':	
    if ($args->{'full'}) {
        require Statistics::Descriptive;
        
        $self->{'bunches'} = Statistics::Descriptive::Full->new();
        $self->{'bunches'}->add_data(@bunches);
       
        $self->{'spaces'} = Statistics::Descriptive::Full->new();
        $self->{'spaces'}->add_data(@spaces);
    }
	else {
	    $self->{'bunches'} = undef;
		$self->{'spaces'} = undef;
	}

    my $precision = defined $args->{'p_precision'} ? $args->{'p_precision'} : $self->{'p_precision'};
    # Resume bad-mouthing, if possible:
    if (!$n) { # State $state did not occur within the data-set for pot-testing
        $self->{$_} = 0 foreach qw/observed expected z_value obs_dev std_dev variance/;
        $self->{'p_value'} = sprintf('%.' . $precision . 'f', 1);
    }
	elsif ($n > $m) {# Numerator sd must != 0, otherwise deviation will soon be divided by it:
	    croak __PACKAGE__, "::test State $state occurred an improbable number of times: $n among $m observations";
	}
    else { # Do some proper work:
        my ($obs_val, $exp_val, $var) = _calc($n, $m, $self->{'scale'}, \@indices);
        $args->{'ccorr'} ||= $self->{'ccorr'};
        my $ccorr = $m <= 60 && $args->{'ccorr'} ? 1 : 0;
        my $dev = Statistics::Zed->new(ccorr => $ccorr, tails => $args->{'tails'}, distribution => $args->{'dist'}, p_precision => $precision);
        my ($z, $pz, $obs_dev, $std_dev) = $dev->zscore(observed => $obs_val, expected => $exp_val, variance => $var);
        # More lumping:
        $self->{'observed'} = $obs_val;
        $self->{'expected'} = $exp_val;
        $self->{'z_value'} = $z;
        $self->{'p_value'} = $pz;
        $self->{'r_value'} = $dev->z_2_r($z);
        $self->{'obs_dev'} = $obs_dev;
        $self->{'std_dev'} = $std_dev;
        $self->{'variance'} = $var;
    }

    # Knock-off time, go home, but allow for sudden calls to other methods:
    return $self;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub dump {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : Class::OOorNO::coerce_array(@_);

    if ($args->{'text'} and $args->{'text'} > 1) {
        print '-' x 50 . "\n";
        print "Pot test results for state $self->{'state'} among $self->{'samples'} observations:\n";
        print '-' x 50 . "\n";
        print " No. of observations of state = $self->{'count'}\n";
        print " No. of bunches of state = " . $self->{'bunches'}->count() . 
                  ', with a mean length of '. sprintf('%.2f', $self->{'bunches'}->mean()) .
                  ', and a mean spacing of '. sprintf('%.2f', $self->{'spaces'}->mean()) ." between each bunch.\n" if $self->{'bunches'};
        print ' Pot calculated with a range of ' . sprintf('%.2f', $self->{'range'}) . " over a scale of $self->{'scale'}\n";   
        $args->{'testname'} = 'Pot';
        $self->SUPER::_dump_verbose($args);
     }
     else {
        $args->{'testname'} = "Pot($self->{state})";
        $self->SUPER::_dump_sparse($args);
    }
    return $self;
}

#-----------------------------------------------
sub _calc {
#-----------------------------------------------
    my ($n, $m, $scale, $indices, $obs_val, $exp_val, $var, $i, $j, $r) = (@_[0 .. 3], 0, 0, 0);
 
    # Init range parameter:
    $r = exp(-$n / $m * $scale);

    # Calculate observed Pot: Schmidt (2000) Equations 6-7:
    for $i (1 .. ($n - 1)) {
        for $j (0 .. $i - 1) {
            $obs_val += $r**abs($indices->[$i] - $indices->[$j]);
        }
    }

    # Calculate expected Pot: Schmidt (2000) Equation 8:
    if ($m > 1) {
        $exp_val = $n * ($n - 1) * $r * ( $m - 1 / (1 - $r) ) / ( $m * ($m - 1) * (1 - $r) );
    }
    
	# Calculate variance: Schmidt (2000) Equation 9a:
    if ($m) {
        $var = ( $r**2 * $n**2 * (1 - $n / $m)**2 )
              /
           ( $m * (1 - $r**2) );
    }

    return ($obs_val, $exp_val, $var);
}

1;

__END__

=pod

=head1 NAME

Statistics::Sequences::Pot - Helmut Schmidt's test of force-like runs of a discrete state within a random distribution

=head1 VERSION

This is documentation for version 0.02 of Statistics::Sequences::Pot, released 27 June 2008.

=head1 SYNOPSIS

 use Statistics::Sequences::Pot;

 $pot = Statistics::Sequences::Pot->new();

 # Load an array (reference) of data (of strings or numbers) into the pot object:
 $pot->load([qw/2 0 8 5 3 5 2 3 1 1 9 4 4 1 5 5 6 5 8 7 5 3 8 5 6/]);

 # Test the relative runs of a specific state (e.g., "5") among these data:
 $pot->test(state => 5);

 # Print out the Pot-statistic, and a z-test of its significance:
 $pot->dump();

 # Prints: State 5: Pot = 4.31, z = -0.19, p = 0.42539

 # or be discretely informed re individual stats, post-test, e.g., :

 print "Observed Pot = $pot->{'observed'} for $pot->{'count'} occurrences
  of $pot->{'state'} among $pot->{'samples'};
  probability = $pot->{'p_value'}\n";

 # Prints: Observed Pot = 4.3098698035002 for 7 occurrences of 5 among 25; probability = 0.42539

=head1 DESCRIPTION

The Pot statistic measures the bunching relative to the spacing of a single state within a series of other states, conceived by Helmut Schmidt as a targeted "potential" energy (or Pot) that dissipates exponentially between states. It's not limited to considering only clusters of I<consecutive> states (or bunches), as is the case with the more familiar Runs test of sequences.

Say you're interested in the occurrence of the state B<3> within an array of digits: note how, in the following arrays, there are increasing breaks between the B<3>s (separated by 0, 1 and then 2 other states):

 4, 7, 3, 3
 3, 4, 3, 7
 3, 8, 1, 3

The occurrence of B<3> is, with the Pot-test, of exponentially declining interest across these sequences, given the increasing breaks by other states between the occurrences of 3. The statistic does not ignore these ever remoter occurrences of the state of interest; it accounts for increased spacing between them as if there were an exponentially declining force, a I<pot>ential towards B<3>, within the data-stream (up to a theoretical or empirical asymptote that may be specified).

Running the Pot-test involves testing its significance as a standard "z" score; Schmidt (2000) provided data demonstrating Pot's conformance with the normal distribution. This will, of course, be improved by repeated sampling, and by pooling observations into blocks.

=head1 METHODS

Methods are essentially as described in L<Statistics::Sequences>. See this manpage for how to handle non-dichotomous data, e.g., numerical data, or those with more than two categories; the relevant methods are not described here.

=head2 new

 $pot = Statistics::Sequences::Runs->new();

Returns a new Runs object. Expects/accepts no arguments but the classname.

=head2 load

 $pot->load(@data);
 $pot->load(\@data);
 $pot->load('sample1' => \@data1, 'sample2' => \@data2)
 $pot->load({'sample1' => \@data1, 'sample2' => \@data2})

Loads data anonymously or by name. See L<load|Statistics::Sequences/load> in the Statistics::Sequences manpage.

=head2 test

 $pot->test(state => 'a'[, scale => 1, full => 0]);

I<Aliases:> C<perform_pot_test>

Runs the Pot-test on a specified state, lumps the Pot object with stats, and returns itself. Note that data must already be loaded into the Pot object before testing, otherwise, expect a C<croak>. The test works with the following required, and then some optional, parameters, each as C<name =E<gt> value> pairs. 

=over 4

=item B<state> => I<string>

The state within the data whose bunching is to be tested. This is the only required parameter to L<test|test>, which will surely C<croak> if no state is specified, or if a state occurs more often than there are data. If the state does not exist in the data, nix is defined.

=item B<scale> => I<numeric> E<gt>= 1

Optionally, the scale of the range parameter, which should be greater than or equal to 1. Default = 1; values less than 1 are effected as 1. For info on how to set this parameter, see the L</DESCRIPTION> above, and also the explanation of C<observed> among the statistical L</ATTRIBUTES>.

=item B<full> => I<boolean>

Optionally, standard descriptive statistics regarding the bunches, and the spaces between them, will be available. See C<bunches> and C<spaces>, below. Default = 0, given the independence of Pot from the central tendencies of bunching and spacing. The facility is provided for purposes of hypothesis-testing, and appraising the character of Pot.

=back

=head2 dump

 $pot->dump(flag => '1|0', text => '0|1|2');

Print Pot-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details.

=head1 ATTRIBUTES

Once calling L<test|test>, the pot object is lumped with the following attributes, each of which may be accessed as C<$pot-E<gt>{'ATTRIBUTE'}>.

=over

=item B<samples>

The size of the submitted data (or number of "trials").

=item B<count>

The number of times the target state occurred in the submitted array.

=item B<observed>

A measure of the number and size of bunchings of the state that occurred within the array.
This is based on Schmidt (2000), Equations 6-7, and his Appended program. The formula is:

=for html <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>I</i>,<i>J</i>=1..<i>N</i><br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;SUM&nbsp;&nbsp;<i>r</i><sup>|<i>n</i>(<i>I</i>) - <i>n</i>(<i>J</i>)|</sup><br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>I</i>&lt;<i>J</i></p>

=item

where 

=for html <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>r</i> = <i>e</i><sup>-<i>N</i>/<i>MS</i></sup></p>

=item

is the number of observations, and I<S> (for scale) is a constant determining the range I<r> of the potential energy between pairs of I<I> and I<J> states.

In most situations, should all states be equiprobable, or their probability be proportionate to their number, I<r> would reflect the average distance, or delay, between I<successive> states, equal to the number of all observations divided by the number of states. For example, if there were 10 possible states, and 100 observations have been made, then the probability of re-occurrence of any one of the 10 states within any slot will be equal to 100/10, with I<S> = 1, i.e., expecting that any one of the states would mostly occur by a spacing of 10, and then by an exponentially declining tendency toward consecutive occurrence. In this way, with I<S> = 1, Pot can be considered to be a measure of "short-range bunching," as Schmidt called it. Bunching over a larger range than this minimally expected range can be measured with I<S> > 1. This is specified, optionally, as the argument named I<scale> to L<test|test>. Hypothesis-testing might be made with respect to various values of the I<scale> parameter.

=item B<expected>

The theoretically expected value of Pot, given I<N> states among I<M> observations, and I<r> range of clustering within these observations. It is calculated as follows, given the above definitions.

=for html <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pot = ((<i>N</i>(<i>N</i> - 1))/(<i>M</i>(<i>M</i> - 1))) . (<i>r</i>/1 - <i>r</i>) . (<i>M</i> - (1/(1 - <i>r</i>)))</p>

=item B<obs_dev>

The observed value of Pot less the expected value of Pot, and hence the amount by which the observed value deviates from that expected by chance.

=item B<std_dev>

The standard deviation, based on the theoretically expected variance of Pot, which is given by:

=for html <p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Variance = (<i>r</i><sup>2</sup>/ (1 - <i>r</i><sup>2</sup>) . (<i>N</i> / <i>M</i>) . (1 - (<i>N</i> / <i>M</i>))<sup>2</sup> . <i>N</i></p>

=item B<range>

The range of observations over which Pot was assessed, simply being the product of the I<scale> and number of observations, divided by the number of states.

=item B<z_value>

The standard score, based on dividing the observed deviation by the standard deviation.

=item B<p_value>

The probability associated with the absolute value of the I<z>-score.

=item B<bunches>

I<Only provided if passing> C<full =E<gt> 1> I<to> L<test|test>. 

A L<Statistics::Descriptive::Full|Statistics::Descriptive> object, loaded with the lengths of each bunch of the state. Statistics such as the I<count>, I<mean>, I<mode> and I<range> of the observed bunches can be called, and the ordered list of bunch sizes can itself be retrieved. E.g.,

 print $pot->{'bunches'}->mean() ."\n";
 @bunch_lengths = $pot->{'bunches'}->get_data();

=item B<spaces>

I<Only provided if passing> C<full =E<gt>> 1 I<to> L<test|test>.

A L<Statistics::Descriptive::Full|Statistics::Descriptive> object, loaded with the lengths of the intervals, or spaces, between each bunch of the state. Statistics such as the I<count>, I<mean>, I<mode> and I<range> of the observed spaces can be called. E.g., 

 print $pot->{'spaces'}->mode() ."\n";

=back

Note that L<test|test> returns the Pot object (itself), so you could get "immediate" access to any of the above by, for example:

 my $n_zeroes = $pot->test(state => 0)->{'count'};

 print 'median spaces = ' . $pot->test(state => 0, full => 1)->{'spaces'}->median() . "\n";

=head1 EXAMPLES 

1. Using Pot as a test of bunching of a particular state within a collection of quasi-random observations.

 use Statistics::Sequences::Pot;
 use strict;
 my ($i, @data) = ();

 # Init an array of random data with integers ranging from 0 to 15:
 for ($i = 0; $i < 960; $i++) {
   $data[$i] = int(rand(16));
 }

 # Assess degree of bunching within these data with respect to a randomly selected target state:
 my $state = int(rand(16));

 my $pot = Statistics::Sequences::Pot->new();
 $pot->load(\@data)->test(state => $state);

 # Access the results of this analysis:
 print "The probability of obtaining as much bunching of $state as observed is $pot->{'p_value'}\n";
 # or:
 print "For state $pot->{'state'} occurring $pot->{'count'} times among $pot->{'samples'}, ".
    "the observed value of Pot was $pot->{'observed'}.\n" .
    "The expected value of Pot was $pot->{'expected'}\n" .
    "The deviation from expectation was $pot->{'obs_dev'}\n" .
    "The standard deviation was $pot->{'std_dev'}\n" .
    "Z = $pot->{'z_value'}\n" .
    "Probability of this z-value = $pot->{'p_value'}\n";
 # or print the lot, and more, in English:
 $pot->dump(text => 2);

 # See what else was happening, having already given test() the data to test:
 foreach (0 .. 15) {
    next if $_ == $state;
    $pot->test(state => $_)->dump();
  }

2. Using Pot as a test of randomness of an array of dichotomous observations. Note: alphabetic strings as the elements of the array; reuse of loaded data; recycling of loaded data between tests; internal storage of the state; and exploitation of the module for semi-Pot purposes.

  use Statistics::Sequences::Pot;
  use strict;
  my ($i, @data) = ();

  # Init an array of random data with values of either 'hit' or 'miss':
  my @categories = (qw/hit miss/);  
  for ($i = 0; $i < 640; $i++) {
    $data[$i] = $categories[int(rand(@categories))];
  }

  # Make a pot object and load up the data:
  my $pot = Statistics::Sequences::Pot->new();
  $pot->load(\@data);

  # Run and dump a couple analyses, on each possible state:
  $pot->test(state => 'hit');
  $pot->dump();
  $pot->test(state => 'miss');
  $pot->dump();

  # Be randomly redundant:
  $pot->test(state => $categories[int(rand(@categories))], full => 1);

  print "Randomly selected state was a $pot->{'state'}, and this occurred $pot->{'count'} times, " .
        "most frequently bunching by a length of " . $pot->{'bunches'}->mode() . "\n";

  # Prints, e.g.:
  ## State hit: Pot = 252.04, z = 0.83, p = 0.20298
  ## State miss: Pot = 246.43, z = 0.82, p = 0.20721
  ## Randomly selected state was a miss, and this occurred 315 times, most frequently bunching by a length of 1

=head1 REFERENCES

Schmidt, H. (2000). A proposed measure for psi-induced bunching of randomly spaced events. I<Journal of Parapsychology, 64,> 301-316.

=head1 SEE ALSO

L<http://www.fourmilab.ch/rpkp/> for Schmidt's many papers on the physical conceptualisation and properties of psi.

L<Statistics::Descriptive|Statistics::Descriptive> : The present module adds data to "Full" objects of this package in order to access descriptives re bunches and spaces.

L<Statistics::Distributions|Statistics::Distributions> : The present module uses the C<uprob()> method of this package for determining the probability associated with the I<z>-value.

L<Statistics::Frequency|Statistics::Frequency> : the C<proportional_frequency()> method in this module could be informative when working with data of the kind used here.

=head1 BUGS/LIMITATIONS

No computational bugs as yet identfied. Hopefully this will change, given time.

Limitations of the code, perhaps, concern the non-unique storage of data arrays (compared to, say, C<Statistics::DependantTTest>, but see C<Statistics::TTest>). This would require a unique name for each array of data, and explicit reference to one or another array with each L<test|test> (when, perhaps, you'd have only one data-set, after all). In any case, the data are accepted as array references.

Limitations of the actual Pot statistic may be considered to be its newness, not having had the opportunity to be critiqued by peers, and that some experimentation might be required to find an optimal C<scale>.

=head1 REVISION HISTORY

=over 4

=item v 0.02

June 2008

Initital release via PAUSE.

See CHANGES in installation dist for revisions.

=back

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2008 Roderick Garton

rgarton@utas_DOT_edu_DOT_au

This program is free software. This module is free software. It may be used, redistributed and/or modified under the stame terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut

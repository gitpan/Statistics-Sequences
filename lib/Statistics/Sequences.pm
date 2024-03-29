package Statistics::Sequences;
use strict;
use warnings FATAL => 'all';
use Carp qw(croak cluck);
use Statistics::Data 0.08;
use base qw(Statistics::Data);
use Scalar::Util qw(looks_like_number);
our $VERSION = '0.12';

=pod

=head1 NAME

Statistics::Sequences - Manage sequences (ordered list of literals) for testing their runs, joins, turns, trinomes, potential energy, etc.

=head1 VERSION

This is documentation for Version 0.12 of Statistics::Sequences.

=head1 SYNOPSIS

 use Statistics::Sequences 0.12;
 $seq = Statistics::Sequences->new();
 my @data = (1, 'a', 'a', 1); # ordered list of literal scalars (numbers, strings), as permitted by specific test
 $seq->load(\@data); # or @data or dataname => \@data
 print $seq->observed(stat => 'runs'); # expected, variance, z_value, p_value - assuming sub-module Runs.pm is installed
 print $seq->test(stat => 'vnomes', length => 2); # - - assuming sub-module Vnomes.pm is installed
 $seq->dump(stat => 'runs', values => {observed => 1, z_value => 1, p_value => 1}, exact => 1, tails => 1);
 # see also Statistics::Data for inherited methods

=head1 DESCRIPTION

Loading, updating and accessing data as ordered list of literal scalars (numbers, strings) for statistical tests of their sequential structure via L<Statistics::Sequences::Joins|Statistics::Sequences::Joins>, L<Statistics::Sequences::Pot|Statistics::Sequences::Pot>, L<Statistics::Sequences::Runs|Statistics::Sequences::Runs>, L<Statistics::Sequences::Turns|Statistics::Sequences::Turns> and L<Statistics::Sequences::Vnomes|Statistics::Sequences::Vnomes>. Note that none of these sub-modules are installed by default; to use this module as intended, install one or more of these sub-modules.

To access the tests, L<use|perlfunc/use> this base module to create a Statistics::Sequences object with L<new|new>, then L<load|load> data into it and access each test by calling the L<test|test> method, specifying the B<stat> attribute: either joins, pot, runs, turns or vnomes, where the relevant sub-module is installed. This allows running several tests on the same data, as the data are immediately available to each test (of joins, pot, runs, turns or vnomes). See the L<SYNOPSIS|Statistics::Sequences/SYNOPSIS> for a simple example. 

Alternatively, L<use|perlfunc/use> each sub-module directly, and restrict analyses to the sub-module's test; this module is used implicitly as their base. That is, to perform a test of one type (e.g., runs), L<use|perlfunc/use> the relevant sub-package, load data via its constructor; see the SYNOPSIS for the particular test, i.e., L<Joins|Statistics::Sequences::Joins/SYNOPSIS>, L<Pot|Statistics::Sequences::Pot/SYNOPSIS>, L<Runs|Statistics::Sequences::Runs/SYNOPSIS>, L<Turns|Statistics::Sequences::Turns/SYNOPSIS> or L<Vnomes|Statistics::Sequences::Vnomes/SYNOPSIS>. You won't be able to access other tests of the same data by this approach, unless you create another object for that test, and then specifically pass the data from the earlier object into the new one.

=head1 SUBROUTINES/METHODS

=head2 new

 $seq = Statistics::Sequences->new();

Returns a new Statistics::Sequences object (inherited from L<Statistics::Data|Statistics::Data>) by which all the methods for caching, reading and testing data can be accessed, including each of the methods for performing the L<Runs-|Statistics::Sequences::Runs>, L<Joins-|Statistics::Sequences::Joins>, L<Pot-|Statistics::Sequences::Pot>, L<Turns-|Statistics::Sequences::Turns> or L<Vnomes-|Statistics::Sequences::Vnomes>tests.

Sub-packages also have their own new method - so, e.g., L<Statistics::Sequences::Runs|Statistics::Sequences::Runs>, can be individually imported, and its own L<new|new> method can be called, e.g.:

 use Statistics::Sequences::Runs;
 $runs = Statistics::Sequences::Runs->new();

In this case, data are not automatically shared across packages, and only one test (in this case, the Runs-test) can be accessed through the class-object returned by L<new|new>.

=head2 load, add, access, unload

All these operations on the basic data are inherited from L<Statistics::Data|Statistics::Data> - see this doc for details of these and other possible methods.

B<Dichotomous data>: Both the runs- and joins-tests expect dichotomous data: a binary or binomial or Bernoulli sequence, but with whatever characters to symbolize the two possible events. They test their "loads" to make sure the data are dichotomous. To reduce numerical and categorical data to a dichotomous level, see the L<pool|Statistics::Data::Dichotomize/pool>, L<match|Statistics::Data::Dichotomize/match>, L<split|Statistics::Data::Dichotomize/split, cut>, L<swing|Statistics::Data::Dichotomize/swing>, L<shrink (boolwin)|Statistics::Data::Dichotomize/shrink, boolwin> and other methods in L<Statistics::Data::Dichotomize|Statistics::Data::Dichotomize>.

=head2 observed, observation

 $v = $seq->observed(stat => 'joins|pot|runs|turns|vnomes', %args); # gets data from cache, with any args needed by the stat
 $v = $seq->observed(stat => 'joins|pot|runs|turns|vnomes', data => [qw/blah bing blah blah blah/]); # just needs args for partic.stats
 $v = $seq->observed(stat => 'joins|pot|runs|turns|vnomes', label => 'myLabelledLoadedData'); # just needs args for partic.stats

Return the observed value of the statistic for the L<load|Statistics::Sequences/load>ed data, or data sent with this call, eg., how many runs in the sequence (1, 1, 0, 1). See the particular statistic's manpage for any other arguments needed or optional. 

=cut

sub observed { return _feedme('observed', @_); } *observation = \&observed;

=head2 expected, expectation

 $v = $seq->expected(stat => 'joins|pot|runs|turns|vnomes', %args); # gets data from cache, with any args needed by the stat
 $v = $seq->expected(stat => 'joins|pot|runs|turns|vnomes', data => [qw/blah bing blah blah blah/]); # just needs args for partic.stats

Return the expected value of the statistic for the L<load|Statistics::Sequences/load>ed data, or data sent with this call, eg., how many runs should occur in a 4-length sequence of two possible events. See the statistic's manpage for any other arguments needed or optional.

=cut

sub expected { return _feedme('expected', @_); } *expectation = \&expected;

=head2 variance

 $seq->variance(stat => 'joins|pot|runs|turns|vnomes', %args); # gets data from cache, with any args needed by the stat
 $seq->variance(stat => 'joins|pot|runs|turns|vnomes', data => [qw/blah bing blah blah blah/]); # just needs args for partic.stats

Returns the expected range of deviation in the statistic's observed value for the given number of trials. 

=cut

sub variance { return _feedme('variance', @_); }

=head2 obsdev, observed_deviation

 $v = $seq->obsdev(stat => 'joins|pot|runs|turns|vnomes', %args); # gets data from cache, with any args needed by the stat
 $v = $seq->obsdev(stat => 'joins|pot|runs|turns|vnomes', data => [qw/blah bing blah blah blah/]); # just needs args for partic.stats

Returns the deviation of (difference between) observed and expected values of the statistic for the loaded/given sequence (I<O> - I<E>). 

=cut

sub obsdev {
    return observed(@_) - expected(@_);
}
*observed_deviation = \&obsdev;

=head2 stdev, standard_deviation

 $v = $seq->stdev(stat => 'joins|pot|runs|turns|vnomes', %args); # gets data from cache, with any args needed by the stat
 $v = $seq->stdev(stat => 'joins|pot|runs|turns|vnomes', data => [qw/blah bing blah blah blah/]); # just needs args for partic.stats

Returns square-root of the variance.

=cut

sub stdev {
    return sqrt variance(@_);
}
*standard_deviation = \&stdev;

=head2 z_value, zscore

 $v = $seq->zscore(stat => 'joins|pot|runs|turns|vnomes', %args); # gets data from cache, with any args needed by the stat
 $v = $seq->zscore(stat => 'joins|pot|runs|turns|vnomes', data => [qw/blah bing blah blah blah/]); # just needs args for partic.stats

Return the deviation ratio: observed deviation to standard deviation. Use argument B<ccorr> for continuity correction.

=cut

sub zscore { return _feedme('zscore', @_); } *z_value = \&zscore;

=head2 p_value, test

 $p = $seq->test(stat => 'runs');
 $p = $seq->test(stat => 'joins');
 $p = $seq->test(stat => 'turns');
 $p = $seq->test(stat => 'pot', state => 'a value appearing in the data');
 $p = $seq->test(stat => 'vnomes', length => 'an integer greater than zero and less than sample-size');

Returns the probability of observing so many runs, joins, etc., versus those expected, relative to the expected variance.

When using a Statistics::Sequences class-object, this method requires naming which test to perform, i.e., runs, joins, pot or vnomes. This is I<not> required when the class-object already refers to one of the sub-modules, as created by the C<new> method within L<Statistics::Sequences::Runs|Statistics::Sequences::Runs/new>, L<Statistics::Sequences::Joins|Statistics::Sequences::Joins/new>, L<Statistics::Sequences::Pot|Statistics::Sequences::Pot/new>, L<Statistics::Sequences::Turns|Statistics::Sequences::Turns/new> and L<Statistics::Sequences::Vnomes|Statistics::Sequences::Vnomes/new>.

=head3 Common options

Options common to all the sub-package tests are as follows.

=over 8

=item data => 'I<string>'

Optionally specify the name of the data to be tested. By default, this is not required: the data tested are those that were last loaded, either anonymously, or as returned by one of the L<Statistics::Data::Dichotomize|Statistics::Data::Dichotomize> methods. Otherwise, I<if the data are already ready for testing in a dichotomous format>, data that were previously loaded by name can be individually tested. For example, here are two sets of data that are loaded by name, and then a single test of one of them is performed.

 @chimps = (qw/banana banana cheese banana cheese banana banana banana/);
 @mice = (qw/banana cheese cheese cheese cheese cheese cheese cheese/);
 $seq->load(chimps => \@chimps, mice => \@mice);
 $p = $seq->test(stat => 'runs', data => 'chimps');

=item ccorr => I<boolean>

Specify whether or not to perform the continuity-correction on the observed deviation. Default is false. Relevant only for those tests relying on a I<Z>-test. See L<Statistics::Zed|Statistics::Zed>.

=item tails => I<1>|I<2>

Specify whether the I<z>-value is calculated for both sides of the normal (or chi-square) distribution (2, the default for most tested data) or only one side (the default for data prepared with the B<swing> method.

=back

=head3 Test-specific required settings and options

Some sub-package tests need to have parameters defined in the call to L<test|test>, and/or have specific options, as follows.

B<Joins> : The Joins test I<optionally> allows the setting of a probability value; see C<test|test> in the  L<Statistics::Sequences::Joins|Statistics::Sequences::Joins/test> manpage.

B<Pot> : The Pot test I<requires> the setting of a state to be tested; see C<test> in the  L<Statistics::Sequences::Pot|Statistics::Sequences::Pot/test> manpage.

B<Vnomes> : The Serial test for v-nomes requires a length, i.e., the value of I<v>; see C<test> in the L<Statistics::Sequences::Vnomes|Statistics::Sequences::Vnomes/test> manpage..

B<Runs>, B<Turns> : There are presently no specific requirements nor options for the Runs- and Turns-tests. 

=cut

sub p_value { return _feedme('p_value', @_); } *test = \&p_value;

=head2 stats_hash

 $href = $seq->stats_hash(stat => 'runs', values => {observed => 1, expected => 1, variance => 1, z_value => 1, p_value => 1});

Returns a hashref with values for any of the descriptives and probability value relevant to the specified B<stat>istic. Include other required or optional arguments relevant to any of the values requested, e.g., B<ccorr> if getting a z_value, B<tails> and B<exact> if getting a p_value, B<state> if testing pot, B<prob> if testing joins, ... B<precision_s>, B<precision_p> ... 

=cut

sub stats_hash {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my @methods = keys %{$args->{'values'}};
    my (%stats_hash) = ();
    no strict 'refs';
    foreach my $method(@methods) {
        if ($args->{'values'}->{$method} == 1) {
            eval {$stats_hash{$method} = $self->$method($args);};
            croak "Method $method is not defined or correctly called for " . __PACKAGE__ if $@;
        }
    }
    if (! scalar keys %stats_hash) { # get default stats:
        foreach my $method(qw/observed p_value/) {
            eval {$stats_hash{$method} = $self->$method($args);};
            croak "Method $method is not defined or correctly called for " . __PACKAGE__ if $@;
        }
    }
    return \%stats_hash;
}

=head2 dump

 $seq->dump(stat => 'runs|joins|pot ...', values => {}, format => 'string|table', flag => '1|0', precision_s => 'integer', precision_p => 'integer');

I<Alias>: B<print_summary>

Print results of the last-conducted test to STDOUT. By default, if no parameters to C<dump> are passed, a single line of test statistics are printed. Options are as follows.

=over 8

=item values => hashref

Hashref of the statistical parameters to dump. Default is observed value and p-value for the given B<stat>.

=item flag => I<boolean>

If true, the I<p>-value associated with the I<z>-value is appended with a single asterisk if the value if below .05, and with two asterisks if it is below .01.

If false (default), nothing is appended to the I<p>-value.

=item format => 'table|labline|csv'

Default is 'csv', to print the stats hash as a comma-separated string (no newline), e.g., '4.0000,0.8596800". If specifying 'labline', you get something like "observed = 4.0000, p_value = 0.8596800\n". If specifying "table", this is a dump from L<Text::SimpleTable|Text::SimpleTable> with the stat methods as headers and column length set to the maximum required for the given headers, level of precision, flag, etc. For example, with B<precision_s> => 4 and B<precision_p> => 7, you get:

 .-----------+-----------.
 | observed  | p_value   |
 +-----------+-----------+
 | 4.0000    | 0.8596800 |
 '-----------+-----------'

=item verbose => 1|0

If true, includes a title giving the name of the statistic, details about the hypothesis tested (if B<p_value> => 1 in the B<values> hashref), et al. No effect if B<format> is not defined or equals 'csv'.

=item precision_s => 'I<non-negative integer>'

Precision of the statistic values (observed, expected, variance, z_value).

=item precision_p => 'I<non-negative integer>'

Specify rounding of the probability associated with the I<z>-value to so many digits. If zero or undefined, you get everything available.

=back

=cut

sub dump {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $stats_hash = $self->stats_hash($args);
    $args->{'format'} ||= 'csv';
    my @standard_methods = (qw/observed expected variance obsdev stdev z_value p_value/);
    my ($maxlen, @strs, @headers, @wanted_methods) = (0);
    foreach my $method(@standard_methods) { # set up what has been requested in a meaningful order:
        push(@wanted_methods, $method) if defined $stats_hash->{$method};
    }
    foreach my $method(keys %{$stats_hash}) {
        push(@wanted_methods, $method) if ! grep/$method/, @wanted_methods;
    }
    foreach my $method(@wanted_methods) {
        my $val = delete $stats_hash->{$method};
        my $len;
        if ($method eq 'p_value') {
            $val = _precisioned($args->{'precision_p'}, $val);
            $val .= ($val < .05 ? ($val < .01 ? q{**} : q{*} ) : q{}) if $args->{'flag'};
        }
        else {
            if (ref $val) {
                if (ref $val eq 'HASH') {
                    my %vals = %{$val};
                    $val = q{};
                    my $delim = $args->{'format'} eq 'table' ? "\n" : q{,};
                    my ($str, $this_len) = (q{});
                    while (my($k, $v) = each %vals) {
                        $str = "'$k' = $v";
                        $this_len = length($str);
                        $len = $this_len if not defined $len or $this_len > $len;
                        $val .= $str;
                        $val .= $delim;
                    }
                    if ($args->{'format'} ne 'table') {
                        chop $val;
                        $val = '(' . $val . ')';
                    }
                }
                else {
                    $val = join q{, }, @{$val};
                }
            }
            elsif(looks_like_number($val)) {
                $val = _precisioned($args->{'precision_s'}, $val);
            }
        }
        push @headers, $method;
        push(@strs, $val);
        $len = length $val if ! defined $len;
        $maxlen = $len if $len > $maxlen;
    }
    if ($args->{'format'} eq 'table') {
        $maxlen = 8 if $maxlen < 8;
        my $title = $args->{'verbose'} ? ucfirst($args->{'stat'}) . " statistics\n" : q{};
        print $title or croak 'Cannot print title for data-table';
        my @hh = ();
        push( @hh, [$maxlen, $_]) foreach @headers;
        require Text::SimpleTable;
        my $tbl = Text::SimpleTable->new(@hh);
        $tbl->row(@strs);
        print $tbl->draw or croak 'Cannot print data-table';
    }
    elsif ($args->{'format'} eq 'labline') {
        my @hh;
        for (my $i = 0; $i <= $#strs; $i++) {
            $hh[$i] = "$headers[$i] = $strs[$i]";
        }
        my $str = join(q{, }, @hh);
        if ($args->{'verbose'}) {
            $str = ucfirst($args->{'stat'}) . ': ' . $str;
        }
        print {*STDOUT} $str, "\n" or croak 'Cannot print data-string';
    }
    else { # csv
        print join(q{,}, @strs) or croak 'Cannot print data-string'
    }
    return;
}
*print_summary = \&dump;

=head2 dump_data

 $seq->dump_data(delim => "\n");

Prints to STDOUT a space-separated line of the tested data - as dichotomized and put to test. Optionally, give a value for B<delim> to specify how the datapoints should be separated. Inherited from L<Statistics::Data|Statistics::Data/dump_data>.

=cut

# PRIVATMETHODEN

sub _feedme {
    my $method = shift;
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $statname = $args->{'stat'} || q{};
    my $class = __PACKAGE__ . q{::} . ucfirst($statname);
    eval "require $class";
    croak __PACKAGE__, " error: Requested sequences module '$class' is not valid/available. You might need to install '$class'" if $@;
    my ($val, $nself) = (q{}, {});
    #my $nself = {};
    bless($nself, $class);#$nself = $class->new();
    $nself->{$_} = $self->{$_} foreach keys %{$self};
    no strict 'refs';
    eval '$val = $nself->$method($args)'; # but does not trap "deep recursion" if method not defined
    croak __PACKAGE__, " error: Method '$method' is not defined for $class" if $@;
    $self->{'stat'} = $statname;
    return $val;
}

sub _precisioned {
    return $_[0] ? sprintf(q{%.} . $_[0] . 'f', $_[1]) : (defined $_[1] ? $_[1] : q{}); # don't lose any zero
}

=head1 BUNDLING

This module C<use>s its sub-modules implicitly - so a bundled program using this module might need to explicitly C<use> its sub-modules if these need to be included in the bundle itself.

=head1 AUTHOR

Roderick Garton, C<< <rgarton at cpan.org> >>

=head1 LICENSE AND COPYRIGHT

=over 4

=item Copyright (c) 2006-2013 Roderick Garton

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut

1; # end of Statistics::Sequences

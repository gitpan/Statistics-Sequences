package Statistics::Sequences;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak cluck);
use vars qw($VERSION);
use Statistics::Descriptive;
use Statistics::Zed 0.02;
use Statistics::Lite qw(:funcs);
use Scalar::Util qw(looks_like_number);

$VERSION = '0.041';
$| = 1;

=pod

=head1 NAME

Statistics::Sequences - Tests of sequential structure in the form of runs, joins, bunches, etc.

=head1 SYNOPSIS

  use Statistics::Sequences;
  $seq = Statistics::Sequences->new();
  $seq->load([1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1]); # dichotomous values
  $seq->test(what => 'runs')->dump(); # or 1st argument to test is 'joins' or 'pot'
  # (prints:)
  # Z = -0.30, p = 0.76206

=head1 DESCRIPTION

Loading and preparing data for statistical tests of their sequential structure via L<Statistics::Sequences::Runs|Statistics::Sequences::Runs>, L<Statistics::Sequences::Joins|Statistics::Sequences::Joins>, and L<Statistics::Sequences::Pot|Statistics::Sequences::Pot>. Examples of the use of each test are given in these pages.

In general, to access the tests, you L<use|perlfunc/use> this base module to directly create a Statistics::Sequences object with the L<new|new> method, L<load|load> data into it, and then access each test by calling the L<test|test> method. This approach is useful for running several tests on the same data, as the data are immediately available to each test (of runs, pot and joins). See the L<SYNOPSIS|Statistics::Sequences/SYNOPSIS> for a simple example.

If you only want to perform a test of one type (e.g., runs), you might want to simply L<use> the relevant sub-package, create a class object specific to it, and load data specfically for its use; see the SYNOPSIS for the particular test, i.e., L<Runs|Statistics::Sequences::Runs/SYNOPSIS>, L<Joins|Statistics::Sequences::Joins/SYNOPSIS> or L<Pot|Statistics::Sequences::Pot/SYNOPSIS>. You won't be able to access other tests by this approach, unless you create another object for that test, and then specifically pass the data from the earlier object into the new one.

Note also that there are methods to anonymously or nominally cache data, and that data might need to be reduced to a dichotomous format, before a valid test can be run. Several dichotomising methods are provided, once data are loaded, and accessible via the generic or specific class objects, as above.

=head1 METHODS

=head2 Interface

The package provides an object-oriented interface for performing the L<Runs-|Statistics::Sequences::Runs>, L<Joins-|Statistics::Sequences::Joins> and L<Pot-|Statistics::Sequences::Pot>tests of sequences. 

Most methods are named with aliases, should you be used to referring to Perl statistics methods by one or another of the conventions. Present conventions are mostly based on those used in Juan Yun-Fang's modules, e.g., L<Statistics::ChisqIndep|Statistics::ChisqIndep>.

=head3 new

 $seq = Statistics::Sequences->new();

Returns a new Statistics::Sequences object by which all the methods for caching, dichotomising, and testing data can be accessed, including each of the methods for performing the L<Runs-|Statistics::Sequences::Runs>, L<Joins-|Statistics::Sequences::Joins> and L<Pot-|Statistics::Sequences::Pot>tests. The parameters C<corr>, C<tails> and C<precision_p> can be usefully set here, during construction, to be used by all tests.

Any one of the sub-packages, such as L<Statistics::Sequences::Runs|Statistics::Sequences::Runs>, can be individually imported, and its own L<new|new> method can be called, e.g.:

 use Statistics::Sequences::Runs;
 $runs = Statistics::Sequences::Runs->new();

In this case, data are not automatically shared across packages, and only one test (in this case, the Runs-test) can be accessed through the class-object returned by L<new|new>.

=cut

#-----------------------------------------------
sub new {
#-----------------------------------------------
    my ($proto, %args) = (shift, @_);
    my $class = ref($proto) || $proto;
    my $self= {};
    bless($self, $class);
    return $self;
}

=head2 Caching data

=head3 load

 $seq->load(@data); # Anonymous load
 $seq->load(\@data); # Anonymously referenced load
 $seq->load(blues => \@blue_scores, reds => \@red_scores); # Named loads
 $seq->load({blues => \@blue_scores, reds => \@red_scores}); # Same, but referenced

I<Aliases:> C<load_data>

Cache an anonymous list of data as an array-reference, or named data-sets as a hash reference, accessible as C<$seq-E<gt>{'data'}>, and available over any number of tests. Each call to L<load|load> removes whatever might have been previously loaded. Sending nothing deletes all loaded data (by C<undef>fing C<$seq-E<gt>{'data'}>); sending another list makes another set of data available for testing.

Anonymous and named loading, and function aliases, are provided given the variety of such methods throughout the Statistics modules. Telling the difference between an unreferenced array (for anonymous loading) and an unreferenced hash (for nominal loading) is simply performed on the basis of the second element: if it's a reference, the list is taken as a hash, otherwise as an array. Inelegant, but accommodating.

=cut

#-----------------------------------------------        
sub load {
#-----------------------------------------------        
    my $self = shift;

    $self->unload();

    my %tmp_dat = ();
    
    if (ref($_[0])) {
       $self->{'data'} = $_[0]; # datawhat => \@data } or \@data
    }
    elsif (ref($_[1])) { # (datawhat => \@data)
        $self->{'data'} = {@_};
    }
    else {
        $self->{'data'} = \@_; # (@data)
    }

    $self->{'_tested'} = 0;

    return $self;
}

=head3 add

 $seq->add_data($char1, $char2)
 $seq->add_data([$char1, $char2])
 $seq->add_data({ reds => 1})

I<Aliases:> C<add_data>, C<append>, C<append_data>

Just push any value(s) or so along, without clobbering what's already in there (as L<load_data|load_data> would).

=cut

#-----------------------------------------------        
sub add {
#-----------------------------------------------        
    my $self = shift;

    $self->{'testdata'} = undef;
    $self->{'data'} = undef;

    my %tmp_dat = ();
    
    if (ref($_[0])) {
       if (ref($_[0]) eq 'HASH') { 
             %tmp_dat = %{$_[0]};
       }
       else { 
            push @{$self->{'data'}}, $_[0];
       }
    }
    elsif (scalar @_) {
        if (ref($_[1])) {
            %tmp_dat = @_;
        }
        else {
            push @{$self->{'data'}}, @_;
        }
    }
        
    while (my($k, $v) = each %tmp_dat ) {
        push @{ $self->{'data'}->{$k}}, $v;
    }    

    $self->{'_tested'} = 0;
 
    return $self;
}

=head3 read

 $seq->read()

Alias: C<get_data>

Return the hash of data (just return $seq->{'data'}).

=cut

#-----------------------------------------------
sub read {
#-----------------------------------------------
    my $self = shift;
    return $self->{'data'};
}

=head3 unload

 $seq->unload()

Alias: C<clear_data>, C<delete_data>

Empty, clear, clobber what's in there. This is always performed ahead of any C<load>.

=cut

sub unload {
    my $self = shift;
    $self->{$_} = undef foreach qw/testdata data test observed_stdev df/;
    $self->{$_} = 0 foreach qw/samplings expected observed variance std_dev obs_dev z_value _tested/;
    $self->{$_} = 1 foreach qw/p_value _ccorr/;
    $self->{$_} = 2 foreach qw/_tails/;
}

=head2 Dichotomising data

Both the runs- and joins-tests expect dichotomous data, i.e., as if there were only two categorical variables. Numerical and multi-valued categorical data, once loaded, can be "reduced" to this format by the following methods, namely, L<cut|cut>, L<swing|swing>, L<pool|pool> and L<match|match>. Or supply data in this format. Both the runs- and joins-test will C<croak> if more (or less) than two states are found in the data.

Each method stores the data in the class object as an array-reference named "testdata", accessed so:

 print 'dichotomous data: ', @{$seq->{'testdata'}}, "\n";

=head3 Numerical data: Single-sample dichotomisation

=head4 cut

 $seq->cut(value => 'median', equal => 'gt'); # cut anonymously cached data at a central tendency
 $seq->cut(value => 23); # cut anonymously cached data at a specific value
 $seq->cut(value => 'mean', data => 'blues'); # cut named data at its average

I<This method is only suitable for numerical data.>

Reduce loaded data to two categories by cutting it about a certain value. For example, the following raw data, when cut for values greater than or equal to 5, yield the subsequent dichotomous series.

 @raw_data = (4, 3, 3, 5, 3, 4, 5, 6, 3, 5, 3, 3, 6, 4, 4, 7, 6, 4, 7, 3);
 @cut_data = (0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0);

The following options may be specified.

=over 8

=item value => 'mean|median|mode|\d+'

Specify the value at which the data will be cut. This could be the mean, median or mode (as calculated by L<Statistics::Lite|Statistics::Lite>), or a numerical value within the range of the data. The default is the I<median>. The cut-value, as specified by C<value>, can be retrieved thus:

 print $seq->{'cut_value'};

=item equal => 'I<gt>|I<lt>|I<0>'

This option specifies how to cut the data should the cut-value (as specified by C<value>) be present in the data. The default value is 0: observations equal to the cut-value are skipped. If C<equal =E<gt> I<gt>>: all data-values I<greater than or equal to> the cut-value will form one group, and all data-values less than the cut-value will form another. To cut with all values I<less than or equal to> in one group, and higher values in another, set this parameter to I<lt>.

=item data => 'I<string>'

Specify which named cached data-set to cut.

=back

=cut

#-----------------------------------------------
sub cut {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $dat = $self->_rawdata_aref($args->{'data'});
    
    # Ensure data can be numerically cut before sending off to Statistics-Lite for doing so (doesn't work if simply using Lite::sum):
    foreach (@{$dat}) {
        croak __PACKAGE__, '::cut All data must be numerical for dichotomizing about a cut-value' if ! looks_like_number($_);
    }
 
    # Get a cut-value:
    if ((!$args->{'value'}) or ($args->{'value'} =~ /^[a-z]+$/i and $args->{'value'} !~ /^m(e(an|dian)|ode)/)) {
        $args->{'value'} = 'median';
    }

    if ($args->{'value'} =~ /^[a-z]+$/i) {
        no strict 'refs';
        $self->{'cut_value'} = &{$args->{'value'}}(@{$dat}); # using Statistics-Lite methods
    }
    else {
        $self->{'cut_value'} = $args->{'value'};
    }

    # Find the number of observations above and below the cut_value:
    my @seqs = ();
    foreach (@{$dat}) {
        if ($_ > $self->{'cut_value'}) {
            push @seqs, 1;
        }
        elsif ($_ < $self->{'cut_value'}) {
            push @seqs, 0;
        }
        elsif ($_ == $self->{'cut_value'}) {
             if (defined $args->{'equal'}) {
                if ($args->{'equal'} eq 'gt') {
                    push @seqs, 1;
                }
                elsif ($args->{'equal'} eq 'lt') {
                    push @seqs, 0;
                }
                else {
                    next;
                }
            }
            else {
                next;
            }
        }
    }

    $self->{'testdata'} = \@seqs;

    return $self;
}

=head4 swing

 $seq->swing();
 $seq->swing(data => 'reds'); # if more than one are loaded, or a single one was loaded with a name

This is another transformation that, like the L<cut|cut> method, can be used to produce a dichotomous sequence from a single set of numerical data. You essentially test the degree of consistency of the rises and falls in the data. Each element in the named data-set is subtracted from its successor, and the result is replaced with a 1 if the difference represents an increase, or 0 if it represents a decrease. For example, the following numerical series produces the subsequent dichotomous series.

 @values = (qw/3 4 7 6 5 1 2/);
 @dicho =  (qw/1 1 0 0 0 1/);

Dichotomously, the data can be seen as commencing with an ascending run of length 2, followed by a descending run of length 3, and finishing with a short increase. Note that the number of resulting observations is less than the original number.

Note that the critical region of the distribution lies (only) in the upper-tail; a one-tailed test of significance is appropriate.

=over 8

=item equal => 'I<gt>|I<lt>|I<rpt>|I<0>'

The default result when the difference between two successive values is zero is to skip the observation, and move onto the next succession (C<equal =E<gt> 0>). Alternatively, you may wish to repeat the result for the previous succession; skipping only a difference of zero should it occur as the first result (C<equal =E<gt> 'rpt'>). Or, a difference greater than or equal to zero is counted as an increase (C<equal =E<gt> 'gt'>), or a difference less than or equal to zero is counted as a decrease. For example, 

 @values =    (qw/3 3 7 6 5 2 2/);
 @dicho_def = (qw/1 0 0 0/); # First and final results (of 3 - 3, and 2 - 2) are skipped
 @dicho_rpt = (qw/1 0 0 0 0/); # First result (of 3 - 3) is skipped, and final result repeats the former
 @dicho_gt =  (qw/1 1 0 0 0 1/); # Greater than or equal to zero is an increase
 @dicho_lt =  (qw/0 1 0 0 0 0/); # Less than or equal to zero is a decrease

=back

=cut

#-----------------------------------------------
sub swing {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $dat = $self->_rawdata_aref($args->{'data'});

    # Replace observations with the succession of rises and falls:
    my ($i, @seqs) = ();
    for ($i = 0; $i < (scalar @{$dat} - 1); $i++) {
        croak __PACKAGE__, '::swing All data must be numerical for dichotomizing into ups and downs' if ! looks_like_number($dat->[$i]);
        my $res = $dat->[($i + 1)] - $dat->[$i];
        if ($res > 0) {
            push @seqs, 1;
        }
        elsif ($res < 0) {
            push @seqs, 0;
        }
        else {
            if (defined $args->{'equal'}) {
                if ($args->{'equal'} eq 'rpt') {
                    push @seqs, $seqs[-1] if scalar @seqs; 
                }
                elsif ($args->{'equal'} eq 'gt') {
                    push @seqs, 1;
                }
                elsif ($args->{'equal'} eq 'lt') {
                    push @seqs, 0;
                }
                else {
                    next;
                }
            }
            else {
                next;
            }
        }
    }
    $self->{'testdata'} = \@seqs;
    #$self->{'tails'} = 1;
    return $self;
}

=head3 Numerical data: Two-sample dichotomisation

=head4 pool

 $seq->pool('data' => ['blues', 'reds']);

Constructs a single series out of two series of cached I<numerical> data as a ranked pool, i.e., by pooling the data from each series according to the magnitude of their values at each trial. This is the typical option when using the Wald-Walfowitz test for determining a difference between two samples. Specifically, the values from both samples are pooled and ordered from lowest to highest, and then clustered into runs according to the sample from which neighbouring values come from. Another run occurs wherever there is a change in the source of the values. A non-random effect of, say, higher or lower values consistently coming from one sample rather than another, would be reflected in fewer runs than expected by chance. See the C<ex/checks.pl> file in the installation distribution for a couple examples.

See also the methods for categorical data where it is ok to ignore any order and intervals in your numerical data.

=cut

#-----------------------------------------------
sub pool {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    croak __PACKAGE__, '::pool Two samples of data for pooling need to be loaded' if !$args->{'data'} or ref $args->{'data'} ne 'ARRAY' or !scalar @{$args->{'data'}} or scalar @{$args->{'data'}} != 2;

    my ($tot, %dat) = (0);

    foreach (@{$args->{'data'}}) {
        croak __PACKAGE__, '::pool Data named \'', $_, '\' for pooling do not exist' if ! @{$self->{'data'}->{$_}};
        $dat{$_} = [sort {$a <=> $b} @{$self->{'data'}->{$_}}];
        $tot += scalar(@{$dat{$_}});
    }
 
    my ($i, $aa, $bb, @testdata) = ();
    while (scalar(@testdata) < $tot) {
        $aa = $dat{$args->{'data'}->[0]}->[0];
        $bb = $dat{$args->{'data'}->[1]}->[0];
        #croak __PACKAGE__, '::pool All data must be numerical for dichotomizing into a ranked pool' if ! looks_like_number($aa) || ! looks_like_number($bb);
        $i = defined $aa && defined $bb ? $aa < $bb ? 0 : 1 : defined $aa ? 0 : 1;
        shift @{$dat{$args->{'data'}->[$i]}};
        push @testdata, $args->{'data'}->[$i];
    }
 
    $self->{'testdata'} = \@testdata;

    return $self;
}

=head3 Categorical data

=head4 match

 $seq->match('data' => ['blues', 'reds']);

Reduce two lists of loaded data to two categories in a single array, according to the match between the elements at each index. Where the data-values are equal at a certain index, they will be represented with a 1; otherwise a 0. Numerical or stringy data can be equated. For example, the following two arrays would be reduced to the third, where a 1 indicates a match of identical values in the two data sources.

 @blues = (qw/1 3 3 2 1 5 1 2 4/);
 @reds =  (qw/4 3 1 2 1 4 2 2 4/);
 @dicho = (qw/0 1 0 1 1 0 0 1 1/);

The following options may be specified.

=over 8

=item data => [qw/blues reds/]

Specify, a referenced array, two named data-sets, as previously passed to L<load|load>. An attempt to match a number of data-sets other than 2 will emit a C<croak>.

=item lag => I<integer> OR [I<integer>, B<I<loop>> (I<boolean>)] (where I<integer> < number of observations I<or> I<integer> > -1 (number of observations) ) 

Match the two data-sets by shifting the first named set ahead or behind the other data-set by C<lag> observations. The default is zero. For example, one data-set might be targets, and another responses to the targets:

 targets   =	cbbbdacdbd
 responses =	daadbadcce

Matched as a single sequence of hits (1) and misses (0) where C<lag> = B<0> yields (for the match on "a" in the 6th index of both arrays):

 0000010000

If C<lag> is set to B<+1>, however, each response is associated with the target one ahead of the trial for which it was observed; i.e., each target is shifted to its +1 index. So the first element in the above responses (I<d>) would be associated with the second element of the targets (I<b>), and so on. Now, matching the two data-sets with a B<+1> lag gives two hits, of the 4th and 7th elements of the responses to the 5th and 8th elements of the targets, respectively:

 000100100

Note that with a lag of zero, there were 3 runs of (a) hit/s or miss/es, but with a lag of 1, there were 5 runs.

Lag values can be negative, so that a B<-2> lag, for instance, will give a hit/miss series of:

 00101010

Here, responses necessarily start at the third element (I<a>), the first hits occurring when the fifth response-element corresponds to the the third target element (I<b>).

In the above example, the last response (I<e>) could not be used, and the number of elements in the hit/miss sequence became n-C<lag> less the original target sequence. This means that the maximum value of lag must be one less the size of the data-sets, or there will be no data.

You can, alternatively, preserve all lagged data by looping any excess to the start or end of the criterion data. The number of observations will then always be the same, regardless of the lag. Matching the data in the example above with a lag of +1, with looping, creates an additional match between the final response and the first target (I<d>):

 1000100100

To effect looping, send a referenced list of the lag and a boolean for the loop, e.g., :

 lag => [-3, 1]

=back

=cut

#-----------------------------------------------
sub match {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    my $dat = $self->_rawdata_aref($args->{'data'});

    croak __PACKAGE__, '::match Names of two loaded data-sets are required for matching' if scalar @{$dat} != 2 or !ref $dat->[0] || !ref $dat->[1];
    
    $dat = $self->lag($args->{'lag'}, $dat->[0], $dat->[1]) if $args->{'lag'};

    # Ensure the criterion data-set is the set with the least observations:
    my $ari = scalar @{$dat->[0]} <= scalar @{$dat->[1]} ? $dat->[0] : $dat->[1];

    my ($i, @seqs) = ();
    for ($i = 0; $i < scalar @{$ari}; $i++) {
        next if !defined $dat->[0]->[$i] || !defined $dat->[1]->[$i];
        $seqs[$i] = $dat->[0]->[$i] eq $dat->[1]->[$i] ? 1 : 0;
    }
    $self->{'testdata'} = \@seqs;
    #$self->{'tails'} = 2;
    return $self;
}

#-----------------------------------------------
sub lag {
#-----------------------------------------------
    my ($self, $lag, $t, $r) = @_;
   
    my $loop;
    if (ref $lag) {
        $loop = $lag->[1];
        $lag = $lag->[0];
    }
    else {
        $loop = 0;
    }
   
    return [$t, $r] if !$lag or abs($lag) >= scalar @{$t};
    
    my @tgt = @{$t};
    my @rsp = @{$r};
    
    if ($lag > 0) {
        foreach (1 .. abs($lag) ) {
            if ($loop) {
                unshift(@rsp, pop @rsp);
            }
            else {
                shift @tgt;
                pop @rsp;
            }
        }
    }
    elsif ($lag < 0) {
        foreach (1 .. abs($lag) ) {
            if ($loop) {
                push(@rsp, shift @rsp);
            }
            else {
                pop @tgt;
                shift @rsp;
            }
        }
    }
    return [\@tgt, \@rsp];
}

=head2 Testing data

=head3 test

 $seq->test(what => 'runs');
 $seq->test(what => 'joins');
 $seq->test(what => 'pot', state => 'a value appearing in the testdata');
 $seq->test(what => 'vnomes', length => 'an integer greater than zero and less than sample-size');

 $runs->test();
 $joins->test(prob => 1/2);
 $pot->test(state => 'circle');
 $vnomes->test(length => 3);

I<Alias:> C<process>

When using a Statistics::Sequences class-object, this method requires naming which test to perform, i.e., runs, joins, pot or vnomes. This is I<not> required when the class-object already refers to one of the sub-modules, as created by the C<new> method within L<Statistics::Sequences::Runs|Statistics::Sequences::Runs/new>, L<Statistics::Sequences::Joins|Statistics::Sequences::Joins/new>, L<Statistics::Sequences::Pot|Statistics::Sequences::Pot/new> and L<Statistics::Sequences::Vnomes|Statistics::Sequences::Vnomes/new>.

Note that simply giving the first argument as the name, unkeyed, is deprecated, and will be removed in the next version; i.e., don't just use, e.g.:

 $seq->test('runs'); # verboten - nicht mehr erlaubt!

=head4 General options

Options to L<test|test> available to all the sub-package tests are as follows.

=over 8

=item data => 'I<string>'

Optionally specify the name of the data to be tested. By default, this is not required: the testdata are those that were last loaded, either anonymously, or as given to one of the dichotomising methods. Otherwise, I<if the data are already ready for testing in a dichotomous format>, data that were previously loaded by name can be individually tested. For example, here are two sets of data that are loaded by name, and then a single test of one of them is performed.

 @chimps = (qw/banana banana cheese banana cheese banana banana banana/);
 @mice = (qw/banana cheese cheese cheese cheese cheese cheese cheese/);
 $seq->load(chimps => \@chimps, mice => \@mice);
 $seq->test(what => 'runs', data => 'chimps')->dump();

=item ccorr => I<boolean>

Specify whether or not to perform the continuity-correction on the observed deviation. Default is false. Relevant only for those tests relying on a I<Z>-test. See L<Statistics::Zed>.

=item tails => I<1>|I<2>

Specify whether the I<z>-value is calculated for both sides of the normal (or chi-square) distribution (2, the default for most testdata) or only one side (the default for data prepared with the B<swing> method.

=back

=head4 Test-specific required settings and options

Some sub-package tests need to have parameters defined in the call to L<test|test>, and/or have specific options, as follows.

B<Joins> : The Joins-test I<optionally> allows the setting of a probability value; see L<test|test> in the  L<Statistics::Sequences::Joins|Statistics::Sequences::Joins/test> manpage.

B<Pot> : The Pot-test I<requires> the setting of a state to be tested; see C<test> in the  L<Statistics::Sequences::Pot|Statistics::Sequences::Pot/test> manpage.

B<Vnomes> : The Seriality test for V-nomes I<requires> a length, i.e., the value of I<v>; see C<test> in the L<Statistics::Sequences::Vnomes|Statistics::Sequences::Vnomes/test> manpage..

B<Runs> : There are presently no specific requirements nor options for the Runs-test. 

=cut

#-----------------------------------------------
sub test {
#-----------------------------------------------
    my ($self, $testname, $args) = (shift);

    if (! ref $_[0] and scalar(@_) % 2) { # uneven: assume first value is the testname:
        $testname = shift;
        $args = {@_};
        cluck __PACKAGE__, "::test Giving an unkeyed test-name as the first argument is deprecated; the name should be hash-keyed as what => '$testname'";
    }
    else {
        $args = ref $_[0] ? $_[0] : {@_};
        $testname = $args->{'what'};
    }

    my $class = __PACKAGE__ . '::' . ucfirst($testname);
    eval "require $class";
    if (!$@) {
        $self->{'test'} = $testname;
        bless($self, $class);
        $self->test($args);
        bless($self, __PACKAGE__);
    }
    else {
        croak __PACKAGE__, "::test Requested test $class is not valid (should be joins, pot, runs or vnomes), or the package is not locatable"
    }
    return $self;
}

=head2 Accessing results

All relevant statistical values are "lumped" into the class-object, and can be retrieved thus:

 $seq->{'observed'} # The observed value of the test-statistic (Runs, Joins, Pot)
 $seq->{'expected'} # The expected value of the test-statistic (Runs, Joins, Pot)
 $seq->{'obs_dev'} # The observed deviation (observed minus expected values), continuity-corrected, if so specified
 $seq->{'std_dev'} # The standard deviation
 $seq->{'variance'} # Variance
 $seq->{'z_value'} # The value of the z-statistic (ratio of observed to standard deviation), where relevant
 $seq->{'p_value'} # The "normal probability" associated with the z-statistic, or other statistic

=head3 dump

 $seq->dump(flag => '1|0', text => '0|1|2', precision_s => 'integer', precision_p => 'integer');

I<Alias:> C<print_summary>

Print results of the last-conducted test to STDOUT. By default, if no parameters to C<dump> are passed, a single line of test statistics are printed. Options are as follows.

=over 8

=item flag => I<boolean>

If true, the I<p>-value associated with the I<z>-value is appended with a single asterisk if the value if below .05, and with two asterisks if it is below .01.

If false (default), nothing is appended to the I<p>-value.

=item text => I<0>|I<1>|I<2>

If set to 1, a single line is printed, beginning with the name of the test, then the observed and expected values of the test-statistic, and the I<z>-value and its associated I<p>-value. The Pot-test, additionally, shows the state tested in parentheses after the test-name, and the test for I<v>-nomes shows the tested length in parentheses after the test-name. For example:

 Joins: expected = 400.00, observed = 360.00, z = -2.83, p = 0.0023389**
 Runs: expected = 398.86, observed = 361.00, z = -2.70, p = 0.0070374**
 Pot(1): expected = 288.51, observed = 303.63, z = 2.64, p = 0.0082920**

If set to anything greater than 1, more verbose info is printed: each of the above bits of info are printed, on separate lines, with some additional info about sample-size, etc.

If set to zero, only the string returned by C<string> is printed, with no terminal line-break.

=item precision_s => 'I<non-negative integer>'

Precision of the z-statistic.

=item precision_p => 'I<non-negative integer>'

Specify rounding of the probability associated with the I<z>-value to so many digits. If zero or undefined, you get everything available.

=back

=cut

#-----------------------------------------------
sub dump {
#-----------------------------------------------
    my ($self) = (shift);
    my $class = __PACKAGE__ . '::' . ucfirst($self->{'test'});

    if (!$@) {
        bless($self, $class);
        $self->dump(@_);
        bless($self, __PACKAGE__);
    }
    else {
        croak __PACKAGE__, '::test Name for valid test (joins, runs, pot) required'
    }
}

=head3 dump_data

 $seq->dump_data(delim => "\n")

Prints to STDOUT a space-separated line of the testdata - as dichotomised and put to test. Optionally, give a value for C<delim> to specify how the datapoints should be separated.

=cut

#-----------------------------------------------
sub dump_data {
#-----------------------------------------------
    my ($self, %args) = @_; 
    my $delim = $args{'delim'} || " ";
    print join($delim, @{$self->{'testdata'}}), "\n";
}

=head3 string

 $seq->string()

Returns a single line giving the I<z>-value and I<p>-value. Accepts the C<precision_s>, C<precision_p> and C<flag> options, as for L<dump|dump>.

=cut

#-----------------------------------------------
sub string {
#-----------------------------------------------
    my ($self) = (shift);
    my $args = ref $_[0] ? $_[0] : [@_];
    my $stat_name = delete $args->{'stat_name'} || 'z_value';
    return if !defined $self->{$stat_name} || !defined $self->{'p_value'};
    $args->{'tails'} ||= $self->{'tails'};
    my $str = $stat_name eq 'z_value' ? 'Z' : $stat_name;
    $str .= ' (' . $self->{'df'} . ')' if defined $self->{'df'};
    $str .= ' = ';
    $str .= _precisioned($args->{'precision_s'}, $self->{$stat_name});
    $str .=  ', ' . ($args->{'tails'}||$self->{'_tails'}||2) . 'p = ';
    $str .= _precisioned($args->{'precision_p'}, $self->{'p_value'});
    $str .= ($self->{'p_value'} < .05 ? ($self->{'p_value'} < .01 ? '**' : '*') : '') if $args->{'flag'};
    return $str;
}

# PRIVATMETHODEN

sub _dump_pass {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : [@_];
    if ($args->{'text'} and $args->{'text'} > 1) {
        $self->_dump_verbose($args);
    }
    else {
        $self->_dump_sparse($args);
    }
}   

sub _dump_sparse {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : [@_];
    if ($args->{'text'}) { # equal 1:
        print $args->{'testname'} . ': ' if defined $args->{'testname'} and length $args->{'testname'};
        printf("observed = %.3f", $self->{'observed'});
        printf $self->{'observed_stdev'} ? (" (%.3f), ", $self->{'observed_stdev'}) : ', ';
        printf("expected = %.3f, ", $self->{'expected'});
        print $self->string($args);
        print "\n";
    }
    else {
        print $self->string($args);
    }
}

sub _dump_verbose {
    my ($self, $args) = @_;
    printf(" Observed %s = %.3f" , $args->{'testname'}, $self->{'observed'});
    printf $self->{'observed_stdev'} ? (" (%.3f)\n", $self->{'observed_stdev'}) : "\n";
    printf(" Expected %s = %.3f\n" , $args->{'testname'}, $self->{'expected'});
    print ' ', $self->string($args);
    print "\n";
}

sub _rawdata_aref { 
    my ($self, $data_name) = @_;
    my $dat;
    
    if ($data_name) {# Has a data array been named by the user?
       if (ref $data_name) {# Yes, but there might be more than one (e.g, for matching or pooling):
           my $i = 0;
           foreach (@{$data_name}) {
              if (ref $self->{'data'} eq 'HASH' and defined $self->{'data'}->{$_}) {
                 $dat->[$i++] = $self->{'data'}->{$_};
              }
              else {
                 croak __PACKAGE__, '  Data named \'', $_, '\' are not loaded';
              }
           }
       }
       elsif (length $data_name) {# One name only, prob. called for cutting:
           if (defined $self->{'data'}->{$data_name}) {
                $dat = $self->{'data'}->{$data_name};
           }
           else {
                croak __PACKAGE__, '  Data named \'', $data_name, '\' are not loaded';
           }
       }
    }
    elsif (ref $self->{'data'} eq 'ARRAY' and scalar(@{$self->{'data'}}) ) {# No named array; get any data stored as "data":
           $dat = $self->{'data'};
    }
    else {
    # if !$args->{'data'} or ref $args->{'data'} ne 'ARRAY' or ! scalar @{$args->{'data'}};
       croak __PACKAGE__, ' No data are accessible for testing - Do they need to be loaded?';
    }

    return $dat;
}

sub _testdata_aref {
    my ($self, $args) = @_;

    if ($args->{'data'} and ref($self->{'data'}->{$args->{'data'}}) && scalar @{$self->{'data'}->{$args->{'data'}}}) {
        $self->{'testdata'} = $self->{'data'}->{$args->{'data'}};
    }
    elsif(!ref $self->{'testdata'} or ref $self->{'testdata'} ne 'ARRAY') {
        if (ref($self->{'data'}) eq 'ARRAY') {
            $self->{'testdata'} = $self->{'data'};
        }
        else {
            croak __PACKAGE__, '::test Data for testing have not been specified or loaded';
        }
    }
}

sub _expound { # Get the expectation, variance & observed N-sequences from each test; test, and lump results into class object
    my ($self, $obs_val, $exp_val, $var, $args) = @_;
    $args->{'tails'} ||= 2;
    $args->{'ccorr'} = 1 if ! defined $args->{'ccorr'};
    my $dev = Statistics::Zed->new(
        ccorr => $args->{'ccorr'}, 
        tails => $args->{'tails'}, 
        precision_s => $args->{'precision_s'}, 
        precision_p => $args->{'precision_p'},
    );

    my ($z, $pz, $obs_dev, $std_dev) = $dev->zscore(observed => $obs_val, expected => $exp_val, variance => $var) if $var;

    # Lump values into class object:
    $self->{'observed'} = $obs_val;
    $self->{'expected'} = $exp_val;
    $self->{'z_value'} = $z;
    $self->{'p_value'} = $pz;
    $self->{'obs_dev'} = $obs_dev;
    $self->{'std_dev'} = $std_dev;
    $self->{'variance'} = $var;

    $self->{'_tested'} = 1;
    $self->{'_tails'} = $args->{'tails'};
    $self->{'_ccorr'} = $args->{'ccorr'};
}

sub _expire {# Do the least amount possible, if variance = 0 or some other nuisance problem in the tests:
    my ($self, $obs_val, $exp_val, $args);
    $self->{$_} = 0 foreach qw/observed z_value obs_dev std_dev variance/;
    $args->{'precision_p'} ||= 0;
    $self->{'p_value'} = sprintf('%.' . $args->{'precision_p'} . 'f', 1);
    $self->{'observed'} = $obs_val;
    $self->{'expected'} = $exp_val;
}

sub _precisioned {
    return $_[0] ? sprintf('%.' . $_[0] . 'f', $_[1]) : (defined $_[1] ? $_[1] : ''); # don't lose any zero
}

# Aliases:
*load_data = \&load;
*add_data = \&add;
*append = \&add;
*append_data = \&add;
*get_data = \&read;
*clear_data = \*unload;
*delete_data = \*unload;
*process = \&test;
*print_summary = \&dump;

1;
__END__

=head1 REFERENCES

Burdick, D. S., & Kelly, E. F. (1977). Statistical methods in parapsychological research. In B. B. Wolman (Ed.), I<Handbook of Parapsychology> (pp. 81-130). New York, NY, US: Van Nostrand Reinhold. [Description of joins-test, with comparision to runs-test.]

Kelly, E. F. (1982). On grouping of hits in some exceptional psi performers. I<Journal of the American Society for Psychical Research>, I<I76>, 101-142. [Application of runs-test, with discussion of normality issue.]

Schmidt, H. (2000). A proposed measure for psi-induced bunching of randomly spaced events. I<Journal of Parapsychology, 64,> 301-316. [Describes the pot-test.]

Swed, F., & Eisenhart, C. (1943). Tables for testing randomness of grouping in a sequence of alternatives. I<Annals of Mathematical Statistics>, I<14>, 66-87. [Look in C<ex/checks.pl> in the installation dist for a few examples from this paper for testing.]

Wald, A., & Wolfowitz, J. (1940). On a test whether two samples are from the same population. I<Annals of Mathematical Statistics>, I<11>, 147-162. [Describes the runs-test.]

Wishart, J. & Hirshfeld, H. O. (1936). A theorem concerning the distribution of joins between line segments. I<Journal of the London Mathematical Society>, I<11>, 227. [Describes the joins-test.]

Wolfowitz, J. (1943). On the theory of runs with some applications to quality control. I<Annals of Mathematical Statistics>, I<14>, 280-288. [Suggests some ways in which data may be dichotomised for testing runs; implemented here.]

=head1 SEE ALSO

L<Statistics::Burst|Statistics::Burst> : Another test of sequences.

=head1 TO DO/BUGS

Results are dubious if there are only two observations.

Testing not by I<z>-scores, and/or using poisson distribution for low number of observations

Fu's Markovian solution

Multivariate extensions

Sort option for pool method ?

=head1 REVISIONS

The series testing methods (series_init, series_update and series_test) have been moved to Statistics::Zed as of v0.03.

Simply giving the first argument to C<test> as the name of the test, unkeyed, is deprecated, and will be removed in the next version.

See CHANGES file in installation dist.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2010 Roderick Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut

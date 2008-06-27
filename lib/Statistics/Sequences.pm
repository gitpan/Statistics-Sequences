package Statistics::Sequences;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak);
use Class::OOorNO qw(coerce_array);
use Statistics::Descriptive;
use Statistics::Zed 0.01;
use Statistics::Lite qw(:funcs);

our $VERSION = '0.02';

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub new {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my ($proto, %args) = (shift, @_);
    my $class = ref($proto) || $proto;
    my $self= {};

    $self->{'data'} = '';
    $self->{'test'} = '';
    
    foreach (qw/samples expected observed variance std_dev obs_dev/) {
        $self->{$_} = $args{$_} || 0;
    }
    foreach (qw/p_value ccorr/) {
        $self->{$_} = $args{$_} || 1;
    }
    foreach (qw/tails/) {
        $self->{$_} = $args{$_} || 2;
    }
    bless($self, $class);
    return $self;
}

sub clear {
    my $self = shift;
    $self->{'data'} = '';
    $self->{'test'} = '';
    foreach (qw/samples expected observed variance std_dev obs_dev z_value/) {
        $self->{$_} = 0;
    }
    foreach (qw/p_value/) {
        $self->{$_} = 1;
    }
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
sub load {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
    my $self = shift;

    $self->{'testdata'} = undef;
    $self->{'data'} = undef;
    ##$self->{'tails'} = undef;

    my %tmp_dat = ();
    
    if (ref($_[0])) {
       if (ref($_[0]) eq 'HASH') { 
             %tmp_dat = %{$_[0]};
       }
       else { 
            $self->{'data'} = $_[0];
       }
    }
    elsif (scalar @_) {
        if (ref($_[1])) {
            %tmp_dat = @_;
        }
        else {
            $self->{'data'} = \@_;
        }
    }
        
    while (my($k, $v) = each %tmp_dat ) {
        $self->{'data'}->{$k} = $v;
    }    

    return $self;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub get_data {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    return $self->{'data'};
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub cut {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = coerce_array(@_);
    my $dat = $self->_get_data($args->{'data'});

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
    ##$self->{'tails'} = 2;
    return $self;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub match {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = coerce_array(@_);
    my $dat = $self->_get_data($args->{'data'});

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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub pool {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = coerce_array(@_);
    croak __PACKAGE__, '::pool Two samples of data for pooling need to be loaded'
    if !$args->{'data'} or ref $args->{'data'} ne 'ARRAY' or !scalar @{$args->{'data'}} or scalar @{$args->{'data'}} != 2;

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
        $i = defined $aa && defined $bb ? $aa < $bb ? 0 : 1 : defined $aa ? 0 : 1;
        shift @{$dat{$args->{'data'}->[$i]}};
        push @testdata, $args->{'data'}->[$i];
    }
 
    $self->{'testdata'} = \@testdata;

    return $self;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub swing {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = coerce_array(@_);
    my $dat = $self->_get_data($args->{'data'});

    # Replace observations with the succession of rises and falls:
    my ($i, @seqs) = ();
    for ($i = 0; $i < (scalar @{$dat} - 1); $i++) {
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub lag {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub test {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my ($self, $testname) = (shift, shift);
    my $class = __PACKAGE__ . '::' . ucfirst($testname);
    eval "require $class";
    if (!$@) {
        $self->{'test'} = $testname;
        bless($self, $class);
        $self->test(@_);
        bless($self, __PACKAGE__);
    }
    else {
        croak __PACKAGE__, '::test Name for valid test (joins, runs, pot) required'
    }
    return $self;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub series_init {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my ($self) = (shift);
    $self->{'series_stat'}->{$_} = Statistics::Descriptive::Sparse->new() foreach qw/observed expected variance/;
    $self->{'series'}->{$_} = 0 foreach qw/observed z_value r_value obs_dev std_dev variance/;
    $self->{'series'}->{'p_value'} = 1;##$self->{'p_value'} = sprintf('%.' . $self->{'p_precision'} . 'f', 1);
    $self->{'series'}->{'expected'} = 1;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub series_update {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my ($self) = (shift);
    $self->{'series_stat'}->{$_}->add_data($self->{$_}) foreach qw/observed expected variance/;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub series_test {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my ($self) = (shift);
    my $dev = Statistics::Zed->new();
    my $o = $self->{'series_stat'}->{'observed'}->sum();
    my $e = $self->{'series_stat'}->{'expected'}->sum();
    my $v = $self->{'series_stat'}->{'variance'}->sum();
    my ($z, $pz, $obs_dev, $std_dev) = $dev->zscore(
        observed => $o,
        expected => $e,
        variance => $v);
    $self->{'series'}->{'observed'} = $o;
    $self->{'series'}->{'expected'} = $e;
    $self->{'series'}->{'z_value'} = $z;
    $self->{'series'}->{'p_value'} = $pz;
    $self->{'series'}->{'r_value'} = $dev->z_2_r($z);
    $self->{'series'}->{'obs_dev'} = $obs_dev;
    $self->{'series'}->{'std_dev'} = $std_dev;
    $self->{'series'}->{'variance'} = $v;
    return $self; 
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub string {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : Class::OOorNO::coerce_array(@_);

    my $str = 'z = ';
    if ($args->{'s_precision'}) {
        $str .= sprintf('%.' . $args->{'s_precision'} . 'f', $self->{'z_value'});
    }
    else {
        $str .= $self->{'z_value'};
    }
    $str .=  ', ' . $self->{'tails'} . '-p = ';
    if ($args->{'p_precision'}) {
        $str .= sprintf('%.' . $args->{'p_precision'} . 'f', $self->{'p_value'});
    }
    else {
        $str .= $self->{'p_value'};
    }
    $str .= ($self->{'p_value'} < .05 ? ($self->{'p_value'} < .01 ? '**' : '*') : '') if $args->{'flag'};
    return $str;
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub dump_data {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my ($self, %args) = @_; 
    my $delim = $args{'delim'} || " ";
    print join($delim, @{$self->{'testdata'}}), "\n";
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub dump {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

#-----------------------------------------------
sub _dump_pass {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : Class::OOorNO::coerce_array(@_);
    $args->{'text'} = 1 if !defined $args->{'text'} && !$args->{'data'};

    print ref $self->{'testdata'} ? @{$self->{'testdata'}} : ref $self->{'data'} ? @{$self->{'data'}} : '', "\n" if $args->{'data'};

    if ($args->{'text'}) {
        if ($args->{'text'} > 1) {
            $self->_dump_verbose($args);
        }
        else {
            $self->_dump_sparse($args);
        }
    }
}   

#-----------------------------------------------
sub _dump_sparse {
#-----------------------------------------------
    my $self = shift;
    my $args = shift;##ref $_[0] ? $_[0] : Class::OOorNO::coerce_array(@_);
    print $args->{'testname'} . ': ' if defined $args->{'testname'} and length $args->{'testname'};
    print  'expected = ' . sprintf('%.2f', $self->{'expected'}) . ', observed = ' . sprintf('%.2f', $self->{'observed'}) . ', ';
    print $self->string($args);
    print "\n";
}

#-----------------------------------------------
sub _dump_verbose {
#-----------------------------------------------
    my ($self, $args) = @_;
    print ' Observed ' . $args->{'testname'} . ' = ' . sprintf('%.2f', $self->{'observed'}) . "\n" .
          ' Expected ' . $args->{'testname'} . ' = ' . sprintf('%.2f', $self->{'expected'}) . "\n" .
          ' Observed deviation = ' . sprintf('%.2f', $self->{'obs_dev'}) . "\n" .
          ' Standard deviation = ' . sprintf('%.2f', $self->{'std_dev'}) . "\n";
    print $self->string($args);
    print "\n";
}

#-----------------------------------------------
sub _get_data { 
#-----------------------------------------------
    my ($self, $data_name) = @_;
    my $dat;
   
    # Has a data array been named by the user?
    if ($data_name) {
       # Yes, but there might be more than one (e.g, for matching or pooling):
       if (ref $data_name) {
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
       # One name only, prob. called for cutting:
       elsif (length $data_name) {
           if (defined $self->{'data'}->{$data_name}) {
                $dat = $self->{'data'}->{$data_name};
           }
           else {
                croak __PACKAGE__, '  Data named \'', $data_name, '\' are not loaded';
           }
       }
    }
    # No named array; get any data stored as "data":
    elsif (ref $self->{'data'} eq 'ARRAY' and scalar(@{$self->{'data'}}) ) {
           $dat = $self->{'data'};
    }
    else {
    # if !$args->{'data'} or ref $args->{'data'} ne 'ARRAY' or ! scalar @{$args->{'data'}};
       croak __PACKAGE__, ' No data are accessible for testing - Do they need to be loaded?';
    }

    return $dat;
}

#-----------------------------------------------
sub _check_testdata {
#-----------------------------------------------
    my ($self, $args) = @_;
    # Try very hard to avoid doing anything:
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

# Aliases:
*load_data = \&load;
*add_data = \&load;
*process = \&test;
*print_summary = \&dump;

1;
__END__

=pod

=head1 NAME

Statistics::Sequences - Tests of sequential structure in the form of runs, joins, bunches, etc.

=head1 SYNOPSIS

  use Statistics::Sequences;
  $seq = Statistics::Sequences->new();
  $seq->load([1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1]); # dichotomous values
  $seq->test('runs')->dump(); # or 1st argument to test is 'joins' or 'pot'
  # (prints:)
  # Runs: expected = 7.00, observed = 7.00, z = -0.30, p = 0.76206

=head1 DESCRIPTION

Loading and preparing data for statistical tests of their sequential structure via L<Statistics::Sequences::Runs|Statistics::Sequences::Runs>, L<Statistics::Sequences::Joins|Statistics::Sequences::Joins>, and L<Statistics::Sequences::Pot|Statistics::Sequences::Pot>. Examples of the use of each test are given in these pages.

In general, to access the tests, you L<use|perlfunc/use> this base module to directly create a Statistics::Sequences object with the L<new|new> method, L<load|load> data into it, and then access each test by calling the L<test|test> method. This approach is useful for running several tests on the same data, as the data are immediately available to each test (of runs, pot and joins). See the L<SYNOPSIS|Statistics::Sequences/SYNOPSIS> for a simple example.

If you only want to perform a test of one type (e.g., runs), you might want to simply L<use> the relevant sub-package, create a class object specific to it, and load data specfically for its use; see the SYNOPSIS for the particular test, i.e., L<Runs|Statistics::Sequences::Runs/SYNOPSIS>, L<Joins|Statistics::Sequences::Joins/SYNOPSIS> or L<Pot|Statistics::Sequences::Pot/SYNOPSIS>. You won't be able to access other tests by this approach, unless you create another object for that test, and then specifically pass the data from the earlier object into the new one.

Note also that there are methods to anonymously or nominally cache data, and that data might need to be reduced to a dichotomous format, before a valid test can be run. Several dichotomising methods are provided, once data are loaded, and accessible via the generic or specific class objects, as above.

=head1 METHODS

=head2 Interface

The package provides an object-oriented interface for performing the L<Runs-|Statistics::Sequences::Runs>, L<Joins-|Statistics::Sequences::Joins> and L<Pot-|Statistics::Sequences::Pot>tests of sequences. 

Most methods come with aliases, should you be used to referring to Perl statistics methods by one or another of the conventions. Present methods are mostly based on those used in Juan Yun-Fang's modules, e.g., L<Statistics::ChisqIndep|Statistics::ChisqIndep>.

=head3 new

 $seq = Statistics::Sequences->new();

Returns a new Statistics::Sequences object by which all the methods for caching, dichotomising, and testing data can be accessed, including each of the methods for performing the L<Runs-|Statistics::Sequences::Runs>, L<Joins-|Statistics::Sequences::Joins> and L<Pot-|Statistics::Sequences::Pot>tests. The parameters C<corr>, C<tails> and C<p_precision> can be usefully set here, during construction, to be used by all tests.

Any one of the sub-packages, such as L<Statistics::Sequences::Runs|Statistics::Sequences::Runs>, can be individually imported, and its own L<new|new> method can be called, e.g.:

 use Statistics::Sequences::Runs;
 $runs = Statistics::Sequences::Runs->new();

In this case, data are not automatically shared across packages, and only one test (in this case, the Runs-test) can be accessed through the class-object returned by L<new|new>.

=head2 Caching data

=head3 load

 $seq->load(@data); # Anonymous load
 $seq->load(\@data); # Anonymously referenced load
 $seq->load(blues => \@blue_scores, reds => \@red_scores); # Named loads
 $seq->load({blues => \@blue_scores, reds => \@red_scores}); # Same, but referenced

I<Aliases:> C<load_data>, C<add_data>

Cache an anonymous list of data as an array-reference, or named data-sets as a hash reference, accessible as C<$seq-E<gt>{'data'}>, and available over any number of tests. Each call to L<load|load> removes whatever might have been previously loaded. Sending nothing deletes all loaded data (by C<undef>fing C<$seq-E<gt>{'data'}>); sending another list makes another set of data available for testing.

Anonymous and named loading, and function aliases, are provided given the variety of such methods throughout the Statistics modules. Telling the difference between an unreferenced array (for anonymous loading) and an unreferenced hash (for nominal loading) is simply performed on the basis of the second element: if it's a reference, the list is taken as a hash, otherwise as an array. Inelegant, but accommodating.

=head2 Dichotomising data

Both the runs- and joins-tests expect dichotomous data, i.e., as if there were only two categorical variables. Numerical and multi-valued categorical data, once loaded, can be "reduced" to this format by the following methods, namely, L<cut|cut>, L<match|match> and L<pool|pool>. Or supply data in this format. Both the runs- and joins-test will C<croak> if more (or less) than two states are found in the data.

Each method stores the data in the class object as an array-reference named "testdata", accessed so:

 print 'dichotomous data: ', @{$seq->{'testdata'}}, "\n";

=head3 cut

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

=head3 match

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

Match the two data-sets by shifting the first named set ahead or behind the other data-set by I<lag> observations. The default is zero. For example, one data-set might be targets, and another responses to the targets:

 targets   =	cbbbdacdbd
 responses =	daadbadcce

Matched as a single sequence of hits (1) and misses (0) where I<lag> = B<0> yields (for the match on "a" in the 6th index of both arrays):

 0000010000

If C<lag> is set to B<+1>, however, each response is associated with the target one ahead of the trial for which it was meant; i.e., each target is shifted to its +1 index. So the first element in the above responses (I<d>) would be associated with the second element of the targets (I<b>), and so on. Now, matching the two data-sets with a B<+1> lag gives two hits, of the 4th and 7th elements of the responses to the 5th and 8th elements of the targets, respectively:

 000100100

Note that with a lag of zero, there were 3 runs of (a) hit/s or miss/es, but with a lag of 1, there were 5 runs.

Lag values can be negative, so that a B<-2> lag, for instance, will give a hit/miss series of:

 00101010

Here, responses necessarily start at the third element (I<a>), the first hits occurring when the fifth response-element corresponds to the the third target element (I<b>).

In the above example, the last response (I<e>) could not be used, and the number of elements in the hit/miss sequence became n-I<lag> less the original target sequence. This means that the maximum value of lag must be one less the size of the data-sets, or there will be no data.

You can, alternatively, preserve all lagged data by looping any excess to the start or end of the criterion data. The number of observations will then always be the same, regardless of the lag. Matching the data in the example above with a lag of +1, with looping, creates an additional match between the final response and the first target (I<d>):

 1000100100

To effect looping, send a referenced list of the lag and a boolean for the loop, e.g., :

 lag => [-3, 1]

=back

=head3 pool

 $seq->pool('data' => ['blues', 'reds']);

Reduce two sets of cached I<numerical> data to two categories in a single array by pooling the data according to the magnitude of the values at each trial. This is the typical option when using the Wald-Walfowitz test for determining a difference between two samples. Specifically, the values from both samples are pooled and ordered from lowest to highest, and then clustered into runs according to the sample from which neighbouring values come from. Another run occurs wherever there is a change in the source of the values. A non-random effect of, say, higher or lower values consistently coming from one sample rather than another, would be reflected in fewer runs than expected by chance. See the C<ex/checks.pl> file in the installation distribution for a couple examples.

=head3 swing

 $seq->swing();
 $seq->swing(data => 'reds'); # if more than one are loaded, or a single one was loaded with a name

This is another transformation that, like the L<cut|cut> method, can be used to produce a dichotomous sequence from a single set of numerical data. You essentially test the degree of consistency of the rises and falls in the data. Each element in the named data-set is subtracted from its successor, and the result is replaced with a 1 if the difference represents an increase, or 0 if it represents a decrease. For example, the following numerical series produces the subsequent dichotomous series.

 @values = (qw/3 4 7 6 5 1 2/);
 @dicho =  (qw/1 1 0 0 0 1/);

Dichotomously, the data can be seen as commencing with an ascending run of length 2, followed by a descending run of length 3, and finishing with a short increase. Note that the number of resulting observations is less than the original number.

Note that the critical region of the distribution lies (only) in the upper-tail; a one-tailed test of significance is appropriate.

=over 8

=item equal => 'I<gt>|I<lt>|I<rpt>|I<0>'

The default result when the difference between two successive values is zero is to skip the observation, and move onto the next succession (C<equal =E<gt> 0>). Alternatively, you may wish to repeat the result for the previous succession; skipping only a difference of zero should it occur as the first result ((C<equal =E<gt> 'rpt'>)). Or, a difference greater than or equal to zero is counted as an increase (C<equal =E<gt> 'gt'>), or a difference less than or equal to zero is counted as a decrease. For example, 

 @values =    (qw/3 3 7 6 5 2 2/);
 @dicho_def = (qw/1 0 0 0/); # First and final results (of 3 - 3, and 2 - 2) are skipped
 @dicho_rpt = (qw/1 0 0 0 0/); # First result (of 3 - 3) is skipped, and final result repeats the former
 @dicho_gt =  (qw/1 1 0 0 0 1/); # Greater than or equal to zero is an increase
 @dicho_lt =  (qw/0 1 0 0 0 0/); # Less than or equal to zero is a decrease

=back

=head2 Testing data

=head3 test

 $seq->test('runs');
 $seq->test('joins');
 $seq->test('pot', state => 'a value appearing in the testdata');

 $runs->test();
 $joins->test(prob => 1/2);
 $pot->test(state => 'circle');

I<Alias:> C<process>

When using a Statistics::Sequences class-object, this method requires specification of which test to perform, i.e., runs, joins or pot; the name of the test is simply given as the first argument. This is I<not> required when the class-object already refers to one of the sub-modules, as created by the C<new> method within L<Statistics::Sequences::Runs|Statistics::Sequences::Runs/new>, L<Statistics::Sequences::Joins|Statistics::Sequences::Joins/new>, and L<Statistics::Sequences::Pot|Statistics::Sequences::Pot/new>.

=head4 General options

Options to L<test|test> available to all the sub-package tests are as follows.

=over 8

=item data => 'I<string>'

Optionally specify the name of the data to be tested. By default, this is not required: the testdata are those that were last loaded, either anonymously, or as given to one of the dichotomising methods. Otherwise, I<if the data are already ready for testing in a dichotomous format>, data that were previously loaded by name can be individually tested. For example, here are two sets of data that are loaded by name, and then a single test of one of them is performed.

 @chimps = (qw/banana banana cheese banana cheese banana banana banana/);
 @mice = (qw/banana cheese cheese cheese cheese cheese cheese cheese/);
 $seq->load(chimps => \@chimps, mice => \@mice);
 $seq->test('runs', data => 'chimps')->dump();

=item ccorr => I<boolean>

Specify whether or not to perform the continuity-correction on the observed deviation. Default is false. See L<Statistics::Zed>.

=item tails => I<1>|I<2>

Specify whether the I<z>-value is calculated for both sides of the normal distribution (2, the default for most testdata) or only one side (the default for data prepared with the B<swing> method).

=back

=head4 Test-specific required settings and options

Some sub-package tests need to have parameters defined in the call to L<test|test>, and/or have specific options, as follows.

B<Joins> : The Joins-test I<optionally> allows the setting of a probability value; see L<test|test> in the  L<Statistics::Sequences::Joins|Statistics::Sequences::Joins/test> manpage.

B<Pot> : The Pot-test I<requires> the setting of a state to be tested; see C<test> in the  L<Statistics::Sequences::Pot|Statistics::Sequences::Pot/test> manpage.

B<Runs> : There are presently no specific requirements nor options for the Runs-test. 

=head2 Series testing (provisional)

A means to aggregate results from multiple tests is exploratively supported, but only when accessing the tests through the main Sequences package. Three methods are presently used to effect this.

=head3 series_init

Clears any already accumulated data from previous tests.

=head3 series_update

Called once you have performed a test on a sample. It caches the observed, expectation and variance values from the test.

=head3 series_test

Sums the observed, expectation and variance values from all the tests updated to the series since calling L<series_init|series_init>, and produces a I<z_value> from these sums. It returns nothing in particular, but the following statement shows how the series values can be accessed.

 print "Series summed runs: 
    expected = ", $seq->{'series'}->{'expected'}, " 
    observed = ", $seq->{'series'}->{'observed'}," 
    z = $seq->{'series'}->{'z_value'}, $seq->{'tails'}-p = $seq->{'series'}->{'p_value'}\n";

=head2 Accessing results

All relevant statistical values are "lumped" into the class-object, and can be retrieved thus:

 $seq->{'observed'} # The observed value of the test-statistic (Runs, Joins, Pot)
 $seq->{'expected'} # The expected value of the test-statistic (Runs, Joins, Pot)
 $seq->{'obs_dev'} # The observed deviation (observed minus expected values), continuity-corrected, if so specified
 $seq->{'std_dev'} # The standard deviation
 $self->{'variance'} # Variance
 $seq->{'z_value'} # The value of the z-statistic (ratio of observed to standard deviation)
 $seq->{'p_value'} # The "normal probability" associated with the z-statistic

=head3 dump

 $seq->dump(flag => '1|0', text => '0|1|2', s_precision => 'integer', p_precision => 'integer');

I<Alias:> C<print_summary>

Print results of the last-conducted test to STDOUT. By default, if no parameters to C<dump> are passed, a single line of test statistics are printed. Options are as follows.

=over 8

=item flag => I<boolean>

If true, the I<p>-value associated with the I<z>-value is appended with a single asterisk if the value if below .05, and with two asterisks if it is below .01.

If false (default), nothing is appended to the p-value.

=item text => I<0>|I<1>|I<2>

If set to 1 (the default for an empty call to L<dump|dump>, a single line is printed, beginning with the name of the test, then the observed and expected values of the test-statistic, and the z-value and its associated p-value. The Pot-test, additionally, shows the state tested in parentheses after the test-name. For example:

 Joins: expected = 400.00, observed = 360.00, z = -2.83, p = 0.0023389**
 Runs: expected = 398.86, observed = 361.00, z = -2.70, p = 0.0070374**
 Pot(1): expected = 288.51, observed = 303.63, z = 2.64, p = 0.0082920**

If set to anything greater than 1, more verbose info is printed: each of the above bits of info are printed, on separate lines, as well as the observed and standard deviations.

If set to zero, no statistics are printed ... This was useful at one point in development, and might become so again ...

=item s_precision => 'I<non-negative integer>'

Precision of the z-statistic.

=item p_precision => 'I<non-negative integer>'

Specify rounding of the probability associated with the z-value to so many digits. If zero or undefined, you get everything available.

=back

=head3 dump_data

 $seq->dump_data(delim => "\n")

Prints to STDOUT a space-separated line of the testdata - as dichotomised and put to test. Optionally, give a value for I<delim> to specify how the datapoints should be separated.

=head3 string

 $seq->string()

Returns a single line giving the z_value and p_value. Accepts the I<s_precision>, I<p_precision> and I<flag> options, as for L<dump|dump>.

=head1 SEE ALSO

L<Statistics::Burst|Statistics::Burst> : Another test of sequences.

=head1 TO DO/BUGS

Results are dubious if there are only two observations.

Testing not by <z>-scores, and/or using poisson distribution for low number of observations

Multivariate extension to Runs test, at least

Sort option for pool method ?

=head1 REVISION HISTORY

=over 4

=item v 0.02

June 2008

Initital release via PAUSE.

See CHANGES in installation dist for revisions.

=back

=head1 AUTHOR

Roderick Garton, E<lt>rgarton@utas_DOT_edu_DOT_auE<gt>

=head1 COPYRIGHT/LICENSE/DISCLAIMER

Copyright (C) 2007-2008 Roderick Garton 

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself, either Perl version 5.8.8 or, at your option, any later version of Perl 5 you may have available. 

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=cut

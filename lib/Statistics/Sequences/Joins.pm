package Statistics::Sequences::Joins;

use 5.008008;
use strict;
use warnings;
use Carp 'croak';
use vars qw($VERSION @ISA);

use Statistics::Zed 0.02;
use Statistics::Sequences 0.041;
@ISA = qw(Statistics::Sequences);

$VERSION = '0.041';

=pod

=head1 NAME

Statistics::Sequences::Joins - the Wishart-Hirshfeld test of dichotomous sequences

=head1 SYNOPSIS

  use Statistics::Sequences::Joins;
  $joins = Statistics::Sequences::Joins->new();
  $joins->load(qw/0 0 1 0 0 0 1 0 0 0 0 1 1 0 0 0 1 1 1 1 0 0/);
  $joins->test()->dump();

=head1 DESCRIPTION

Superficially, this test is very similar to the L<Runs Test|Statistics::Sequences::Runs>. A I<join> is a point in a sequence of dichotomous data where the values alternate. For example, in the following series

 0 0 1 0 0 0 1 0 0 0 0 1 1 0 0 0 1 1 1 1 0 0

there is a join (of the values 0 and 1) at the indices of 1 and 2, then immediately another join (of the values 1 to 0) at 2 and 3, and then another join at 5 and 6. (Continuing, there are eight joins in total in this series.) The only difference between I<joins> and I<runs>, then, is that joins are counted at the point of alternation, whereas runs are counted for each unique segment (between the alternations). So joins will always be at least 1 less the number of runs.

The test-statistics, however, differ more fundamentally. Runs are tested on the basis of the observed distribution (e.g., of hits and misses), whereas joins are tested on the basis of a theoretically expected distribution.

=head1 METHODS

Methods are essentially as described in L<Statistics::Sequences>. See this manpage for how to handle non-dichotomous data, e.g., numerical data, or those with more than two categories; the relevant methods are not described here.

=head2 new

 $join = Statistics::Sequences::Joins->new();

Returns a new Joins object. Expects/accepts no arguments but the classname.

=head2 load

 $joins->load(@data);
 $joins->load(\@data);
 $joins->load('sample1' => \@data1, 'sample2' => \@data2)
 $joins->load({'sample1' => \@data1, 'sample2' => \@data2})

Loads data anonymously or by name. See L<load|Statistics::Sequences/load> in the Statistics::Sequences manpage.

=head2 test

 $joins->test(prob => 1/3);

Test the currently loaded data for significance of the number of joins, given a probability of each state in the data array of B<prob>. This parameter is optionally defined in the call to C<test>; the default value is 0.5. Valid values are between 0 and 1, inclusive. In a classic signal detection task, for instance, where a state is either a hit or miss, the probability of each state is 0.5. In a classic ESP task, with 5 possible alternatives as targets on each trial, each state has a 1 in 5 chance of occurring, but we still test for a probability of 1/2 as the data will be have to be reduced to a dichotomous hit/miss format in order to be Runs- or Joins-tested.

=cut

#-----------------------------------------------
sub test {
#-----------------------------------------------

   my $self = shift;
   my $args = ref $_[0] ? $_[0] : {@_};

   $self->_testdata_aref($args);
   
   my @data = @{$self->{'testdata'}};
   
   my $p = defined $args->{'prob'} ? delete $args->{'prob'} : 0.5;
   my $q = 1 - $p;

   my $w = defined $args->{'windows'} ? delete $args->{'windows'} : scalar @data;
   
   my ($obs_val, $exp_val, $var, $i, %states) = (0, 0);    

   # Count the number of joins:
   foreach ($i = 0; $i < @data; $i++) {
       $states{$data[$i]}++;
       if ($i and $data[$i] ne $data[$i - 1] ) {
           $obs_val++;
       }
   }

   return undef if scalar keys %states < 2;
   #croak __PACKAGE__, '::test Less than two states were found in the data' if scalar keys %states < 2;
   croak __PACKAGE__, '::test More than two states were found in the data: ' . join(' ', keys(%states)) if scalar keys %states > 2;

   $exp_val = 2 * ($w - 1) * $p * $q;
    
   $var = ( 4 * $w * $p * $q ) * (1 - ( 3 * $p * $q ) ) - ( ( 2 * $p * $q ) * (3 - ( 10 * $p * $q ) ) ); 
 
   if ($var) {
       $self->_expound($obs_val, $exp_val, $var, $args);
   }
   else {
       $self->_expire($obs_val, $exp_val, $args);
   }

   return $self;
}


=head2 dump

 $joins->dump(flag => '1|0', text => '0|1|2');

Print Joins-test results to STDOUT. See L<dump|Statistics::Sequences/dump> in the Statistics::Sequences manpage for details.

=cut

#-----------------------------------------------
sub dump {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $args->{'testname'} = 'Joins';
    $self->SUPER::_dump_pass($args);
    return $self;
 }

1;

__END__

=head1 EXAMPLE

Here the problem is to assess the degree of consistency of ESP scoring from the number of hits obtained in each of 200 runs of 25 trials each. The number of hits expected on the basis of chance is 5 per run. To test for sustained high or low scoring sequences, a join is defined as the point at which a score on one side of this expectation value is followed by a score on the other side. Ignoring scores equalling the expectation value of 5, the probability of a join is 1/2, or 0.5 (the default value to L<test|test>), assuming that, say, a score of 4 is as likely as a score of 6, and anything greater than a deviation of 5 (from 5) is improbable (or impossible). A meaningful result would obtain when the number of joins observed was below the number expected; we ignore the converse (and perverse) situation.

 use Statistics::Sequences;

 # Conduct 200 pseudo-ESP runs:
 my ($i, $hits, $target, $response, @scores);
 foreach ($i = 0; $i < 200; $i++) {
    $hits = 0;
    for (0 .. 24) {
        $target = (qw/star plus wave square circle/)[int(rand(5))];
        $response = (qw/star plus wave square circle/)[int(rand(5))];
        $hits++ if $target eq $response;
    }
    $scores[$i] = $hits;
  }

  my $expected_hits = 5;

  my $seq = Statistics::Sequences->new();
  $seq->load(@scores);
  $seq->cut(value => $expected_hits, equal => 0);
  $seq->test(what => 'joins', tails => 1, ccorr => 1)->dump(text => 1, flag => 1);
  # prints, e.g., Joins: expected = 79.00, observed = 67.00, z = -1.91, p = 0.028109*

=head1 REFERENCES

Burdick, D. S., & Kelly, E. F. (1977). Statistical methods in parapsychological research. In B. B. Wolman (Ed.), I<Handbook of Parapsychology> (pp. 81-130). New York, NY, US: Van Nostrand Reinhold.

Wishart, J. & Hirshfeld, H. O. (1936). A theorem concerning the distribution of joins between line segments. I<Journal of the London Mathematical Society>, I<11>, 227.

=head1 SEE ALSO

L<Statistics::Sequences::Runs|lib::Statistics::Sequences::Runs> : An analoguous and more widely known test. 

L<Statistics::Sequences::Pot|lib::Statistics::Sequences::Pot> : Another concept of sequences.

=head1 BUGS/LIMITATIONS

No computational bugs as yet identfied. Hopefully this will change, given time.

=head1 REVISION HISTORY

See CHANGES in installation dist for revisions.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2010 Roderick Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut

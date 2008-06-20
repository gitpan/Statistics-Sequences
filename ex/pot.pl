use Statistics::Sequences::Pot;
use strict;
   $|=1;
my $pot = Statistics::Sequences::Pot->new();

# Test the relative runs of a specific event (e.g., "5") among these data:
$pot->load([qw/2 0 8 5 3 5 2 3 1 1 9 4 4 1 5 5 6 5 8 7 5 3 8 5 6/]);
$pot->test(event => 5);
print "\n\nTEST 1\n";
$pot->dump(data => 1); # Print out the Pot-statistic, and a z-test of its significance:


# Init an array of random data with integers ranging from 0 to 15:
my @data = ();
my $i = 0;

for ($i = 0; $i < 960; $i++) {
  $data[$i] = int(rand(16));
}

$pot->load(\@data);

# Assess degree of bunching within these data with respect to a randomly selected target event:
my $event = int(rand(16));

$pot->test(event => $event);

# Access the results of this analysis:
$pot->dump(text => 2, flag => 1, data => 1);
  
# See what else was happening, having already given run() the data to test:
foreach (0 .. 15) {
	next if $_ == $event;
    $pot->test(event => $_)->dump(flag => 1);
}

  # Example 3
  print "\n\nTEST 3\n";
  my $i;
  # Init an array of random data with values of either 'hit' or 'miss':
  my @categories = (qw/hit miss/);  
  my @data = ();
  for ($i = 0; $i < 640; $i++) {
    $data[$i] = $categories[int(rand(@categories))];
  }

  # Make a pot object and load up the data:
  $pot = Statistics::Sequences::Pot->new();
  $pot->load(\@data);

  # Run and dump a couple analyses, on each possible event:
  $pot->test(event => 'hit');
  $pot->dump(); # prints, e.g.: Event hit: Pot = 250.07, z = 1.12, p = 0.13142
  $pot->test(event => 'miss');
  $pot->dump(); # prints, e.g.: Event miss: Pot = 252.45, z = 1.16, p = 0.12373

  $pot->test(event => $categories[int(rand(@categories))], full => 1);

  print "Randomly selected event was a $pot->{'event'}, and this occurred $pot->{'n_events'} times, " .
        "most frequently bunching by a length of " . $pot->{'bunches'}->mode() . "\n";

__END__


  use Statistics::Sequences::Pot;
  use strict;
  my @data = ();
  my $total = 960;
  my $z = 0;
  
  # Init an array of random data with integer values ranging from 0 to 15:
  for ($z = 0; $z < $total; $z++) {
    $data[$z] = int(rand(16));
  }

  # Analyse the degree of bunching within these data with respect to a single target event:
  my $event = int(rand(16));
  my $pot = Statistics::Sequences::Pot->new();
  $pot->run(data => \@data, event => $event);
  
  # Access the results of this analysis:
  print "The probability of obtaining as much bunching of $event as observed is $pot->{'p_normal'}\n";
  # or:
  print "For event $pot->{event} occurring $pot->{n_events} times among $pot->{'n_units'},".
     "the observed value of Pot was $pot->{'observed'}.\n" .
     "The expected value of Pot was $pot->{'expected'}\n" .
     "The deviation from expectation was $pot->{'deviation'}\n" .
     "The standard deviation was $pot->{'sd'}\n" .
     "Z = $pot->{'z'}\n" .
     "Probability of this deviation = $pot->{'p_normal'}\n";
  # or print the lot, and more, in English:
  $pot->dump(text => 2);
  
  # Now let's briefly see what else was happening, having already given run() the data to test:
  foreach (0 .. 15) {
  	next if $_ == $event;
	$pot->run(event => $_)->dump();
  }
   
   __END__

__END__
use Statistics::Sequences::Pot;
use Math::Random::Handy;
my $rnd = Math::Random::Handy->new();
my $pcqng = $rnd->pcqng();
my @data = ();
my $z = 0;
my $lim = 6;
#for ($z = 0; $z < $M; $z++) {
#   $data[$z] = $rnd->randomate_int_in($pcqng->RandUniform, 0, 15);
#}

#$pot = Statistics::Sequences::Pot->new;
#$pot->run(\@data, 7);
#$pot->print_summary();
#use Math::Random::MT qw(srand rand);
my $pot = Statistics::Sequences::Pot->new;
my $total = 960;
#srand;
#my $event = int(rand($lim));
my $event = $rnd->randomate_int_in($pcqng->RandUniform, 0, $lim);
print "event = $event\n";
$| = 1;
my @categories = ('hit', 'miss');
my $s = scalar(@categories);
#print "s = $s\n";
  # Build up an array of random data with values of either 'hit' or 'miss':
open F, ">pot2.dat";
#foreach (0 .. 9999) {
  for ($z = 0; $z < $total; $z++) {
	my $n = $rnd->randomate_int_in($pcqng->RandUniform, 0, $lim);
	#my $n = int(rand(@categories));
    #my $n =int(rand($lim));
	#my $n = random_int_in(0, $lim);
	#my $d = $categories[$n];
	#print "d = n = $n $d\n";
	$data[$z] = $n;
	print F $data[$z] . "\n";
  }
  $pot->load(\@data);
  $pot->run(event => $event, cc => 2, scale => 1, full => 1);
  $pot->dump(text => 2);#
  foreach (0 .. $lim) {
  	next if $_ == $event;
	$pot->run(event => $_, cc => 2, scale => 1, full => 1)->dump();
  }
  #  $pot->run(event => 'miss', cc => 1, scale => 1, full => 1);
  #$pot->dump(text => 2);#


#  print F "$pot->{'z'}\t$pot->{'p_normal'}\n";
#}
close F;

 sub random_int_in ($$) {
     my($min, $max) = @_;
      # Assumes that the two arguments are integers themselves!
     return $min if $min == $max;
     ($min, $max) = ($max, $min)  if  $min > $max;
     return $min + int rand(1 + $max - $min);
   }


__END__
N = 44
The observed Pot = 47.65
 The expected Pot is 40.63 
The deviation of Pot from chance is 7.02727904186796 
The standard deviation is 4.2186910373977 which gives z = 1.6657486835544
P = 0.047882




use constant two_pi_sqrt_inverse => 1 / sqrt(8 * atan2(1, 1));
sub gaussian {
	my ($x, $mean, $variance) = @_;
	return two_pi_sqrt_inverse * exp( -( ($x - $mean) ** 2) / (2 * $variance) ) / sqrt($variance);
}
	
	my $prob2 = gaussian($Pot, $Pot0, ($Upper/$Lower));
	print  "prob 2 = $prob2\n";
	


						#'with ' . 
						#    ($dev > 0 ? 'wider' : 'narrower' )
						#. ' spacing, ' .	
						

						

__END__

	
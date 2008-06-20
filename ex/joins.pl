 use Statistics::Sequences;
use Math::Random::Handy;

my $rnd = Math::Random::Handy->new();
my $fn = $rnd->init_rnd_src('mt');

 # Conduct 200 pseudo-ESP runs:
 my ($i, $hits, $target, $response, @scores, $above);
 foreach ($i = 0; $i < 1000; $i++) {
    $hits = 0;
    for (0 .. 24) {
        $target = (qw/star plus wave square circle/)[int(rand(5))];
        $response = (qw/star plus wave square circle/)[$fn->(0, 4)];
        $hits++ if $target eq $response;
    }
    $scores[$i] = $hits;#($hits - (int(rand(25))));
    $above++ if $hits >= 10;
  }
  
  my $expected_hits = 5;
  print "above or equal 10 = $above\n";
  #print join(', ', @scores), "\n";
  
  my $seq = Statistics::Sequences->new();
  $seq->load(@scores);
  $seq->cut(at => $expected_hits, equal => 0);
  #$seq->updown();
  $seq->test('joins', tails => 1)->dump(data => 0, text => 1, flag => 1);
  $seq->test('runs')->dump(data => 0, text => 1, flag => 1);
  $seq->test('pot', event => 1)->dump(data => 0, text => 1, flag => 1);
  
  #my @dat = ();
  #for (0 .. 200) {
  #      push @dat, (qw/star plus wave square circle/)[int(rand(5))];
  #}
  #$seq->load(@dat);
  #$seq->test('joins', prob => 1/5, tails => 2)->dump(data => 0, text => 1, flag => 1);
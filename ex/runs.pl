use Statistics::Sequences;
use Statistics::Sequences::Joins;
use Statistics::Sequences::Runs;
use Statistics::Sequences::Pot;

my $runs = Statistics::Sequences::Runs->new();
my $pot = Statistics::Sequences::Pot->new();
my $joins = Statistics::Sequences::Joins->new();

$| = 1;

my ($i, @targs, @resps, @scores);

use Statistics::Sequences::Runs;
$runs = Statistics::Sequences::Runs->new();
# $runs->load(qw/2 0 8 5 3 5 2 3 1 1/);
$runs->load(qw/1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 1 1 0 1/)->test()->dump();

for ($i = 0; $i < 25; $i++) {
    $targs[$i] = (qw/a b c d e/)[int(rand(5))];
    $resps[$i] = (qw/a b c d e/)[int(rand(5))];
    $scores[$i] = (qw/1 2 3/)[int(rand(3))];
}

print "\n" . '-=' x 20 . "\n";
print "Anonymous load and cut test\n";
#print '-=' x 20 . "\n";
$runs->load(\@scores);
$runs->cut(at => 'mean', equal => 0);
$runs->test()->dump();
print '-=' x 20 . "\n\n";

print "\n" . '-=' x 20 . "\n";
print "Nominal load and match test\n";
#print '-=' x 20 . "\n";
#$runs->load();
$runs->load(targs => \@targs, resps => \@resps);
$runs->match('data' => ['targs', 'resps'], lag => 1, loop => 1); # match 2 samples as hits/misses
$runs->test()->dump();

$joins->load(targs => \@targs, resps => \@resps);
$joins->match(data => ['targs', 'resps']); # match 2 samples
$joins->test(prob => 1/5)->dump();#, windows => 100

#$joins->load();
#$joins->load();
$joins->load($runs->{'testdata'});
$joins->test(prob => 1/5)->dump();#, windows => 100

$pot->load($runs->{'testdata'});
$pot->test(event => 1)->dump(data => 1);

print "\n" . '-=' x 20 . "\n";

 my @scores = (0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0);
 my $runs = Statistics::Sequences::Runs->new();
 $runs->load(\@scores);
 $runs->test();
 print " probability of obtaining these $runs->{'observed'} runs in 25 trials is $runs->{'p_value'}\n";
  
 my ($i, @targets, @responses);
 for ($i = 0; $i < 25; $i++) {
    $targets[$i] = (qw/star plus wave square circle/)[int(rand(5))];
    $responses[$i] = (qw/star plus wave square circle/)[int(rand(5))];
}
 
 my $runs = Statistics::Sequences::Runs->new();
 $runs->load(targets => \@targets, responses => \@responses);
 $runs->match(data => [qw/targets responses/]);
 $runs->test()->dump(text => 2, data => 1, flag => 1);
 print " probability of obtaining $runs->{'observed'} runs is $runs->{'p_value'}\n";
    
 # But what if the responses were actually guessed for the target on the trial on ahead?
 $runs->match(data => [qw/targets responses/], lag => 1)->test()->dump(data => 1);
 print "With precognitive responses to the target for the next trial,\n$runs->{'observed'} runs in 24 trials were produced when $runs->{'expected'} were expected,\nfor which the probability is $runs->{'p_value'}\n"; 
    
1;
use strict;
use Statistics::Sequences;
$|=1;
my ($i, @targs, @resps, @blues, @reds);

for ($i = 0; $i < 25; $i++) {
    $targs[$i] = (qw/a b c d e/)[int(rand(5))];
    $resps[$i] = (qw/a b c d e/)[int(rand(5))];
    $blues[$i] = int(rand(26)) + 10;
    $reds[$i] = int(rand(26)) + 10;
}

print "One sequence\n";
my $seq = Statistics::Sequences->new();
$seq->load([1, 1, 1,1,1,1,1,1]);#0, 0, 1, 0, 1, 1, 0, 0, 0, 1
$seq->test('runs')->dump(text => 1, flag => 1);
$seq->test('joins', prob => 1)->dump(text => 1);
$seq->test('pot', event => 1)->dump(text => 1);
#print "targets =\t", @targs, "\nresponses =\t", @resps, "\n";

#$seq->load(@scores);
#$seq->cut(point => 'median');
print "\nTargets and responses matched\n";
$seq->load({'targs' => \@targs, 'resps' => \@resps});
$seq->match(data => ['targs', 'resps'], lag => 0, loop => 0); # match 2 samples as hits/misses
$seq->test('runs')->dump(text => 1);
$seq->test('joins', prob => 1/5)->dump(data => 0, text => 1, flag => 1);
$seq->test('pot', event => 1)->dump(text => 1);
$seq->test('pot', data => 'resps', event => 'a')->dump(text => 1);
#$seq->match('data' => ['targs', 'resps'], lag => 2, loop => 1); # match 2 samples as hits/misses
#$seq->test('runs', dump => 1);

#$seq->dump();
#$seq->pool('data' => ['targs', 'resps'], lag => -1, loop => 0); # match 2 samples as hits/misses

print "\nChimps and Mice\n";
 my @chimps = (qw/banana banana cheese banana cheese banana banana banana/);
 my @mice = (qw/banana cheese cheese cheese cheese cheese cheese cheese/);
 $seq->load(chimps => \@chimps, mice => \@mice);
 $seq->test('joins', data => 'chimps', prob => 1/2)->dump(text => 1);
 $seq->test('joins', data => 'mice', prob => 1/2)->dump();
 
print "\nBlues and Reds Pooled\n";
 $seq->load(blues => \@blues, reds => \@reds);
 $seq->pool(data => [qw/blues reds/]);
 $seq->test('runs', tails => 2)->dump(data => 0, text => 1);
 $seq->test('pot', event => 'blues')->dump(data => 0, text => 1);
 $seq->test('pot', event => 'reds')->dump(data => 0, text => 1);
 
print "\nSwing Blues\n";
#$seq->load(blues => \@blues);
$seq->swing(data => 'blues', equal => 0);
$seq->test('runs')->dump(data => 1, text => 1);
$seq->test('pot', event => 1, dist => 'norm', tails => 1)->dump(data => 0, text => 1);

my @dt = (qw/1 1 1 1 1 0 0 0 0 0/);
$seq->load(\@dt);
$seq->test('runs')->dump();

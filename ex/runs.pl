my ($i, @targs, @resps, @scores);

use Statistics::Sequences::Runs;
my $runs = Statistics::Sequences::Runs->new();

for ($i = 0; $i < 100; $i++) {
    $targs[$i] = (qw/a b c d e/)[int(rand(5))];
    $resps[$i] = (qw/a b c d e/)[int(rand(5))];
    $scores[$i] = (qw/1 2 3/)[int(rand(3))];
}

print "\n" . '-=' x 20 . "\n";
print "ANONYMOUS LOAD & CUT TEST\n";
print "- a 100-sample list of scores made up of 3 states - the digits 1, 2, 3 - has been randomly generated\n";
print "- they will now be tested for their runs about the mean\n\n";
print "\nuse Statistics::Sequences::Runs;\n";
print "my \$runs = Statistics::Sequences::Runs->new();\n";
print "\$runs->load([\@scores]);\n";
print "\$runs->cut(value => 'mean', equal => 0);\n";
print "\$runs->test()->dump();\n";
print "\n-This leads to the following print to STDOUT:\n";
$runs->load(\@scores);
$runs->cut(value => 'mean', equal => 0);
$runs->test()->dump();
print "\n-But what, actually, was the mean?\n";
print "print \$runs->{'cut_value'};\n";
print $runs->{'cut_value'}, "\n";
print '-=' x 20 . "\n\n";

print "\n" . '-=' x 20 . "\n";
print "NOMINAL LOAD & MATCH TEST\n";
print "- 2 datasets of 100 samplings each have been generated\n";
print "- they each contain an independent sampling of one of 5 alternative states\n";
print "- they will now be matched for the synchrony of their states at each sampling event:\n";
print "\n";
print "\$runs->load(targs => [\@targs], resps => [\@resps]);\n";
print "\$runs->match('data' => ['targs', 'resps']); # match 2 samples as hits/misses\n";
print "\$runs->test()->dump();\n";

$runs->load(targs => \@targs, resps => \@resps);
$runs->match('data' => ['targs', 'resps']); # match 2 samples as hits/misses
$runs->test()->dump();

print "\n\n-Do it again, but this time lag the response data, so each of its states refers to the state on the next target sample:\n";
print "\n";
print "\$runs->load(targs => [\@targs], resps => [\@resps]);\n";
print "\$runs->match('data' => ['targs', 'resps'], lag => 1, loop => 0);\n";
print "\$runs->test()->dump();\n";

$runs->load(targs => \@targs, resps => \@resps);
$runs->match('data' => ['targs', 'resps'], lag => 1, loop => 0); # match 2 samples as hits/misses
$runs->test()->dump();

print "\n- ... and what would the Pot test do with the same data?\n";
print "\n";
print "use Statistics::Sequences::Pot;\n";
print "my \$pot = Statistics::Sequences::Pot->new();\n";
print "\$pot->load(\$runs->{'testdata'});\n";
print "\$pot->test(state => 1)->dump();\n";
require Statistics::Sequences::Pot;
my $pot = Statistics::Sequences::Pot->new();
$pot->load($runs->{'testdata'});
$pot->test(state => 1)->dump();
print "\n";
print "\n" . '-=' x 20 . "\n";

1;
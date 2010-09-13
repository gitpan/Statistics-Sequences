use Statistics::Sequences::Runs;
my $runs = Statistics::Sequences::Runs->new();

print "CHECKING the output of Statistics::Sequences with published examples\n";
print "(Note: all z-scores are calculated with continuity correction, i.e., ccorr => 1)\n";
print "Examples from Swed and Eisenhart (1943):\n\n";

print "1. TESTS of RANDOMNESS\n";

print "(a) ROW of PLANTS\n";
print "-The following represents the healthy (H) vs. diseased (D) state of a sample of plants grown in a row:\n";
my @dat = (qw/H H H H H H H H H H D H D D D D H H H H H H H H H/);
print @dat, "\n";
print "-Swed and Eisenhart report: Runs observed = 5, p = .0183512\n-What does the Runs module produce?\n";
$runs->load(\@dat);
$runs->test(ccorr => 1, tails => 1)->dump(text => 1);
print "With fewer runs than expected, it looks like disease is local and not randomly distributed over the sample.\n";
print "\n\n";

print "(b) ROW of SEATS\n";
print "Swed and Eisenhart provide another example for occupied (O) and empty (E) seats in a row at a lunch counter.\nHave people taken up their seats on a random basis?\n";
@dat = (qw/E O E E O E E E O E E E O E O E/);
print @dat, "\n";
print "They report: Runs = 11, p = .057\n-What does the Runs module produce?\n";
$runs->load(\@dat);
$runs->test(ccorr => 1, tails => 1)->dump(text => 1);
print "\nA customer sits at the fifth empty seat ...\n";
@dat = (qw/E O E E O E O E O E E E O E O E/);
print @dat, "\n";
print "They report: Runs = 13, p = .010\n-What does the Runs module produce?\n";
$runs->load(\@dat);
$runs->test(ccorr => 1, tails => 1)->dump(text => 1);
print "They concluded: \"Both of these cases exhibit too many groups to be considered random arrangements\",\ni.e., there are too many alternations in the seating sequences than would be expected if they were produced by a random process.\n";
print "\n\n";

print "2. TESTS of difference between two samples\n";

print "(a) Gains of CALVES on different food rations\n";
my @calves_a = (qw/1.95 2.17 2.06 2.11 2.24 2.52 2.04 1.95/);
my @calves_b = (qw/1.82 1.85 1.87 1.74 2.04 1.78 1.76 1.86/);
print "Calves A: ", join(" ", @calves_a), "\n";
print "Calves B: ", join(" ", @calves_b), "\n";
print "Reported number of runs by pooling: 4, with p = .0088578\nWhat does Runs module yield?\n";
$runs->load(Ca => \@calves_a, Cb => \@calves_b);
$runs->pool(data => [qw/Ca Cb/]);
#print "testdata: ", join(" ", @{$runs->{'testdata'}}), "\n";
$runs->test(ccorr => 1, tails => 1)->dump(text => 2);
print "Look at the data, as ordered, to see what's happening (Ca = Calves A, Cb = Calves B):\n";
$runs->dump_data(delim => ", ");
print "So it looks like Calves B are consistently gaining less on their ration than the portly Calves A.\nIf the calves have been randomly assigned to these rations, it looks like the ration given to the A's\nis a better gainer than the other.";
print "\n\n";

print "(b) Attidue scores between unmarried and married persons (example from Sarantakos (1993, p. 396))\n";
my @unmarried = (qw/14 13 11 9 9 8 8 5/);
my @married = (qw/12 10 10 7 7 15 6 6/);
print "unmarried: ", join(" ", @unmarried), "\n";
print "married: ", join(" ", @married), "\n";
print "Reported number of runs by pooling: 8; and table of critical values reveals no sig. diff.\nWhat does the Runs module yield?\n";
$runs->load(u => \@unmarried, m => \@married);
$runs->pool(data => [qw/u m/]);
$runs->test(ccorr => 1, tails => 2)->dump(text => 2);
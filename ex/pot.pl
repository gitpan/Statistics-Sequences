use Statistics::Sequences::Pot;
use strict;
my $pot = Statistics::Sequences::Pot->new();
$|=1;

print "\n\nTEST 1\n";
print "-Init an array of random data with integers ranging from 0 to 15 ...\n";
my @data = ();
my $i = 0;

for ($i = 0; $i < 960; $i++) {
  $data[$i] = int(rand(16));
}
print "-load these into the pot object like so:\n\n";
print "use Statistics::Sequences::Pot;\n";
print "use strict;\n";
print "my \$pot = Statistics::Sequences::Pot->new();\n";
print "\$pot->load([\@data]);\n\n";
$pot->load(\@data);

print "-Assess degree of bunching within these data with respect to a randomly selected target state ...\n\n";
print "my \$state = int(rand(16));\n";
print "\$pot->test(state => \$state);\n\n";
my $state = int(rand(16));
$pot->test(state => $state);

print "-Access the results of this analysis:\n\n";
print "\$pot->dump(text => 2, flag => 1);\n\n";
$pot->dump(text => 2, flag => 1);
  
print "\n-See what else was happening, having already given run() the data to test:\n\n";
print "foreach (0 .. 15) {\n";
print "\tnext if \$_ == \$state;\n";
print "\t\$pot->test(state => \$_)->dump(flag => 1);\n";
print "}\n";
print "\n-This gets the following printed to SDTOUT:\n\n";
foreach (0 .. 15) {
	next if $_ == $state;
    $pot->test(state => $_)->dump(flag => 1);
    print "\n";
}
	
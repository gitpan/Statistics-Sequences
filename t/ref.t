use strict;
use warnings;
use Test::More tests => 6;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::Sequences') };

my $seq = new_ok('Statistics::Sequences');

my %refdat = (
    test1 => { observed => 5, expected => 9, z_value => -2.29, p_value => .010973},
    test2 => { observed => 4, expected => 9, z_value => -2.33, p_value => .009930},
);

my @dat = (qw/H H H H H H H H H H D H D D D D H H H H H H H H H/);

# Load by array-ref:
eval { $seq->load(\@dat);};
ok(!$@);

eval {$seq->unload();};
ok(!$@);

# check minimal integration with Statistics::Data parent's load(), add() and read() methods:
$seq->load(coinflip => \@dat, otherdat => [1, 2, 3]);
my $data = $seq->read(label => 'coinflip');
ok(join('', @$data) eq join('',@dat), "Failed to read data");
$seq->add(otherdat => [4]);
$data = $seq->read(label => 'otherdat');
my $sum = 0;
$sum += $_ foreach @$data;
ok($sum == 10, "Failed to read data");

sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}

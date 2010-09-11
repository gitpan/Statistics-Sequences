use strict;
use warnings;
use Test::More tests => 61;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::Sequences') };

my $seq = Statistics::Sequences->new(
	ccorr => 1,
    tails => 1,
  	precision_s => 2,
    precision_p => 6,
);
isa_ok($seq, 'Statistics::Sequences');

my %refdat = (
    test1 => { observed => 5, expected => 9, z_value => -2.29, p_value => .010973},
    test2 => { observed => 4, expected => 9, z_value => -2.33, p_value => .009930},
    chimps => {observed  => 4, expected => 3.50, z_value => 0, p_value => 1.00000},
    mice => {observed  => 1, expected => 3.50, z_value => -1.512, p_value => 0.13057},
);

my @dat = (qw/H H H H H H H H H H D H D D D D H H H H H H H H H/);

# Load by array-ref:
eval { $seq->load(\@dat);};
ok(!$@);
eval { $seq->test(what => 'runs', tails => 1);};
ok(!$@);

foreach (qw/observed expected z_value p_value/) {
    ok(defined $seq->{$_} );
    ok(equal($seq->{$_}, $refdat{'test1'}->{$_}), "$_  $seq->{$_} = $refdat{'test1'}->{$_}");
}

eval {$seq->unload();};
ok(!$@);

# Load by array:
eval { $seq->load(@dat);};
ok(!$@);
eval { $seq->test(what => 'runs', tails => 1);};
ok(!$@);

foreach (qw/observed expected z_value p_value/) {
    ok(defined $seq->{$_} );
    ok(equal($seq->{$_}, $refdat{'test1'}->{$_}), "$_  $seq->{$_} = $refdat{'test1'}->{$_}");
}

my @calves_a = (qw/1.95 2.17 2.06 2.11 2.24 2.52 2.04 1.95/);
my @calves_b = (qw/1.82 1.85 1.87 1.74 2.04 1.78 1.76 1.86/);

# Load by hash:
eval { $seq->load(Ca => \@calves_a, Cb => \@calves_b);};
ok(!$@);
eval { $seq->pool(data => [qw/Ca Cb/]);};
ok(!$@);
eval { $seq->test(what => 'runs', tails => 1);};
ok(!$@);

foreach (qw/observed expected z_value p_value/) {
    ok(defined $seq->{$_} );
    ok(equal($seq->{$_}, $refdat{'test2'}->{$_}), "$_  $seq->{$_} = $refdat{'test2'}->{$_}");
}

$seq->unload();

# Load by hash-ref:
eval { $seq->load({Ca => \@calves_a, Cb => \@calves_b});};
ok(!$@);
eval { $seq->pool(data => [qw/Ca Cb/]);};
ok(!$@);
eval { $seq->test(what => 'runs');};
ok(!$@);

foreach (qw/observed expected z_value p_value/) {
    ok(defined $seq->{$_} );
    ok(equal($seq->{$_}, $refdat{'test2'}->{$_}), "$_  $seq->{$_} = $refdat{'test2'}->{$_}");
}

# Access data-series by name for unique stats (using joins), also changing s_precision:
my @chimps = (qw/banana banana cheese banana cheese banana banana banana/);
my @mice = (qw/banana cheese cheese cheese cheese cheese cheese cheese/);
$seq->load({chimps => \@chimps, mice => \@mice});
$seq->match(data => ['chimps', 'mice']);
$seq->test(what => 'joins', data => 'chimps', prob => 1/2, precision_s => 3, tails => 2);
foreach (qw/observed expected z_value p_value/) {
   ok(defined $seq->{$_} );
   ok(equal($seq->{$_}, $refdat{'chimps'}->{$_}), "$_  $seq->{$_} = $refdat{'chimps'}->{$_}");
}
$seq->test(what => 'joins', data => 'mice', prob => 1/2, precision_s => 3, tails => 2);
foreach (qw/observed expected z_value p_value/) {
   ok(defined $seq->{$_} );
   ok(equal($seq->{$_}, $refdat{'mice'}->{$_}), "$_  $seq->{$_} = $refdat{'mice'}->{$_}");
}

sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}

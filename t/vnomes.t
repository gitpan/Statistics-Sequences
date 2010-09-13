use strict;
use warnings;
use Test::More tests => 12;
use constant EPS => 1e-3;

BEGIN { use_ok('Statistics::Sequences::Vnomes') };

my $seq = Statistics::Sequences::Vnomes->new();
isa_ok($seq, 'Statistics::Sequences::Vnomes');

my %refdat = (
    nist => {
        psisq => 0.8,
        p_value => 0.67032,
    },
	gatlin => {
        psisq => 36.2909090909091,
        p_value => .45509,}
);

# Gatlin data:
my @data = (qw/G A A T A C A G C C T G T C G G T T C T C C G A T G G C A A G T A C T T T A C T G G T T C A G A A T G C G C C C G T A A T C C T T G C A C A G A G G A T C T T A C G C A G T G A A C C G G C T C C G T G G A T T A G C A A C/);

eval {
    $seq->load_data(\@data);
};
ok(!$@, $@);

$seq->test(length => 3, delta => 2, circularize => 1, states => [qw/A C G T/], precision_s => 3);

foreach (qw/psisq p_value/) {
   ok(defined $seq->{$_} );
   ok(about_equal($seq->{$_}, $refdat{'gatlin'}->{$_}), "$_  $seq->{$_} = $refdat{'gatlin'}->{$_}");
}

# NIST data:
@data = (qw/0 0 1 1 0 1 1 1 0 1/);

eval {
    $seq->load_data(\@data);
};
ok(!$@, $@);

$seq->test(length => 3, delta => 2, states => [0, 1], precision_s => 3);

foreach (qw/psisq p_value/) {
   ok(defined $seq->{$_} );
   ok(about_equal($seq->{$_}, $refdat{'nist'}->{$_}), "$_  $seq->{$_} = $refdat{'nist'}->{$_}");
}

sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}

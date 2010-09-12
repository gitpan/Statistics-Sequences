use strict;
use warnings;
use Test::More tests => 15;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::Sequences::Turns') };

my $seq = Statistics::Sequences::Turns->new();
isa_ok($seq, 'Statistics::Sequences::Turns');

my %refdat = (
        std_dev => 3.04,
        z_value => -0.054717,
        p_value => 0.95636,
        variance => 9.278,
        observed => 35,
        expected => 34.667,
);

# Gatlin data:
my @data = (15.2, 16.9, 15.3, 14.9, 15.7, 15.1, 16.7, 16.3, 16.5, 13.3, 16.5, 15.0, 15.9, 15.5, 16.9, 16.4, 14.9, 14.5, 16.6, 15.1, 14.6, 16.0, 16.8, 16.8, 15.5, 17.3, 15.5, 15.5, 14.2, 15.8, 15.7, 14.1, 14.8, 14.4, 15.6, 13.9, 14.7, 14.3, 14.0, 14.5, 15.4, 15.3, 16.0, 16.4, 17.2, 17.8, 14.4, 15.0, 16.0, 16.8, 16.9, 16.6, 16.2, 14.0, 18.1, 17.5);

eval {
    $seq->load_data(\@data);
};
ok(!$@, $@);

$seq->test(precision_s => 3, tails => 2);

foreach (qw/observed expected variance std_dev z_value p_value/) {
   ok(defined $seq->{$_} );
   ok(about_equal($seq->{$_}, $refdat{$_}), "$_  $seq->{$_} = $refdat{$_}");
}

sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}

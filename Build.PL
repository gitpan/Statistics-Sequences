use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Sequences',
    license             => 'perl',
    dist_author         => 'Roderick Garton <rgarton@cpan.org>',
    dist_version_from   => 'lib/Statistics/Sequences.pm',
    build_requires => {
        'Test::More' => 0,
        'Class::OOorNo' => 0,
        'Statistics::Lite' => 3.2,
        'Statistics::Zed' => 0.01,
        'Statistics::Distributions' => 1.02,
        'Statistics::ChisqIndep' => 0.1,
	'Algorithm::Combinatorics' => '0.23',
	'Scalar::Util' => '1.18',
    },
    add_to_cleanup      => [ 'Statistics-Sequences-0.031-*' ],
);

$builder->create_build_script();

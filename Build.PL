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
        'Statistics::Lite' => 3.2,
        'Statistics::Zed' => 0.04,
        'Statistics::Distributions' => 1.02,
    	'Algorithm::Combinatorics' => '0.23',
	'Scalar::Util' => '1.18',
        'Math::Cephes' => '0.43',
    },
    add_to_cleanup      => [ 'Statistics-Sequences-0.051-*' ],
);

$builder->create_build_script();

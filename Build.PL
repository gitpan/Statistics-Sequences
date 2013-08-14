use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Sequences',
    license             => 'perl',
    dist_author         => 'Roderick Garton <rgarton@cpan.org>',
    dist_version_from   => 'lib/Statistics/Sequences.pm',
    dist_abstract       => 'Manage sequences (ordered list of literals) for testing their runs, joins, turns, trinomes, potential energy, etc.',
    configure_requires => {
        'Module::Build' => '0.4',
    },
    build_requires => {
        'Test::More' => 0,
        'Statistics::Data' => '0.08',
        'Scalar::Util' => '1.18',
    },
    requires => {
        'Statistics::Data' => '0.08',
        'Scalar::Util' => '1.18',
    },
    add_to_cleanup      => [ 'Statistics-Sequences-0.12-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();

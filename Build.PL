use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Sequences',
    license             => 'perl',
    dist_author         => 'Roderick Garton <rgarton@cpan.org>',
    dist_version_from   => 'lib/Statistics/Sequences.pm',
    configure_requires => {
        'Module::Build' => '0.36',
    },
    build_requires => {
        'Test::More' => 0,
        'Statistics::Data' => '0.02',
    },
    add_to_cleanup      => [ 'Statistics-Sequences-0.10-*' ],
);

$builder->create_build_script();

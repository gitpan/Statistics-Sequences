use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Statistics::Sequences',
    AUTHOR              => 'Roderick Garton <rgarton@cpan.org>',
    VERSION_FROM        => 'lib/Statistics/Sequences.pm',
    ABSTRACT_FROM       => 'lib/Statistics/Sequences.pm',
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
        'Statistics::Data' => '0.02',
        'Scalar::Util' => '1.18',
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Statistics-Sequences-0.11-*' },
);

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
        'Statistics::Lite' => 3.2,
        'Statistics::Zed' => 0.01,
        'Statistics::Distributions' => 1.02,
    	'Algorithm::Combinatorics' => '0.23',
	    'Scalar::Util' => '1.18',
        'Math::Cephes' => '0.43',
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Statistics-Sequences-0.042-*' },
);

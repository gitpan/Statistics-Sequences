use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Statistics::Sequences',
    VERSION_FROM      => 'lib/Statistics/Sequences.pm', # finds $VERSION
    PREREQ_PM         => {'Class::OOorNo' => 0, 'Statistics::Lite' => 0, 'Statistics::Deviation' => 0, 'Statistics::Distributions' => 0},
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Statistics/Sequences.pm', # retrieve abstract from module
       AUTHOR         => 'rgarton <rgarton@cpan.org>') : ()),
);

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Statistics-Sequences.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
BEGIN { 
    use_ok('Statistics::Sequences')
};

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

use Statistics::Basic::Correlation;

plan tests => 5;

my  $corr = new Statistics::Basic::Correlation([1..10], [1..10]);

ok( $corr->query == 1.0 );

    $corr->insert( 11, 7 );
ok( $corr->query == ( (129/20) / (sqrt(609/100) * sqrt(165/20))));

    $corr->set_vector( [11..13], [11..13] );
ok( $corr->query == 1.0 );

    $corr->ginsert( 13, 12 );
ok( $corr->query == ( (1/2) / (sqrt(11/16) * sqrt(1/2)) ));

my  $j = new Statistics::Basic::Correlation;
    $j->set_vector( [11..13], [11..13] );
ok( $j->query == 1.0 );

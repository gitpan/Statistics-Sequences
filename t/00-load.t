#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Statistics::Sequences' );
}

diag( "Testing Statistics::Sequences $Statistics::Sequences::VERSION, Perl $], $^X" );

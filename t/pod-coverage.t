#!perl -T

use Test::More;
eval "use Test::Pod::Coverage 1.04";
plan skip_all => "Test::Pod::Coverage 1.04 required for testing POD coverage" if $@;
all_pod_coverage_ok({trustme => ['add_data', 'append', 'append_data', 'clear_data', 'delete_data', 'get_data', 'load_data', 'print_summary', 'process']});


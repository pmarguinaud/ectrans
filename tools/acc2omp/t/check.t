#!/usr/bin/env perl
# Tests of acc2omp.pl --check; run with: prove tools/acc2omp/t
use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);

my $tool = "$Bin/../acc2omp.pl";
chdir($Bin) or die;

my $out = `$^X $tool --check check_subset.F90 2>&1`;
is($?, 0, 'subset: exit status 0');
is($out, '', 'subset: no error');

$out = `$^X $tool --check --pairs check_subset.F90 2>&1`;
isnt($?, 0, 'subset with --pairs: KERNELS inside #ifdef ACCGPU is checked');

$out = `$^X $tool --check check_errors.F90 2>&1`;
isnt($?, 0, 'errors: exit status not 0');
open(my $fh, '<', 'check_errors.ref') or die;
my $ref = do { local $/; <$fh> };
is($out, $ref, 'errors: messages');

done_testing();

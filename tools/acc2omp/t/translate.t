#!/usr/bin/env perl
# Tests of acc2omp.pl translation mode; run with: prove tools/acc2omp/t
use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);

my $tool = "$Bin/../acc2omp.pl";
chdir($Bin) or die;

sub slurp { my $f = shift; open(my $fh, '<', $f) or die "$f: $!"; local $/; return <$fh>; }
sub nlines { my $s = shift; return ($s =~ tr/\n//); }

my $out = `$^X $tool translate.F90 2>&1`;
is($?, 0, 'translate: exit status 0');
is($out, &slurp('translate.ref'), 'translate: expected output');
is(&nlines($out), &nlines(&slurp('translate.F90')), 'translate: same number of lines');

my $tmp = "translate.out.$$";
`$^X $tool -o $tmp translate.F90 2>&1`;
is($?, 0, 'translate -o: exit status 0');
is(&slurp($tmp), &slurp('translate.ref'), 'translate -o: same output');
unlink($tmp);

$out = `$^X $tool check_subset.F90 2>&1`;
is($?, 0, 'translate check_subset: exit status 0');
is(&nlines($out), &nlines(&slurp('check_subset.F90')), 'translate check_subset: same number of lines');

$out = `$^X $tool -o $tmp check_errors.F90 2>&1`;
isnt($?, 0, 'translate errors: exit status not 0');
is($out, &slurp('check_errors.ref'), 'translate errors: same messages as --check');
ok(! -e $tmp, 'translate errors: no output file');

done_testing();

#!/usr/bin/env perl
# Tests of unpair.pl; run with: prove tools/acc2omp/t
use strict;
use warnings;
use Test::More;
use File::Copy;
use FindBin qw($Bin);

my $tool = "$Bin/../unpair.pl";
chdir($Bin) or die;

sub slurp { my $f = shift; open(my $fh, '<', $f) or die "$f: $!"; local $/; return <$fh>; }

my $out = `$^X $tool --list unpair.F90 2>&1`;
is($?, 0, 'list: exit status 0');
is($out, &slurp('unpair_list.ref'), 'list: blocks and removability');

my @ok = map { m/^unpair\.F90:(\d+): .*: removable$/ ? ($1) : () } split(/\n/, $out);
my $tmp = "unpair.tmp.$$.F90";
copy('unpair.F90', $tmp) or die;
`$^X $tool --dry-run @{[ map { "$tmp:$_" } @ok ]} 2>&1`;
is(&slurp($tmp), &slurp('unpair.F90'), 'dry-run: file unchanged');
`$^X $tool @{[ map { "$tmp:$_" } @ok ]} 2>&1`;
is($?, 0, 'remove: exit status 0');
is(&slurp($tmp), &slurp('unpair.ref'), 'remove: expected file (rules R1 + R2)');
is((&slurp($tmp) =~ tr/\n//), (&slurp('unpair.F90') =~ tr/\n//), 'remove: same number of lines');

copy('unpair.F90', $tmp) or die;
$out = `$^X $tool $tmp:42 $tmp:6 2>&1`;
isnt($?, 0, 'not removable: exit status not 0');
like($out, qr/:42: error: pair not removable: code at line 43/, 'not removable: message');
is(&slurp($tmp), &slurp('unpair.F90'), 'not removable: file unchanged');
unlink($tmp);

# the result translates: the OpenMP directives are on the lines of the hand-written ones
$out = `$^X $Bin/../acc2omp.pl unpair.ref 2>&1`;
my @l = split(/\n/, $out, -1);
like($l[10], qr/^!\$OMP TARGET DATA\b/, 'translated: TARGET DATA on line 11');
like($l[15], qr/^!\$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO/, 'translated: kernel on line 16');

done_testing();
